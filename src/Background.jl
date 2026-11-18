"""
    TabulatedBackground

Background `b(|Q|, ω)` tabulated on `|Q|` shells (edges `qedges`, Å⁻¹) and the energy bins of the histogram it was
estimated from (`ω`, `ωedges`). `value` and `σ` are `nshell × nω`. `count` holds the source bins of each measured
cell and `source` records how each cell was obtained: `:measured`, `:omega` (interpolated along energy), `:q`
(interpolated or extended along `|Q|`), `:fallback`, or `:none`.

Calling `b(q, ω)` returns `(value, σ)` for the shell containing `q` and the energy bin containing `ω`, and
`(NaN, NaN)` outside the table, so it can be passed as `background` to [`cut_path`](@ref), [`to_smag`](@ref) and
[`magnetic_histogram`](@ref). Bins outside the table are then left out.
"""
struct TabulatedBackground
    qedges::Vector{Float64}
    ω::Vector{Float64}
    ωedges::Vector{Float64}
    value::Matrix{Float64}
    σ::Matrix{Float64}
    count::Matrix{Int}
    source::Matrix{Symbol}
end

function _bin_index(edges::AbstractVector{Float64}, x::Real)
    (isfinite(x) || x == last(edges)) && first(edges) <= x <= last(edges) || return 0
    return min(searchsortedlast(edges, x), length(edges) - 1)
end

function (b::TabulatedBackground)(q::Real, ω::Real)
    c, j = _bin_index(b.qedges, q), _bin_index(b.ωedges, ω)
    (c == 0 || j == 0) && return (NaN, NaN)
    return (b.value[c, j], b.σ[c, j])
end

function Base.show(io::IO, b::TabulatedBackground)
    n = length(b.source)
    share(s) = round(100 * count(==(s), b.source) / n; digits=1)
    print(io, "TabulatedBackground(", size(b.value, 1), " shells × ", size(b.value, 2), " energies; measured ",
          share(:measured), "%, along ω ", share(:omega), "%, along |Q| ", share(:q), "%, fallback ", share(:fallback), "%)")
end

"""
    momentum_magnitude(h; reciprocal_basis=nothing) -> Array{Float64}

Return `|Q|` in Å⁻¹ for every momentum bin, from the stored HKL and the histogram lattice (or `reciprocal_basis`,
columns `2π a*`, `2π b*`, `2π c*`). The array has the momentum axes of `h` in order, without the energy axis.
Stored indices give the measured momentum transfer even when the indexing needs correction.
"""
function momentum_magnitude(h::NormalizedHistogram; reciprocal_basis=nothing)
    g = _momentum_geometry(h)
    B = _reciprocal_basis(h, reciprocal_basis)
    return map(i -> norm(B * _stored_hkl(g, i)), CartesianIndices(g.shape))
end

function _site(site)
    length(site) in (2, 3) || throw(ArgumentError("a site is (position, halfwidth) or (position, halfwidth, (ωmin, ωmax))"))
    position, width = site[1], site[2]
    p = _hkl_point(position)
    P = length(p)
    w = width isa Real ? ntuple(_ -> Float64(width), P) : Tuple(Float64.(width))
    length(w) == P && all(x -> isfinite(x) && x >= 0, w) || throw(ArgumentError("invalid exclusion half-width"))
    energy = length(site) == 3 ? (Float64(site[3][1]), Float64(site[3][2])) : (-Inf, Inf)
    energy[1] <= energy[2] || throw(ArgumentError("a site energy window must have ωmin ≤ ωmax"))
    return (p, SVector{P,Float64}(w), energy)
end

function _near_site(q, p::SVector{P}, w::SVector{P}) where {P}
    for d in 1:P
        Δ = q[d] - p[d]
        abs(Δ - round(Δ)) <= w[d] + 8eps(max(1.0, abs(q[d]))) || return false
    end
    return true
end

function _qband(band)
    centre, width = band
    isfinite(centre) && isfinite(width) && width >= 0 || throw(ArgumentError("invalid |Q| exclusion band"))
    return (Float64(centre), Float64(width))
end

_mode_energies(x::Real) = (Float64(x),)
_mode_energies(x) = x

"""
    exclusion_mask(h; sites=(), dispersion=nothing, window=1.0, qmin=0, qmax=Inf, qbands=(),
                   indices=nothing, reciprocal_basis=nothing) -> BitArray

Mark histogram bins to leave out of a background or a scale estimate (`true` = excluded).

- `sites`: pairs `(position, halfwidth)` in HKL, repeated over all reciprocal-lattice translations. A bin is
  excluded when every component lies within `halfwidth` of the site; 2D positions test `(H, K)` only. Use them for
  Bragg peaks, nuclear zone centres and magnetic zone centres. A third element `(ωmin, ωmax)` limits the site to
  that energy range, e.g. `((0.5, 0.5), 0.075, (-Inf, 0.5))` removes the elastic tail of a magnetic Bragg peak
  while keeping the inelastic signal at the same wave vector.
- `dispersion`: a function of the bin's `(H, K, L)` returning one mode energy or a collection of them; bins
  with `|ω - ω_mode| ≤ window` are excluded.
- `qmin`, `qmax`, `qbands`: limits on `|Q|` (Å⁻¹) and bands `(centre, halfwidth)` around powder lines.

Positions and the dispersion use `indices` (e.g. [`corrected_indices`](@ref)) when given, and the stored HKL
otherwise. `|Q|` always uses the stored HKL. Combine masks with `.|`.
"""
function exclusion_mask(h::NormalizedHistogram; sites=(), dispersion=nothing, window::Real=1.0, qmin::Real=0.0,
                        qmax::Real=Inf, qbands=(), indices=nothing, reciprocal_basis=nothing)
    isfinite(window) && window >= 0 || throw(ArgumentError("window must be finite and nonnegative"))
    qmin <= qmax || throw(ArgumentError("qmin must not exceed qmax"))
    g = _momentum_geometry(h)
    if indices !== nothing
        axes(indices) == map(Base.OneTo, g.shape) || throw(DimensionMismatch("indices must match the momentum grid"))
    end
    parsed = [_site(s) for s in sites]
    bands = [_qband(b) for b in qbands]
    needs_q = qmin > 0 || isfinite(qmax) || !isempty(bands)
    Q = needs_q ? momentum_magnitude(h; reciprocal_basis) : nothing
    momentum = falses(g.shape)
    windows = [NTuple{2,Float64}[] for _ in 1:prod(g.shape)]
    modes = Vector{Vector{Float64}}(undef, prod(g.shape))
    for (k, i) in enumerate(CartesianIndices(g.shape))
        hkl = indices === nothing ? Tuple(_stored_hkl(g, i)) : Tuple(Float64.(indices[i]))
        excluded = false
        for (p, w, e) in parsed
            _near_site(hkl, p, w) || continue
            isinf(e[1]) && isinf(e[2]) && e[1] < 0 < e[2] ? (excluded = true) : push!(windows[k], e)
        end
        if Q !== nothing
            q = Q[i]
            excluded |= q < qmin || q > qmax || any(((c, w),) -> abs(q - c) <= w, bands)
        end
        momentum[i] = excluded
        modes[k] = dispersion === nothing ? Float64[] : collect(Float64, _mode_energies(dispersion(hkl)))
    end
    ω = h.centres[g.energy]
    ids = LinearIndices(g.shape)
    mask = falses(size(h.D))
    for i in CartesianIndices(h.D)
        m = _momentum_index(g, i)
        e = ω[i[g.energy]]
        mask[i] = momentum[m] || any(((a, b),) -> a <= e <= b, windows[ids[m]]) ||
                  any(x -> abs(e - x) <= window, modes[ids[m]])
    end
    return mask
end

function _check_exclude(exclude, h)
    exclude === nothing && return nothing
    exclude isa AbstractArray{Bool} && size(exclude) == size(h.D) ||
        throw(DimensionMismatch("exclude must be a Boolean array with the histogram shape"))
    return exclude
end

_excluded(::Nothing, i) = false
_excluded(mask, i) = mask[i]

# Sums of D, D2 and N, and source-bin counts, per (shell, energy).
function _accumulate_cells!(sums, h, exclude, shell::F) where {F}
    D, D2, N = h.D, h.D2, h.N
    e = energy_axis(h)
    for i in CartesianIndices(D)
        n = N[i]
        n > 0 && !_excluded(exclude, i) || continue
        c = shell(i)
        c > 0 || continue
        j = i[e]
        sums.D[c, j] += D[i]
        sums.N[c, j] += n
        D2 === nothing || (sums.V[c, j] += D2[i])
        sums.count[c, j] += 1
    end
    return sums
end

_cell_sums(nshell, nω) = (; D=zeros(nshell, nω), N=zeros(nshell, nω), V=zeros(nshell, nω), count=zeros(Int, nshell, nω))

function _measured_cells(sums, min_bins, has_variance)
    value, σ = fill(NaN, size(sums.N)), fill(NaN, size(sums.N))
    source = fill(:none, size(sums.N))
    for k in eachindex(value)
        sums.count[k] >= min_bins && sums.N[k] > 0 || continue
        value[k] = sums.D[k] / sums.N[k]
        σ[k] = has_variance ? sqrt(sums.V[k]) / sums.N[k] : 0.0
        source[k] = :measured
    end
    return value, σ, source
end

function _fill_nearest!(v, σ, source, tag)
    ok = findall(isfinite, v)
    isempty(ok) && return v
    for j in eachindex(v)
        isfinite(v[j]) && continue
        k = ok[argmin(abs.(ok .- j))]
        v[j], σ[j], source[j] = v[k], σ[k], tag
    end
    return v
end

"""
    background_profile(h; exclude=nothing, fill=false, min_bins=1) -> TabulatedBackground

Momentum-independent background `b(ω) = ΣD/ΣN` over the measured bins not in `exclude`, with
`σ = sqrt(ΣD2)/ΣN`. The result has a single `|Q|` shell covering all momenta. With `fill = true`, energies with
fewer than `min_bins` source bins take the nearest measured energy.

    background_profile(path_D, path_N; dim, exclude=nothing, ranges=nothing, fill=false, min_bins=1)

Accumulate the same sums one slice of axis `dim` at a time (see [`eachslice_mdnorm`](@ref)). Here `exclude` is
`nothing` or a function of the slice returning its mask.
"""
function background_profile(h::NormalizedHistogram; exclude=nothing, fill::Bool=false, min_bins::Integer=1)
    mask = _check_exclude(exclude, h)
    e = energy_axis(h)
    sums = _accumulate_cells!(_cell_sums(1, size(h.D, e)), h, mask, _ -> 1)
    return _profile_table(sums, h.centres[e], h.edges[e], h.D2 !== nothing, fill, min_bins)
end

function background_profile(path_D::AbstractString, path_N::AbstractString; dim::Integer, exclude=nothing,
                            ranges=nothing, fill::Bool=false, min_bins::Integer=1)
    state = Ref{Any}(nothing)
    eachslice_mdnorm(path_D, path_N; dim, ranges) do slice, _
        e = energy_axis(slice)
        if state[] === nothing
            state[] = (; sums=_cell_sums(1, size(slice.D, e)), ω=slice.centres[e], ωedges=slice.edges[e],
                       variance=slice.D2 !== nothing)
        end
        mask = exclude === nothing ? nothing : _check_exclude(exclude(slice), slice)
        _accumulate_cells!(state[].sums, slice, mask, _ -> 1)
    end
    s = state[]
    return _profile_table(s.sums, s.ω, s.ωedges, s.variance, fill, min_bins)
end

function _profile_table(sums, ω, ωedges, variance, fill, min_bins)
    min_bins >= 1 || throw(ArgumentError("min_bins must be positive"))
    value, σ, source = _measured_cells(sums, min_bins, variance)
    fill && _fill_nearest!(view(value, 1, :), view(σ, 1, :), view(source, 1, :), :omega)
    return TabulatedBackground([0.0, Inf], collect(Float64, ω), collect(Float64, ωedges), value, σ, sums.count, source)
end

function _fill_interior!(v, σ, source)
    ok = findall(isfinite, v)
    for a in 1:length(ok)-1
        l, r = ok[a], ok[a+1]
        for x in l+1:r-1
            w = (x - l) / (r - l)
            v[x] = (1 - w) * v[l] + w * v[r]
            σ[x] = max(σ[l], σ[r])
            source[x] = :omega
        end
    end
    return v
end

function _fill_shells!(v, σ, source)
    ok = findall(isfinite, v)
    isempty(ok) && return false
    for c in eachindex(v)
        isfinite(v[c]) && continue
        il, ih = findlast(<(c), ok), findfirst(>(c), ok)
        if il === nothing || ih === nothing
            k = ok[il === nothing ? ih : il]
            v[c], σ[c] = v[k], σ[k]
        else
            l, r = ok[il], ok[ih]
            w = (c - l) / (r - l)
            v[c] = (1 - w) * v[l] + w * v[r]
            σ[c] = max(σ[l], σ[r])
        end
        source[c] = :q
    end
    return true
end

function _fallback_column(fallback::TabulatedBackground, ω)
    size(fallback.value, 1) == 1 || throw(ArgumentError("fallback must be a single-shell background profile"))
    return fallback(0.0, ω)
end
_fallback_column(fallback::Real, ω) = (Float64(fallback), 0.0)
_fallback_column(::Nothing, ω) = (NaN, NaN)

"""
    background_shells(h; exclude=nothing, dq=0.15, qmax=nothing, min_bins=40, interpolate=:omega_then_q,
                      smooth=0.0, smooth_from=-Inf, fallback=nothing, reciprocal_basis=nothing) -> TabulatedBackground

Background `b(|Q|, ω) = ΣD/ΣN` in `|Q|` shells of width `dq` (Å⁻¹) from 0 to `qmax`, over the measured bins not in
`exclude`. A (shell, energy) cell counts as measured when it holds at least `min_bins` source bins. Missing cells
are filled in order:

1. `interpolate = :omega_then_q` interpolates linearly along energy inside each shell, between measured cells
   only (for example across an excluded dispersion); `:q` skips this step and `:none` skips both.
2. Along `|Q|`, linearly between measured shells, and from the nearest measured shell at the ends.
3. Energies with no measured shell take `fallback` (a [`background_profile`](@ref), a number, or `nothing`).

Filled cells take the larger neighbouring σ. Finally, if `smooth > 0`, a boxcar of width `smooth` (meV) is applied
along energy to all energies at or above `smooth_from`. `qmax` defaults to the largest measured `|Q|`.
"""
function background_shells(h::NormalizedHistogram; exclude=nothing, dq::Real=0.15, qmax=nothing,
                           min_bins::Integer=40, interpolate::Symbol=:omega_then_q, smooth::Real=0.0,
                           smooth_from::Real=-Inf, fallback=nothing, reciprocal_basis=nothing)
    interpolate in (:omega_then_q, :q, :none) || throw(ArgumentError("interpolate must be :omega_then_q, :q or :none"))
    isfinite(dq) && dq > 0 || throw(ArgumentError("dq must be positive"))
    min_bins >= 1 || throw(ArgumentError("min_bins must be positive"))
    isfinite(smooth) && smooth >= 0 || throw(ArgumentError("smooth must be finite and nonnegative"))
    mask = _check_exclude(exclude, h)
    g = _momentum_geometry(h)
    Q = momentum_magnitude(h; reciprocal_basis)
    top = qmax === nothing ? maximum((Q[_momentum_index(g, i)] for i in CartesianIndices(h.D) if h.N[i] > 0); init=0.0) :
          Float64(qmax)
    isfinite(top) && top > 0 || throw(ArgumentError("qmax must be positive"))
    nshell = max(1, ceil(Int, top / dq - 8eps(top / dq)))
    qedges = [collect((0:nshell-1) .* Float64(dq)); top]
    shell = i -> (q = Q[_momentum_index(g, i)]; q <= top ? min(floor(Int, q / dq) + 1, nshell) : 0)
    e = g.energy
    ω, ωedges = h.centres[e], h.edges[e]
    sums = _accumulate_cells!(_cell_sums(nshell, length(ω)), h, mask, shell)
    value, σ, source = _measured_cells(sums, min_bins, h.D2 !== nothing)
    if interpolate === :omega_then_q
        for c in 1:nshell
            _fill_interior!(view(value, c, :), view(σ, c, :), view(source, c, :))
        end
    end
    for j in eachindex(ω)
        if interpolate === :none || !_fill_shells!(view(value, :, j), view(σ, :, j), view(source, :, j))
            any(isfinite, view(value, :, j)) && continue
            b, s = _fallback_column(fallback, ω[j])
            isfinite(b) || continue
            value[:, j] .= b
            σ[:, j] .= s
            source[:, j] .= :fallback
        end
    end
    smooth > 0 && _smooth_energy!(value, ω, smooth, smooth_from)
    return TabulatedBackground(qedges, collect(Float64, ω), collect(Float64, ωedges), value, σ, sums.count, source)
end

function _smooth_energy!(value, ω, width, from)
    js = findall(>=(from), ω)
    length(js) > 1 || return value
    hw = max(0, round(Int, width / (2 * mean(diff(ω)))))
    hw == 0 && return value
    for c in axes(value, 1)
        row = value[c, js]
        for (a, j) in enumerate(js)
            window = filter(isfinite, view(row, max(1, a - hw):min(length(js), a + hw)))
            isempty(window) || (value[c, j] = mean(window))
        end
    end
    return value
end

"""
    magnetic_histogram(h, g; background=nothing, exclude=nothing, reciprocal_basis=nothing) -> NormalizedHistogram

Return a histogram whose ratio is `S_mag`: numerator `D - N b`, normalization `N g`, and variance `D2 + (N σ_b)²`.
`g` and `background` take the forms accepted by [`to_smag`](@ref). Bins in `exclude`, or with a nonfinite
conversion or background, become unmeasured. [`rebin`](@ref), [`integrate`](@ref), [`fold_to_zone`](@ref) and
[`cut_path`](@ref) then combine bins as `Σ(D - N b)/Σ(N g)`. Weights derived from `N` (for example
[`weighted_average`](@ref)) become weights in `N g`.
"""
function magnetic_histogram(h::NormalizedHistogram, g; background=nothing, exclude=nothing, reciprocal_basis=nothing)
    _check_bin_array(g, h, "g")
    _check_bin_array(background, h, "background")
    mask = _check_exclude(exclude, h)
    geom = _uses_background_geometry(background) ? _momentum_geometry(h) : nothing
    norms = _background_norms(h, geom, background, reciprocal_basis)
    coords = geom === nothing ? i -> (0.0, 0.0) :
        i -> (norms[_momentum_index(geom, i)], h.centres[geom.energy][i[geom.energy]])
    D, D2, N = _magnetic_arrays(h.D, h.D2, h.N, coords, g, background, mask)
    return NormalizedHistogram(D, D2, N, h.edges, h.centres, h.names, h.units, h.frames, h.meta)
end

function _magnetic_arrays(D0, D20::V, N0, coords::F, g::G, background::B, mask) where {V,F,G,B}
    D, N = zeros(size(D0)), zeros(size(D0))
    D2 = D20 === nothing ? nothing : zeros(size(D0))
    for i in CartesianIndices(D0)
        _excluded(mask, i) && continue
        q, ω = coords(i)
        bin = _conversion_bin(D0, N0, g, background, i, q, ω)
        bin === nothing && continue
        numerator, weight, var, _ = bin
        weight > 0 || continue
        D[i], N[i] = numerator, weight
        D2 === nothing || (D2[i] = D20[i] + var)
    end
    return D, D2, N
end

_gauss_shape(x, fwhm) = exp(-4log(2) * (x / fwhm)^2)

# Width and centre of a Gaussian elastic line from the row-summed profile: half-maximum estimate, then a grid.
function _elastic_gaussian(ω, profile)
    k = argmax(profile)
    base = minimum(profile)
    half = base + (profile[k] - base) / 2
    lo, hi = k, k
    while lo > 1 && profile[lo-1] >= half; lo -= 1; end
    while hi < length(ω) && profile[hi+1] >= half; hi += 1; end
    Δ = abs(ω[min(k + 1, end)] - ω[max(k - 1, 1)]) / 2
    Γ0 = max(abs(ω[hi] - ω[lo]) + Δ, 2Δ)
    best = (Inf, Γ0, Float64(ω[k]))
    for Γ in range(0.5Γ0, 1.5Γ0; length=101), c in range(ω[k] - Γ0 / 2, ω[k] + Γ0 / 2; length=51)
        use = [j for j in eachindex(ω) if abs(ω[j] - c) <= 1.5Γ]
        length(use) >= 4 || continue
        X = [_gauss_shape.(ω[use] .- c, Γ) ones(length(use))]
        r = sum(abs2, profile[use] .- X * (X \ profile[use]))
        r < best[1] && (best = (r, Γ, c))
    end
    return best[2], best[3]
end

"""
    subtract_elastic(I, ω; σ=nothing, valid=nothing, lineshape=nothing, fwhm=nothing, centre=nothing,
                     window=nothing, min_bins=4)
        -> (; I, σ, model, amplitude, σ_amplitude, fwhm, centre)

Remove the elastic line from every momentum row of `I` (rows: momentum, columns: energy `ω`). Each row is fitted
over `|ω - centre| ≤ window` with `A f(ω - centre) + c`, and `A f` is subtracted at all energies, so its tail no
longer reaches the inelastic window. This keeps the inelastic signal under the tail, which an energy-limited site
exclusion in [`exclusion_mask`](@ref) would discard.

`lineshape` is the instrument's elastic line `f(Δω)` (any normalization); give it whenever it is known, because
a real line usually has heavier tails than a Gaussian and extrapolating a Gaussian under-subtracts exactly where
the tail matters. Without it a Gaussian is used, with `fwhm` and `centre` fitted to the row-summed profile when
not given. `window` defaults to the FWHM (required with a `lineshape`). With `σ`, rows are fitted with weights
`1/σ²` and the amplitude error is added in quadrature to the returned `σ`. Unmeasured bins (`valid = false`) are
neither fitted nor changed.
"""
function subtract_elastic(I::AbstractMatrix{<:Real}, ω::AbstractVector{<:Real}; σ=nothing, valid=nothing,
                          lineshape=nothing, fwhm=nothing, centre=nothing, window=nothing, min_bins::Integer=4)
    size(I, 2) == length(ω) || throw(DimensionMismatch("I must have one column per energy"))
    σ === nothing || size(σ) == size(I) || throw(DimensionMismatch("σ must match I"))
    measured = valid === nothing ? isfinite.(I) : BitMatrix(valid) .& isfinite.(I)
    size(measured) == size(I) || throw(DimensionMismatch("valid must match I"))
    ωf = Float64.(ω)
    if lineshape === nothing
        if fwhm === nothing || centre === nothing
            profile = [sum(I[i, j] for i in axes(I, 1) if measured[i, j]; init=0.0) for j in eachindex(ωf)]
            Γ, c = _elastic_gaussian(ωf, profile)
            fwhm = fwhm === nothing ? Γ : fwhm
            centre = centre === nothing ? c : centre
        end
        Γ = Float64(fwhm)
        f = x -> _gauss_shape(x, Γ)
        window = window === nothing ? Γ : window
    else
        window === nothing && throw(ArgumentError("give the fit window with a lineshape"))
        f = lineshape
        centre = centre === nothing ? 0.0 : centre
    end
    c0 = Float64(centre); W = Float64(window)
    W > 0 || throw(ArgumentError("window must be positive"))
    shape = [Float64(f(x - c0)) for x in ωf]
    fitbins = findall(x -> abs(x - c0) <= W, ωf)
    out = Matrix{Float64}(I); model = zeros(size(I))
    σout = σ === nothing ? nothing : Matrix{Float64}(σ)
    amplitude = zeros(size(I, 1)); σamp = fill(NaN, size(I, 1))
    for i in axes(I, 1)
        use = [j for j in fitbins if measured[i, j] && (σ === nothing || σ[i, j] > 0)]
        length(use) >= min_bins || continue
        w = σ === nothing ? ones(length(use)) : 1 ./ Float64.(σ[i, use])
        X = [shape[use] ones(length(use))] .* w
        y = Float64.(I[i, use]) .* w
        θ = X \ y
        C = inv(X' * X)
        if σ === nothing
            C .*= sum(abs2, y .- X * θ) / max(length(use) - 2, 1)
        end
        amplitude[i] = θ[1]; σamp[i] = sqrt(max(C[1, 1], 0.0))
        for j in eachindex(ωf)
            measured[i, j] || continue
            model[i, j] = θ[1] * shape[j]
            out[i, j] -= model[i, j]
            σout === nothing || (σout[i, j] = hypot(σout[i, j], σamp[i] * shape[j]))
        end
    end
    return (; I=out, σ=σout, model, amplitude, σ_amplitude=σamp, fwhm=lineshape === nothing ? Float64(fwhm) : NaN,
            centre=c0)
end
