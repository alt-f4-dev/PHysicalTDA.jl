"""
    GAMMA_R0_SQ

`(γ r₀)² = 0.2905` barn, the square of the neutron magnetic scattering length.
"""
const GAMMA_R0_SQ = 0.2905

# Sears, Neutron News 3, 29 (1992), via NIST: b_coh (fm); σ_coh, σ_inc, σ_scatt and σ_abs (barn, 2200 m/s).
# Natural abundance; `:D` is ²H. Elements with complex scattering lengths are not tabulated.
const _NEUTRON_XS = Dict{Symbol,NTuple{5,Float64}}(
    :H => (-3.7390, 1.7568, 80.26, 82.02, 0.3326), :D => (6.671, 5.592, 2.05, 7.64, 0.000519),
    :Li => (-1.90, 0.454, 0.92, 1.37, 70.5), :C => (6.6460, 5.551, 0.001, 5.551, 0.0035),
    :N => (9.36, 11.01, 0.5, 11.51, 1.9), :O => (5.803, 4.232, 0.0008, 4.232, 0.00019),
    :F => (5.654, 4.017, 0.0008, 4.018, 0.0096), :Na => (3.63, 1.66, 1.62, 3.28, 0.53),
    :Mg => (5.375, 3.631, 0.08, 3.71, 0.063), :Al => (3.449, 1.495, 0.0082, 1.503, 0.231),
    :Si => (4.1491, 2.163, 0.004, 2.167, 0.171), :P => (5.13, 3.307, 0.005, 3.312, 0.172),
    :S => (2.847, 1.0186, 0.007, 1.026, 0.53), :Cl => (9.5770, 11.5257, 5.3, 16.8, 33.5),
    :K => (3.67, 1.69, 0.27, 1.96, 2.1), :Ca => (4.70, 2.78, 0.05, 2.83, 0.43),
    :Ti => (-3.438, 1.485, 2.87, 4.35, 6.09), :V => (-0.3824, 0.0184, 5.08, 5.1, 5.08),
    :Cr => (3.635, 1.66, 1.83, 3.49, 3.05), :Mn => (-3.73, 1.75, 0.4, 2.15, 13.3),
    :Fe => (9.45, 11.22, 0.4, 11.62, 2.56), :Co => (2.49, 0.779, 4.8, 5.6, 37.18),
    :Ni => (10.3, 13.3, 5.2, 18.5, 4.49), :Cu => (7.718, 7.485, 0.55, 8.03, 3.78),
    :Zn => (5.680, 4.054, 0.077, 4.131, 1.11), :Ga => (7.288, 6.675, 0.16, 6.83, 2.75),
    :Ge => (8.185, 8.42, 0.18, 8.6, 2.2), :As => (6.58, 5.44, 0.06, 5.5, 4.5),
    :Se => (7.970, 7.98, 0.32, 8.3, 11.7), :Br => (6.795, 5.8, 0.1, 5.9, 6.9),
    :Rb => (7.09, 6.32, 0.5, 6.8, 0.38), :Sr => (7.02, 6.19, 0.06, 6.25, 1.28),
    :Y => (7.75, 7.55, 0.15, 7.7, 1.28), :Zr => (7.16, 6.44, 0.02, 6.46, 0.185),
    :Nb => (7.054, 6.253, 0.0024, 6.255, 1.15), :Mo => (6.715, 5.67, 0.04, 5.71, 2.48),
    :Ru => (7.03, 6.21, 0.4, 6.6, 2.56), :Rh => (5.88, 4.34, 0.3, 4.6, 144.8),
    :Pd => (5.91, 4.39, 0.093, 4.48, 6.9), :Ag => (5.922, 4.407, 0.58, 4.99, 63.3),
    :Sn => (6.225, 4.871, 0.022, 4.892, 0.626), :Sb => (5.57, 3.9, 0.007, 3.9, 4.91),
    :Te => (5.80, 4.23, 0.09, 4.32, 4.7), :I => (5.28, 3.5, 0.31, 3.81, 6.15),
    :Cs => (5.42, 3.69, 0.21, 3.9, 29.0), :Ba => (5.07, 3.23, 0.15, 3.38, 1.1),
    :La => (8.24, 8.53, 1.13, 9.66, 8.97), :Ce => (4.84, 2.94, 0.001, 2.94, 0.63),
    :Pr => (4.58, 2.64, 0.015, 2.66, 11.5), :Nd => (7.69, 7.43, 9.2, 16.6, 50.5),
    :Tb => (7.38, 6.84, 0.004, 6.84, 23.4), :Ho => (8.01, 8.06, 0.36, 8.42, 64.7),
    :Er => (7.79, 7.63, 1.1, 8.7, 159.0), :Yb => (12.43, 19.42, 4.0, 23.4, 34.8),
    :Ta => (6.91, 6.0, 0.01, 6.01, 20.6), :W => (4.86, 2.97, 1.63, 4.6, 18.3),
    :Re => (9.2, 10.6, 0.9, 11.5, 89.7), :Os => (10.7, 14.4, 0.3, 14.7, 16.0),
    :Pt => (9.60, 11.58, 0.13, 11.71, 10.3), :Au => (7.63, 7.32, 0.43, 7.75, 98.65),
    :Pb => (9.405, 11.115, 0.003, 11.118, 0.171), :Bi => (8.532, 9.148, 0.0084, 9.156, 0.0338),
    :U => (8.417, 8.903, 0.005, 8.908, 7.57),
)

"""
    neutron_xs(element) -> (b_coh, σ_coh, σ_inc, σ_scatt, σ_abs)

Natural-abundance scattering length `b_coh` (fm) and cross sections (barn; `σ_abs` at 2200 m/s) from Sears,
Neutron News 3, 29 (1992). `element` is a symbol or string such as `:Mn` or `"Rb"`; `:D` selects ²H.
"""
function neutron_xs(element::Union{Symbol,AbstractString})
    key = Symbol(element)
    haskey(_NEUTRON_XS, key) || throw(ArgumentError("no tabulated neutron cross sections for $element"))
    b, coh, inc, scatt, abs = _NEUTRON_XS[key]
    return (; b_coh=b, σ_coh=coh, σ_inc=inc, σ_scatt=scatt, σ_abs=abs)
end

function _parse_formula(formula::AbstractString)
    stack = [Dict{Symbol,Float64}()]
    i = firstindex(formula)
    number(j) = (m = match(r"^\d*\.?\d+|^\d+", SubString(formula, j)); m === nothing ? (1.0, j) :
                 (parse(Float64, m.match), j + ncodeunits(m.match)))
    while i <= lastindex(formula)
        c = formula[i]
        if isspace(c)
            i = nextind(formula, i)
        elseif c == '('
            push!(stack, Dict{Symbol,Float64}())
            i = nextind(formula, i)
        elseif c == ')'
            length(stack) > 1 || throw(ArgumentError("unbalanced parentheses in $formula"))
            group = pop!(stack)
            n, i = number(nextind(formula, i))
            for (el, x) in group
                stack[end][el] = get(stack[end], el, 0.0) + n * x
            end
        elseif isuppercase(c)
            m = match(r"^[A-Z][a-z]?", SubString(formula, i))
            el = Symbol(m.match)
            n, i = number(i + ncodeunits(m.match))
            stack[end][el] = get(stack[end], el, 0.0) + n
        else
            throw(ArgumentError("cannot parse chemical formula $formula"))
        end
    end
    length(stack) == 1 || throw(ArgumentError("unbalanced parentheses in $formula"))
    isempty(stack[1]) && throw(ArgumentError("empty chemical formula"))
    return stack[1]
end

"""
    formula_xs(formula) -> (σ_coh, σ_inc, σ_scatt, σ_abs, composition)

Cross sections (barn) per formula unit, summed over the atoms of `formula` (for example `"Rb2MnF4"` or
`"Cu(NO3)2"`) with [`neutron_xs`](@ref). `composition` maps each element to its count.
"""
function formula_xs(formula::AbstractString)
    composition = _parse_formula(formula)
    totals = zeros(4)
    for (el, n) in composition
        xs = neutron_xs(el)
        totals .+= n .* (xs.σ_coh, xs.σ_inc, xs.σ_scatt, xs.σ_abs)
    end
    return (; σ_coh=totals[1], σ_inc=totals[2], σ_scatt=totals[3], σ_abs=totals[4], composition)
end

"""
    elastic_fwhm(ω, y; window=(-3, 3), centres=-0.3:0.01:0.3, widths=exp.(range(log(0.03), log(1.5); length=120)))
    elastic_fwhm(b::TabulatedBackground; kwargs...)
    elastic_fwhm(h::NormalizedHistogram; exclude=nothing, kwargs...)

Full width at half maximum (meV) of the elastic line: least-squares fit of a constant plus a Gaussian over the
grid of `centres` and standard deviations `widths`, using the finite points with `ω` in `window`. The histogram
method fits its [`background_profile`](@ref).
"""
function elastic_fwhm(ω::AbstractVector{<:Real}, y::AbstractVector{<:Real}; window=(-3.0, 3.0),
                      centres=-0.3:0.01:0.3, widths=exp.(range(log(0.03), log(1.5); length=120)))
    length(ω) == length(y) || throw(DimensionMismatch("ω and y must have the same length"))
    keep = [j for j in eachindex(ω) if window[1] <= ω[j] <= window[2] && isfinite(y[j])]
    length(keep) >= 3 || throw(ArgumentError("too few finite points in the elastic window"))
    x, v = Float64.(ω[keep]), Float64.(y[keep])
    best, width = Inf, NaN
    X = ones(length(x), 2)
    for μ in centres, σ in widths
        @. X[:, 2] = exp(-(x - μ)^2 / (2σ^2))
        r = sum(abs2, v .- X * (X \ v))
        r < best && ((best, width) = (r, σ))
    end
    return 2 * sqrt(2 * log(2)) * width
end

function elastic_fwhm(b::TabulatedBackground; kwargs...)
    size(b.value, 1) == 1 || throw(ArgumentError("expected a single-shell background profile"))
    return elastic_fwhm(b.ω, vec(b.value); kwargs...)
end

elastic_fwhm(h::NormalizedHistogram; exclude=nothing, kwargs...) = elastic_fwhm(background_profile(h; exclude); kwargs...)

"""
    incoherent_scale(profile::TabulatedBackground; formula=nothing, σ_inc=nothing, magnetic_ions=1, fwhm=nothing,
                     min_window=(-1, 1), window_fwhm=(4, 2.5), flank_width=nothing)
    incoherent_scale(h::NormalizedHistogram; exclude=nothing, kwargs...)

Absolute scale from the incoherent elastic line. The integral `I_el = Σ (b(ω) - floor) Δω` runs over the window
`(-max(-min_window[1], window_fwhm[1] fwhm), max(min_window[2], window_fwhm[2] fwhm))`. The floor is the lower of
the mean levels in the two flanks of width `flank_width` (default `max(1, fwhm)` meV) outside the window; the
lower floor gives the larger integral and so the smaller scale. With `σ_inc` per formula unit (barn), given
directly or through `formula` ([`formula_xs`](@ref)),

    κ = σ_inc / (4π I_el (γ r₀)² magnetic_ions)

converts the ratio `(D/N)/g` into `S_mag` per magnetic ion, where `magnetic_ions` is their number per formula unit.
`fwhm` defaults to [`elastic_fwhm`](@ref) of the profile. Elastic scattering other than the sample's incoherent
scattering (container, mount) makes `κ` too small; see [`container_corrected_scale`](@ref).

Returns `(κ, integral, floor, window, fwhm, σ_inc, below_fwhm)`, where `below_fwhm` is the share of the integral
below `-fwhm`.
"""
function incoherent_scale(b::TabulatedBackground; formula=nothing, σ_inc=nothing, magnetic_ions::Real=1, fwhm=nothing,
                          min_window=(-1.0, 1.0), window_fwhm=(4.0, 2.5), flank_width=nothing)
    size(b.value, 1) == 1 || throw(ArgumentError("expected a single-shell background profile"))
    (formula === nothing) != (σ_inc === nothing) || throw(ArgumentError("give exactly one of formula and σ_inc"))
    xs = σ_inc === nothing ? formula_xs(formula).σ_inc : Float64(σ_inc)
    xs > 0 && magnetic_ions > 0 || throw(ArgumentError("σ_inc and magnetic_ions must be positive"))
    width = fwhm === nothing ? elastic_fwhm(b) : Float64(fwhm)
    isfinite(width) && width > 0 || throw(ArgumentError("fwhm must be positive"))
    window = (-max(-min_window[1], window_fwhm[1] * width), max(min_window[2], window_fwhm[2] * width))
    fw = flank_width === nothing ? max(1.0, width) : Float64(flank_width)
    ω, y, Δω = b.ω, vec(b.value), diff(b.ωedges)
    function flank(lo, hi)
        v = [y[j] for j in eachindex(ω) if lo < ω[j] <= hi && isfinite(y[j])]
        return isempty(v) ? Inf : mean(v)
    end
    level = min(flank(window[2], window[2] + fw), flank(window[1] - fw, window[1]))
    isfinite(level) || throw(ArgumentError("no measured flank outside the elastic window"))
    inside = [j for j in eachindex(ω) if window[1] <= ω[j] <= window[2] && isfinite(y[j])]
    integral = sum((y[j] - level) * Δω[j] for j in inside; init=0.0)
    integral > 0 || throw(ArgumentError("the elastic line does not rise above its floor"))
    below = sum((y[j] - level) * Δω[j] for j in inside if ω[j] < -width; init=0.0) / integral
    κ = xs / (4π * integral * GAMMA_R0_SQ * magnetic_ions)
    return (; κ, integral, floor=level, window, fwhm=width, σ_inc=xs, below_fwhm=below)
end

incoherent_scale(h::NormalizedHistogram; exclude=nothing, kwargs...) =
    incoherent_scale(background_profile(h; exclude); kwargs...)

"""
    container_corrected_scale(κ, share, σ_share; nσ=2, max_share=0.95) -> (κ, κ_lo)

Remove a share `share ± σ_share` of the flat elastic line that belongs to another scatterer (for example
aluminium, from [`powder_line_flat`](@ref)) from an incoherent scale: `κ/(1 - share)`, with the low end
`κ/(1 - max(share - nσ σ_share, 0))`. Returns `NaN` for both when `share ≥ max_share`.
"""
function container_corrected_scale(κ::Real, share::Real, σ_share::Real; nσ::Real=2, max_share::Real=0.95)
    σ_share >= 0 && nσ >= 0 || throw(ArgumentError("σ_share and nσ must be nonnegative"))
    isfinite(share) && share < max_share || return (; κ=NaN, κ_lo=NaN)
    return (; κ=κ / (1 - share), κ_lo=κ / (1 - max(share - nσ * σ_share, 0.0)))
end

# Radial profile accumulators: Σ D and Σ N per (|Q| bin, energy bin) over the energy window.
function _accumulate_radial!(acc, h, mask, qedges, energy_window, reciprocal_basis)
    g = _momentum_geometry(h)
    Q = momentum_magnitude(h; reciprocal_basis)
    e = g.energy
    ω = h.centres[e]
    nq = length(qedges) - 1
    for i in CartesianIndices(h.D)
        n = h.N[i]
        n > 0 && !_excluded(mask, i) || continue
        j = i[e]
        energy_window[1] <= ω[j] <= energy_window[2] || continue
        c = searchsortedlast(qedges, Q[_momentum_index(g, i)])
        1 <= c <= nq || continue
        acc.D[c, j] += h.D[i]
        acc.N[c, j] += n
    end
    return acc
end

function _radial_result(acc, qedges, ω, ωedges, energy_window)
    js = [j for j in eachindex(ω) if energy_window[1] <= ω[j] <= energy_window[2]]
    isempty(js) && throw(ArgumentError("no energy bins in the elastic window"))
    Δω = diff(ωedges)
    R = map(axes(acc.N, 1)) do c
        all(j -> acc.N[c, j] > 0, js) || return NaN
        return sum(acc.D[c, j] / acc.N[c, j] * Δω[j] for j in js)
    end
    Q = [(qedges[c] + qedges[c+1]) / 2 for c in 1:length(qedges)-1]
    return (; Q, R, qedges)
end

"""
    radial_elastic_profile(h; energy_window, dq=0.02, qmax=5.0, exclude=nothing, reciprocal_basis=nothing)
    radial_elastic_profile(path_D, path_N; dim, energy_window, ranges=nothing, exclude=nothing, kwargs...)

Energy-integrated elastic intensity against `|Q|`: `R(Q) = Σ_ω (ΣD/ΣN) Δω` over the bins in `|Q|` bins of width
`dq`, for the energies in `energy_window`. `|Q|` includes every momentum axis, so an `(H, K, E, L)` workspace gives
the three-dimensional `|Q|` needed for powder lines. A `|Q|` bin missing any energy of the window is `NaN`. The
file method accumulates one slice of axis `dim` at a time; its `exclude` is a function of the slice.

Returns `(Q, R, qedges)`, ready for [`fit_powder_lines`](@ref).
"""
function radial_elastic_profile(h::NormalizedHistogram; energy_window, dq::Real=0.02, qmax::Real=5.0, exclude=nothing,
                                reciprocal_basis=nothing)
    qedges = _radial_edges(dq, qmax)
    e = energy_axis(h)
    acc = (; D=zeros(length(qedges) - 1, size(h.D, e)), N=zeros(length(qedges) - 1, size(h.D, e)))
    _accumulate_radial!(acc, h, _check_exclude(exclude, h), qedges, energy_window, reciprocal_basis)
    return _radial_result(acc, qedges, h.centres[e], h.edges[e], energy_window)
end

function radial_elastic_profile(path_D::AbstractString, path_N::AbstractString; dim::Integer, energy_window,
                                ranges=nothing, dq::Real=0.02, qmax::Real=5.0, exclude=nothing, reciprocal_basis=nothing)
    qedges = _radial_edges(dq, qmax)
    state = Ref{Any}(nothing)
    eachslice_mdnorm(path_D, path_N; dim, ranges) do slice, _
        e = energy_axis(slice)
        if state[] === nothing
            n = length(qedges) - 1
            state[] = (; acc=(; D=zeros(n, size(slice.D, e)), N=zeros(n, size(slice.D, e))), ω=slice.centres[e],
                       ωedges=slice.edges[e])
        end
        mask = exclude === nothing ? nothing : _check_exclude(exclude(slice), slice)
        _accumulate_radial!(state[].acc, slice, mask, qedges, energy_window, reciprocal_basis)
    end
    s = state[]
    return _radial_result(s.acc, qedges, s.ω, s.ωedges, energy_window)
end

function _radial_edges(dq, qmax)
    isfinite(dq) && dq > 0 && isfinite(qmax) && qmax > dq || throw(ArgumentError("invalid |Q| binning"))
    return collect(0.0:Float64(dq):Float64(qmax))
end

"""
    fit_powder_lines(Q, R, τs; window=0.45, offsets=-0.03:0.001:0.03, widths=range(0.015, 0.2; length=40),
                     min_coverage=0.7)

Fit Gaussian powder lines at `τ (1 + ε)`, sharing the fractional offset `ε` and the width `σ_Q`, on a linear
baseline, over `minimum(τs) - window ≤ Q ≤ maximum(τs) + window`. Line areas are linear parameters; `ε` and `σ_Q`
are searched over `offsets` and `widths`. Returns `nothing` when fewer than `min_coverage` of the window's points
are finite. Otherwise returns `(ε, σQ, area, σarea, Q, model)` with least-squares area uncertainties.
"""
function fit_powder_lines(Q::AbstractVector{<:Real}, R::AbstractVector{<:Real}, τs; window::Real=0.45,
                          offsets=-0.03:0.001:0.03, widths=range(0.015, 0.2; length=40), min_coverage::Real=0.7)
    length(Q) == length(R) || throw(DimensionMismatch("Q and R must have the same length"))
    τ = collect(Float64, τs)
    isempty(τ) && throw(ArgumentError("supply at least one powder line"))
    lo, hi = minimum(τ) - window, maximum(τ) + window
    inside = [i for i in eachindex(Q) if lo <= Q[i] <= hi]
    keep = filter(i -> isfinite(R[i]), inside)
    (isempty(inside) || length(keep) < min_coverage * length(inside) || length(keep) <= length(τ) + 2) && return nothing
    x, y = Float64.(Q[keep]), Float64.(R[keep])
    xm = mean(x)
    best = (Inf, NaN, NaN, Float64[], zeros(0, 0))
    X = hcat(ones(length(x)), x .- xm, zeros(length(x), length(τ)))
    for ε in offsets, σ in widths
        for (k, t) in enumerate(τ)
            @. X[:, 2+k] = exp(-(x - t * (1 + ε))^2 / (2σ^2)) / (sqrt(2π) * σ)
        end
        θ = X \ y
        r = sum(abs2, y .- X * θ)
        r < best[1] && (best = (r, ε, σ, θ, copy(X)))
    end
    r, ε, σ, θ, Xb = best
    cov = (r / max(length(y) - size(Xb, 2), 1)) * inv(Xb' * Xb)
    return (; ε, σQ=σ, area=θ[3:end], σarea=sqrt.(diag(cov)[3:end]), Q=x, model=Xb * θ)
end

function _cubic_multiplicity(hkl)
    v = Int.(collect(hkl))
    images = Set{NTuple{3,Int}}()
    for p in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1)), s in Iterators.product((1, -1), (1, -1), (1, -1))
        push!(images, (s[1] * v[p[1]], s[2] * v[p[2]], s[3] * v[p[3]]))
    end
    return length(images)
end

"""
    powder_line_ratio(τ, multiplicity, F2, cell_volume, σ_inc_cell)
    powder_line_ratio(element, hkl; a)

Ratio of a powder line's `|Q|`-integrated elastic intensity to the flat incoherent elastic level of the same
material, for an ideal powder with `e^{-2W} ≈ 1`:

    (2π)³ m |F|² / (v₀ τ² σ_inc,cell)    (Å⁻¹)

with `|F|²` and `σ_inc,cell` per unit cell (barn), `v₀` in Å³ and `τ` in Å⁻¹. The second method treats a
single-element fcc metal with lattice parameter `a` (Å), using [`neutron_xs`](@ref) and the cubic multiplicity of
`hkl`; for aluminium at low temperature, `a ≈ 4.03` Å. Divide a fitted line area by the ratio to obtain that
material's flat elastic level.
"""
function powder_line_ratio(τ::Real, multiplicity::Real, F2::Real, cell_volume::Real, σ_inc_cell::Real)
    all(x -> isfinite(x) && x > 0, (τ, multiplicity, F2, cell_volume, σ_inc_cell)) ||
        throw(ArgumentError("powder-line parameters must be positive"))
    return (2π)^3 * multiplicity * F2 / (cell_volume * τ^2 * σ_inc_cell)
end

function powder_line_ratio(element::Union{Symbol,AbstractString}, hkl; a::Real)
    length(hkl) == 3 && all(x -> x isa Integer, hkl) || throw(ArgumentError("hkl must be three integers"))
    all(iseven, hkl) || all(isodd, hkl) || throw(ArgumentError("$hkl is extinct in an fcc lattice"))
    xs = neutron_xs(element)
    F2 = 16 * xs.b_coh^2 / 100     # fm² to barn
    τ = 2π * norm(hkl) / a
    return powder_line_ratio(τ, _cubic_multiplicity(hkl), F2, a^3, 4 * xs.σ_inc)
end

"""
    powder_line_flat(fit, ratios; detect=3) -> (flat, σflat, flats, σflats, detected)

Convert fitted line areas ([`fit_powder_lines`](@ref)) to flat elastic levels with the ratios of
[`powder_line_ratio`](@ref), and combine the lines detected above `detect` standard deviations by inverse-variance
weighting. With more than one line, `σflat` is at least their standard deviation, since texture need not affect
the lines alike. `flat` and `σflat` are `NaN` when no line is detected.
"""
function powder_line_flat(fit, ratios; detect::Real=3)
    fit === nothing && return (; flat=NaN, σflat=NaN, flats=Float64[], σflats=Float64[], detected=falses(0))
    length(ratios) == length(fit.area) || throw(DimensionMismatch("one ratio per fitted line required"))
    flats, σflats = fit.area ./ ratios, fit.σarea ./ ratios
    detected = BitVector(fit.area .> detect .* fit.σarea)
    any(detected) || return (; flat=NaN, σflat=NaN, flats, σflats, detected)
    w = 1 ./ σflats[detected] .^ 2
    flat = sum(w .* flats[detected]) / sum(w)
    σflat = 1 / sqrt(sum(w))
    count(detected) > 1 && (σflat = max(σflat, std(flats[detected])))
    return (; flat, σflat, flats, σflats, detected)
end

function _sum_rule_window(ωs, window, Δω)
    widths = Δω === nothing ? energy_bin_widths(ωs) : (Δω isa Real ? fill(Float64(Δω), length(ωs)) : collect(Float64, Δω))
    length(widths) == length(ωs) || throw(DimensionMismatch("Δω must match ωs"))
    js = [j for j in eachindex(ωs) if window[1] <= ωs[j] <= window[2]]
    isempty(js) && throw(ArgumentError("no energy bins in the sum-rule window"))
    return js, widths
end

"""
    first_moment_scale(S, σ, ωs, model; T, window, Δω=nothing, model_full=nothing, correction=1, correction_σ=0,
                       min_snr=5, min_fraction=0.97, low_quantile=0.25)

Scale from the first-moment sum rule, momentum by momentum. For each row `i` of the `nq × nω` arrays,

    M₁(i) = Σ_ω ω (1 - e^{-ω/k_B T}) S(i, ω) Δω

over `window`, which is the zero-temperature first moment. With the model `S_model` per magnetic ion on the same
grid, `κ_i = correction · M₁_model(i) / M₁(i)`, where `correction` is the ratio of the true to the model first
moment (for linear spin waves, the correction of the near-neighbour correlator). A row is used when it is finite
across the window, `M₁ > min_snr σ_M₁` (skipped if `σ` is `nothing`), and, if `model_full` (each row's full
model first moment) is given, the window holds at least `min_fraction` of it.

Returns `(κ, κ_lo, κ_q, M1, σM1, M1_model, used)`: `κ` is the median over the rows used, and
`κ_lo = (correction - correction_σ)/correction · quantile(κ_q, low_quantile)`.
"""
function first_moment_scale(S::AbstractMatrix{<:Real}, σ, ωs, model::AbstractMatrix{<:Real}; T::Real, window,
                            Δω=nothing, model_full=nothing, correction::Real=1.0, correction_σ::Real=0.0,
                            min_snr::Real=5.0, min_fraction::Real=0.97, low_quantile::Real=0.25)
    size(S) == size(model) || throw(DimensionMismatch("S and model must have the same shape"))
    σ === nothing || size(σ) == size(S) || throw(DimensionMismatch("σ must match S"))
    size(S, 2) == length(ωs) || throw(DimensionMismatch("S must have one column per energy"))
    model_full === nothing || length(model_full) == size(S, 1) || throw(DimensionMismatch("one model_full per row"))
    T > 0 && correction > 0 && 0 <= correction_σ < correction || throw(ArgumentError("invalid T or correction"))
    js, Δ = _sum_rule_window(ωs, window, Δω)
    kT = KB_MEV_PER_K * T
    weight = [ωs[j] * (1 - exp(-ωs[j] / kT)) * Δ[j] for j in js]
    moment(X, i) = all(j -> isfinite(X[i, j]), js) ? sum(weight[k] * X[i, j] for (k, j) in enumerate(js)) : NaN
    M1 = [moment(S, i) for i in axes(S, 1)]
    σM1 = σ === nothing ? fill(NaN, size(S, 1)) :
          [all(j -> isfinite(σ[i, j]), js) ? sqrt(sum((weight[k] * σ[i, j])^2 for (k, j) in enumerate(js))) : NaN
           for i in axes(S, 1)]
    M1_model = [moment(model, i) for i in axes(S, 1)]
    used = BitVector(map(axes(S, 1)) do i
        isfinite(M1[i]) && M1[i] > 0 && isfinite(M1_model[i]) || return false
        σ === nothing || (isfinite(σM1[i]) && M1[i] > min_snr * σM1[i]) || return false
        model_full === nothing || M1_model[i] >= min_fraction * model_full[i] || return false
        return true
    end)
    κ_q = [used[i] ? correction * M1_model[i] / M1[i] : NaN for i in axes(S, 1)]
    values = κ_q[used]
    isempty(values) && return (; κ=NaN, κ_lo=NaN, κ_q, M1, σM1, M1_model, used)
    κ_lo = (correction - correction_σ) / correction * quantile(values, low_quantile)
    return (; κ=median(values), κ_lo, κ_q, M1, σM1, M1_model, used)
end

"""
    total_moment_scale(S, model, ωs; window, Δω=nothing, correction=1, correction_σ=0)

Scale from a total-moment sum rule over a zone grid. `S` (measured, relative) and `model` (per magnetic ion) have
the energy as their last axis; `S` is `NaN` in unmeasured cells. Over the cells and energies of `window` where `S`
is finite,

    κ = correction · Σ model Δω / Σ S Δω,

where `correction` is the ratio of the true to the model weight in the channel. `S` on a zone grid, with every bin
at all its symmetry images, follows from [`magnetic_histogram`](@ref), [`fold_to_zone`](@ref) and
[`intensity`](@ref).

Returns `(κ, κ_lo, covered_fraction)`: `covered_fraction` is the share of the model weight in the window that
falls in measured cells, and `κ_lo = (correction - correction_σ)/correction · κ`.
"""
function total_moment_scale(S::AbstractArray{<:Real}, model::AbstractArray{<:Real}, ωs; window, Δω=nothing,
                            correction::Real=1.0, correction_σ::Real=0.0)
    size(S) == size(model) || throw(DimensionMismatch("S and model must have the same shape"))
    size(S, ndims(S)) == length(ωs) || throw(DimensionMismatch("the last axis of S must be energy"))
    correction > 0 && 0 <= correction_σ < correction || throw(ArgumentError("invalid correction"))
    js, Δ = _sum_rule_window(ωs, window, Δω)
    measured_model = all_model = measured = 0.0
    for i in CartesianIndices(S)
        j = i[ndims(S)]
        j in js || continue
        m = model[i] * Δ[j]
        all_model += m
        isfinite(S[i]) || continue
        measured_model += m
        measured += S[i] * Δ[j]
    end
    κ = measured > 0 ? correction * measured_model / measured : NaN
    return (; κ, κ_lo=(correction - correction_σ) / correction * κ, covered_fraction=measured_model / all_model)
end

"""
    compare_scales(estimates; reference=nothing)

Tabulate scale estimates given as named tuples `(name, κ)` or `(name, κ, κ_lo)`. Each row gives `used`, the low end
when it is finite and `κ` otherwise, and `ratio = κ / κ_reference`, where `reference` names one estimate (`NaN`
without one). See [`select_scale`](@ref).
"""
function compare_scales(estimates; reference=nothing)
    rows = map(estimates) do e
        κ = Float64(e.κ)
        lo = hasproperty(e, :κ_lo) ? Float64(e.κ_lo) : NaN
        (; name=String(e.name), κ, κ_lo=lo, used=isfinite(lo) ? lo : κ)
    end
    ref = reference === nothing ? NaN : begin
        k = findfirst(r -> r.name == String(reference), rows)
        k === nothing && throw(ArgumentError("no scale named $reference"))
        rows[k].κ
    end
    return [merge(r, (; ratio=r.κ / ref)) for r in rows]
end

"""
    select_scale(table; among=nothing)

Return the row of [`compare_scales`](@ref) with the smallest finite `used` value, optionally among the names in
`among`. The smallest scale gives the smallest `S_mag` and so the most conservative lower bound on the QFI.
"""
function select_scale(table; among=nothing)
    names = among === nothing ? nothing : Set(String.(collect(among)))
    rows = [r for r in table if isfinite(r.used) && (names === nothing || r.name in names)]
    isempty(rows) && throw(ArgumentError("no finite scale to select"))
    return rows[argmin([r.used for r in rows])]
end
