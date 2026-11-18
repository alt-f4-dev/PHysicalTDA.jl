# Branch identity through unresolved junctions of two ridges.
#
# Near a junction two branches are the eigenvalues of a two-level effective Hamiltonian
#     H(q) = [e₁(q) g; g e₂(q)],   E± = m ± r,   m = (e₁ + e₂)/2,  d = (e₁ − e₂)/2,  r = √(d² + g²),
# with smooth (polynomial) diabatic energies e₁, e₂ and a constant coupling g. The upper eigenvector has the weight
# c² = (1 + d/r)/2 on |1⟩, so the eigenbranches carry W± = A₁c² + A₂(1 − c²) and A₁(1 − c²) + A₂c² for diabatic
# amplitudes A₁, A₂ (no interference term). Each branch is a Gaussian of FWHM √(Γω² + (∂E/∂q Γq,eff)²), the
# momentum resolution projected along the local slope. The model is fitted to the intensity map itself.

const _JUNCTION_LABELLINGS = (:auto, :diabatic, :adiabatic)

"""
    JunctionFit

A two-level fit of one junction ([`fit_junction`](@ref)).

- `q`: momentum range of the fitted `columns`.
- `class`: connectivity of the best fit, `:cross` (the diabatic energies pass through each other), `:avoid`
  (a significant gap; the branches keep their energy order) or `:bounce` (no exchange across the window).
- `verdict`: `class`, or `:ambiguous` when another class is within `Δχ2_min` of the best fit.
- `Δχ2`: χ² margin to the best fit of another class, divided by the best fit's reduced χ² when that exceeds 1.
- `gap`, `σgap`: the fitted gap `2|g|` and its error; `χ2ν`: reduced χ² of the best fit.
- `Γq`: the fitted momentum resolution projected on the cut; `Γω`: the energy resolution used.
"""
struct JunctionFit
    q::NTuple{2,Float64}
    columns::UnitRange{Int}
    window::NTuple{2,Float64}
    class::Symbol
    verdict::Symbol
    Δχ2::Float64
    gap::Float64
    σgap::Float64
    χ2ν::Float64
    Γq::Float64
    Γω::Float64
    params::Vector{Float64}
    σparams::Vector{Float64}
    qc::Float64
    halfwidth::Float64
    degree::Int
    ωmean::Float64
    classes::Vector{Symbol}
    χ2::Vector{Float64}
end

function Base.show(io::IO, f::JunctionFit)
    r(x) = round(x; sigdigits=3)
    print(io, "JunctionFit(q ∈ [", r(f.q[1]), ", ", r(f.q[2]), "]: ", f.verdict, ", Δχ² = ", r(f.Δχ2), ", gap = ",
          r(f.gap), " ± ", r(f.σgap), ", χ²ν = ", r(f.χ2ν), ")")
end

_junction_layout(de) = (e1=1:de+1, e2=de+2:2de+2, g=2de+3, a1=2de+4:2de+6, a2=2de+7:2de+9, b=2de+10:2de+12, rq=2de+13)
_jpoly(c, x) = (s = 0.0; for k in length(c):-1:1; s = s * x + c[k]; end; s)
_jdpoly(c, x) = (s = 0.0; for k in length(c):-1:2; s = s * x + (k - 1) * c[k]; end; s)
_jline(δ, Γω, Γs) = (s2 = (Γω^2 + Γs^2) / (8log(2)); exp(-δ^2 / (2s2)) / sqrt(2π * s2))

function _junction_model!(out, p, X, W, Γω, Γq, sx, ω̄, de)
    L = _junction_layout(de)
    e1, e2, A1, A2 = view(p, L.e1), view(p, L.e2), view(p, L.a1), view(p, L.a2)
    g, b0, bx, bω = p[L.g], p[L.b[1]], p[L.b[2]], p[L.b[3]]
    Γs = Γq * exp(p[L.rq]) / sx
    @inbounds for k in eachindex(X)
        x, ω = X[k], W[k]
        a, b = _jpoly(e1, x), _jpoly(e2, x)
        m, d = (a + b) / 2, (a - b) / 2
        r = sqrt(d^2 + g^2)
        c2 = r > 0 ? (1 + d / r) / 2 : 0.5
        dm, dd = (_jdpoly(e1, x) + _jdpoly(e2, x)) / 2, (_jdpoly(e1, x) - _jdpoly(e2, x)) / 2
        t = r > 0 ? d * dd / r : 0.0
        w1 = exp(A1[1] + A1[2] * x + A1[3] * x^2)
        w2 = exp(A2[1] + A2[2] * x + A2[3] * x^2)
        out[k] = b0 + bx * x + bω * (ω - ω̄) +
                 (w1 * c2 + w2 * (1 - c2)) * _jline(ω - m - r, Γω, (dm + t) * Γs) +
                 (w1 * (1 - c2) + w2 * c2) * _jline(ω - m + r, Γω, (dm - t) * Γs)
    end
    return out
end

# Levenberg–Marquardt with a forward-difference Jacobian; `fixed` parameters stay at their start values. The
# covariance is scaled by the reduced χ² when the model misfits.
function _junction_lm(p0, X, W, Y, S, Γω, Γq, sx, ω̄, de; fixed=Int[], iters=200, rtol=1e-7)
    model!(out, p) = _junction_model!(out, p, X, W, Γω, Γq, sx, ω̄, de)
    p = copy(p0); free = setdiff(eachindex(p), fixed)
    f = similar(Y); model!(f, p)
    χ2 = sum(abs2, (f .- Y) ./ S); λ = 1e-2
    J = zeros(length(Y), length(free)); fp = similar(Y)
    function jacobian!()
        for (j, i) in enumerate(free)
            h = 1e-6 * max(abs(p[i]), 1e-3)
            p[i] += h; model!(fp, p); p[i] -= h
            @views J[:, j] .= (fp .- f) ./ (h .* S)
        end
    end
    for _ in 1:iters
        jacobian!()
        A = J' * J; grad = J' * ((f .- Y) ./ S)
        accepted = false
        for _ in 1:12
            δ = -((A + λ * Diagonal(diag(A) .+ 1e-12)) \ grad)
            trial = copy(p); trial[free] .+= δ
            model!(fp, trial)
            χ2new = sum(abs2, (fp .- Y) ./ S)
            if isfinite(χ2new) && χ2new < χ2
                gain = χ2 - χ2new
                p, χ2 = trial, χ2new; f .= fp
                λ = max(λ / 5, 1e-10); accepted = true
                gain < rtol * χ2 && (accepted = false)   # converged
                break
            end
            λ *= 8
        end
        accepted || break
    end
    jacobian!()
    ν = max(length(Y) - length(free), 1)
    C = pinv(J' * J) .* max(χ2 / ν, 1.0)
    σ = zeros(length(p)); σ[free] .= sqrt.(max.(diag(C), 0.0))
    return (; p, χ2, σ, ν)
end

"""
    fit_junction(I, σ, qs, ωs, columns, window, seeds; resolution, degree=3, fit_resolution=true, nσ=3,
                 Δχ2_min=4, valid=nothing) -> JunctionFit

Fit the two-level junction model to `I[columns, ω ∈ window]` with errors `σ`. `seeds` are the lower and upper
ridge on each side, `(left_low, left_high, right_low, right_high)`, each with fields `q`, `ω` and `I` (peak
height), from resolved columns next to the junction; they only start the fit. `resolution = (Γq, Γω)` (FWHM).

The diabatic energies are polynomials of `degree` in the momentum, the amplitudes exponentials of quadratics, and
the background linear in momentum and energy. With `fit_resolution`, the scale of `Γq` is fitted (the projection
of the momentum resolution on a cut is rarely known exactly). Fits start from both pairings of the seeds, with
the coupling fixed at zero, at half the seed gap, and at `Γω`. A fit is `:avoid` when `2|g| ≥ nσ σ(2g)` and
`2|g| ≥ 0.1 Γω`; otherwise `:cross` or `:bounce` from the sign of `e₁ − e₂` at the ends of the window. The
lowest χ² wins; its `verdict` is `:ambiguous` when a fit of another class is within `Δχ2_min` (after dividing
by the reduced χ² of the best fit when that exceeds 1).
"""
function fit_junction(I::AbstractMatrix{<:Real}, σ::AbstractMatrix{<:Real}, qs, ωs, columns::AbstractUnitRange,
                      window, seeds; resolution, degree::Integer=3, fit_resolution::Bool=true, nσ::Real=3.0,
                      Δχ2_min::Real=4.0, valid=nothing)
    size(I) == size(σ) == (length(qs), length(ωs)) || throw(DimensionMismatch("I, σ, qs and ωs must agree"))
    first(columns) >= 1 && last(columns) <= length(qs) && length(columns) >= 3 ||
        throw(ArgumentError("columns must lie inside the momentum axis and span at least three columns"))
    degree >= 1 || throw(ArgumentError("degree must be at least 1"))
    Γq, Γω = Float64.(resolution)
    Γq > 0 && Γω > 0 || throw(ArgumentError("resolution must be positive"))
    length(seeds) == 4 || throw(ArgumentError("seeds are (left_low, left_high, right_low, right_high)"))
    measured = valid === nothing ? trues(size(I)) : BitMatrix(valid)
    qc = (qs[first(columns)] + qs[last(columns)]) / 2
    sx = abs(qs[last(columns)] - qs[first(columns)]) / 2
    X, W, Y, S = Float64[], Float64[], Float64[], Float64[]
    for i in columns, j in eachindex(ωs)
        window[1] <= ωs[j] <= window[2] && measured[i, j] && isfinite(I[i, j]) && σ[i, j] > 0 || continue
        push!(X, (qs[i] - qc) / sx); push!(W, ωs[j]); push!(Y, I[i, j]); push!(S, σ[i, j])
    end
    de = Int(degree)
    L = _junction_layout(de)
    length(Y) > 4(2de + 13) || throw(ArgumentError("too few measured bins in the junction window"))
    ω̄ = sum(W) / length(W)
    xs(a, b) = [(a.q .- qc) ./ sx; (b.q .- qc) ./ sx]
    function energy_seed(a, b)
        x = xs(a, b)
        M = [x[i]^k for i in eachindex(x), k in 0:min(de, length(x) - 1)]
        c = M \ [a.ω; b.ω]
        return [c; zeros(de + 1 - length(c))]
    end
    function amplitude_seed(a, b)
        x = xs(a, b); v = log.(max.([a.I; b.I] .* (1.0645 * Γω), 1e-300))
        c = [ones(length(x)) x] \ v
        return [c[1], c[2], 0.0]
    end
    background = sort(Y)[max(1, round(Int, 0.1length(Y)))]
    lo_l, hi_l, lo_r, hi_r = seeds
    starts = Vector{Vector{Float64}}()
    for (u, v, w, z) in ((lo_l, hi_r, hi_l, lo_r), (hi_l, hi_r, lo_l, lo_r))   # crossing, then energy-ordered pairs
        e1, e2 = energy_seed(u, v), energy_seed(w, z)
        half = abs(e1[1] - e2[1]) / 2
        for g0 in (0.0, max(half, 0.05Γω), Γω)
            push!(starts, [e1; e2; g0; amplitude_seed(u, v); amplitude_seed(w, z); background; 0.0; 0.0; 0.0])
        end
    end
    fits = map(starts) do p0
        fixed = Int[]
        p0[L.g] == 0 && push!(fixed, L.g)
        fit_resolution || push!(fixed, L.rq)
        _junction_lm(p0, X, W, Y, S, Γω, Γq, sx, ω̄, de; fixed)
    end
    function classify(fit)
        g2, σg2 = 2abs(fit.p[L.g]), 2fit.σ[L.g]
        σg2 > 0 && g2 >= nσ * σg2 && g2 >= 0.1Γω && return :avoid
        d(x) = _jpoly(view(fit.p, L.e1), x) - _jpoly(view(fit.p, L.e2), x)
        return sign(d(-1.0)) == sign(d(1.0)) ? :bounce : :cross
    end
    classes = classify.(fits)
    χ2 = [f.χ2 for f in fits]
    k = argmin(χ2)
    best = fits[k]
    others = [χ2[i] for i in eachindex(fits) if classes[i] !== classes[k]]
    Δ = isempty(others) ? Inf : (minimum(others) - best.χ2) / max(best.χ2 / best.ν, 1.0)
    verdict = Δ < Δχ2_min ? :ambiguous : classes[k]
    return JunctionFit((Float64(qs[first(columns)]), Float64(qs[last(columns)])), columns, (Float64(window[1]), Float64(window[2])), classes[k], verdict, Δ,
                       2abs(best.p[L.g]), 2best.σ[L.g], best.χ2 / best.ν, Γq * exp(best.p[L.rq]), Γω, best.p, best.σ,
                       qc, sx, de, ω̄, classes, χ2)
end

"""
    diabatic_energies(f::JunctionFit, q) -> (e₁, e₂, g)

The fitted diabatic energies and coupling at momentum `q`; the eigenbranches are `(e₁+e₂)/2 ∓ √(((e₁−e₂)/2)² + g²)`.
"""
function diabatic_energies(f::JunctionFit, q::Real)
    L = _junction_layout(f.degree); x = (q - f.qc) / f.halfwidth
    return _jpoly(view(f.params, L.e1), x), _jpoly(view(f.params, L.e2), x), f.params[L.g]
end

"""
    junction_model(f::JunctionFit, qs, ωs) -> Matrix

The fitted intensity on the grid `qs × ωs` (for residual maps).
"""
function junction_model(f::JunctionFit, qs, ωs)
    X = [(q - f.qc) / f.halfwidth for q in qs, _ in ωs]; W = [ω for _ in qs, ω in ωs]
    out = similar(X)
    Γq = f.Γq / exp(f.params[_junction_layout(f.degree).rq])
    _junction_model!(vec(out), f.params, vec(X), vec(W), f.Γω, Γq, f.halfwidth, f.ωmean, f.degree)
    return out
end

"""
    BranchResolution

Branches assembled from all child tracks of a [`ChildRidgeResult`](@ref) by [`resolve_branches`](@ref).

- `q`: momenta; `ω`: branch energies (branches × momenta, `NaN` where a branch is absent).
- `from_model`: entries inside a junction, taken from the fitted diabatic energies (after an exchange) or the
  fitted eigenvalues (otherwise).
- `junctions`: every fitted junction; `pairs`: the two branches it joined (labelled by their state left of the
  junction); `unresolved`: the momentum range of its unresolved columns.
"""
struct BranchResolution
    q::Vector{Float64}
    ω::Matrix{Float64}
    from_model::BitMatrix
    junctions::Vector{JunctionFit}
    pairs::Vector{NTuple{2,Int}}
    unresolved::Vector{NTuple{2,Float64}}
end

function Base.show(io::IO, r::BranchResolution)
    j = ["$(f.verdict) at q ∈ [$(round(u[1]; sigdigits=3)), $(round(u[2]; sigdigits=3))]" for (f, u) in zip(r.junctions, r.unresolved)]
    print(io, "BranchResolution($(size(r.ω, 1)) branches, $(length(r.junctions)) junctions", isempty(j) ? "" : ": ",
          join(j, "; "), ")")
end

# Ridge points per momentum column from every child track: energies sorted, points closer than Γω merged.
function _column_points(tracks, nq, Γω)
    P = [Float64[] for _ in 1:nq]
    for i in 1:nq
        e = sort!([Float64(T[b, i]) for T in tracks for b in axes(T, 1) if isfinite(T[b, i])])
        isempty(e) && continue
        group = [e[1]]
        for x in @view e[2:end]
            if x - group[end] < Γω
                push!(group, x)
            else
                push!(P[i], sum(group) / length(group)); group = [x]
            end
        end
        push!(P[i], sum(group) / length(group))
    end
    return P
end

# Predicted energy of a branch at column i from its last `k` assigned columns (linear in q).
function _branch_prediction(ω, qs, b, upto, i; k=3)
    cols = Int[]
    for c in upto:-1:1
        isfinite(ω[b, c]) && push!(cols, c)
        length(cols) == k && break
    end
    isempty(cols) && return NaN
    length(cols) == 1 && return ω[b, cols[1]]
    q, e = qs[cols], ω[b, cols]
    slope = sum((q .- sum(q) / length(q)) .* (e .- sum(e) / length(e))) / max(sum(abs2, q .- sum(q) / length(q)), eps())
    return sum(e) / length(e) + slope * (qs[i] - sum(q) / length(q))
end

"""
    resolve_branches(r::ChildRidgeResult, σ; resolution, valid=nothing, labelling=:auto, min_side=4,
                     max_side=nothing, max_junction=nothing, max_approach=4, max_gap=nothing,
                     window_pad=0.75, kw...) -> BranchResolution

Assemble branches from all child tracks of `r` (independent of how the tracks are split into parents and of the
tracker's labels) and decide their identity at every unresolved junction from a two-level fit to the data.

1. In each momentum column the ridge points of all tracks are merged when closer than `Γω`. Points are linked to
   the linear continuation of the active branches, and of branches that ended within `max_gap` columns (default
   `max_junction`), by nearest energy within `max(4Γω, 2 Δq params.max_slope)`.
2. For every pair of branches, a run of at most `max_junction` columns (default `6Γq`) where they are unresolved
   (closer than `Γω`, or one missing), flanked on both sides by at least `min_side` columns where both are
   resolved and no farther apart than `max_approach Γω`, is a junction. It is fitted with
   [`fit_junction`](@ref) over at most `max_side` flanking columns per side (default `3Γq`) and the seed energies
   `± window_pad`; further keywords go to it. The labels right of the junction are exchanged when the fit and
   the linking of step 1 disagree, and the junction columns are filled from the fit.

`labelling` decides the connectivity at a junction: `:auto` exchanges after a confident `:cross`; `:diabatic`
also after an `:avoid`, so that branches follow the uncoupled modes through a gap (appropriate when the gap's
mixing region, `gap/|Δv|`, is below the momentum resolution); `:adiabatic` never. Ambiguous junctions keep
energy order.
"""
function resolve_branches(r::ChildRidgeResult, σ::AbstractMatrix{<:Real}; resolution, valid=nothing,
                          labelling::Symbol=:auto, min_side::Integer=4, max_side=nothing, max_junction=nothing,
                          max_approach::Real=4.0, max_gap=nothing, window_pad::Real=0.75, kw...)
    labelling in _JUNCTION_LABELLINGS || throw(ArgumentError("labelling must be :auto, :diabatic or :adiabatic"))
    min_side >= 2 || throw(ArgumentError("min_side must be at least 2"))
    c = r.candidates
    I = Matrix{Float64}(c.I); qs = Float64.(c.q); ωs = Float64.(c.ω); nq = length(qs)
    size(σ) == size(I) || throw(DimensionMismatch("σ must match the intensity map"))
    measured = valid === nothing ? trues(size(I)) : BitMatrix(valid)
    size(measured) == size(I) || throw(DimensionMismatch("valid must match the intensity map"))
    I[.!measured] .= NaN
    Γq, Γω = Float64(resolution[1]), Float64(resolution[2])
    dq = abs(qs[2] - qs[1])
    side_max = max_side === nothing ? max(min_side, round(Int, 3Γq / dq)) : Int(max_side)
    run_max = max_junction === nothing ? max(2, round(Int, 6Γq / dq)) : Int(max_junction)
    gap = max_gap === nothing ? run_max : Int(max_gap)
    slope = Float64(c.params.max_slope)

    # 1. link points into branches
    P = _column_points(r.tracks, nq, Γω)
    ω = zeros(0, nq)
    for i in 1:nq
        isempty(P[i]) && continue
        pool = [b for b in axes(ω, 1) if (l = findlast(isfinite, @view ω[b, 1:i-1]); l !== nothing && i - l <= gap + 1)]
        pred = [_branch_prediction(ω, qs, b, i - 1, i) for b in pool]
        cost = sort!([(abs(P[i][m] - pred[a]), a, m) for a in eachindex(pool), m in eachindex(P[i]) if isfinite(pred[a])])
        used_a, used_m = falses(length(pool)), falses(length(P[i]))
        for (d, a, m) in cost
            l = findlast(isfinite, @view ω[pool[a], 1:i-1])
            (used_a[a] || used_m[m] || d > max(4Γω, 2abs(qs[i] - qs[l]) * slope)) && continue
            used_a[a] = used_m[m] = true
            ω[pool[a], i] = P[i][m]
        end
        for m in eachindex(P[i])
            used_m[m] && continue
            ω = vcat(ω, fill(NaN, 1, nq)); ω[end, i] = P[i][m]
        end
    end
    from_model = falses(size(ω))

    # 2. junctions of every pair, in momentum order
    nearest(e) = argmin(abs.(ωs .- e))
    peak(i, e) = (v = I[i, nearest(e)]; isfinite(v) ? v : 0.0)
    resolved(a, b, i) = isfinite(ω[a, i]) && isfinite(ω[b, i]) && abs(ω[a, i] - ω[b, i]) >= Γω
    function junctions_of(a, b)
        runs = Tuple{UnitRange{Int},Vector{Int},Vector{Int}}[]
        i = 1
        while i <= nq
            if !resolved(a, b, i) || i == 1
                i += 1; continue
            end
            # i is resolved; look for an unresolved run starting at i + 1
            j = i + 1
            while j <= nq && !resolved(a, b, j) && j - i - 1 < run_max; j += 1; end
            if j > i + 1 && j <= nq && resolved(a, b, j)
                left = Int[]; k = i
                while k >= 1 && resolved(a, b, k) && length(left) < side_max; push!(left, k); k -= 1; end
                right = Int[]; k = j
                while k <= nq && resolved(a, b, k) && length(right) < side_max; push!(right, k); k += 1; end
                sort!(left)
                close = abs(ω[a, i] - ω[b, i]) <= max_approach * Γω && abs(ω[a, j] - ω[b, j]) <= max_approach * Γω
                # a junction needs both branches present on both flanks, and at least one unresolved column where
                # they merge (both missing everywhere is a data gap, not a junction)
                merged = any(k -> isfinite(ω[a, k]) || isfinite(ω[b, k]), i+1:j-1)
                length(left) >= min_side && length(right) >= min_side && close && merged &&
                    push!(runs, (i+1:j-1, left, right))
            end
            i = j
        end
        return runs
    end
    events = Tuple{Int,Int,Int}[]
    for a in axes(ω, 1), b in a+1:size(ω, 1)
        for (J, _, _) in junctions_of(a, b)
            push!(events, (first(J), a, b))
        end
    end
    sort!(events)
    fits = JunctionFit[]; pairs = NTuple{2,Int}[]; spans = NTuple{2,Float64}[]
    for (_, a, b) in events
        found = junctions_of(a, b)        # recomputed: earlier exchanges may have moved it
        isempty(found) && continue
        J, left, right = found[argmin([abs(first(x[1]) - first(found[1][1])) for x in found])]
        any(f -> f[1] == (a, b) && f[2] == (qs[first(J)], qs[last(J)]), zip(pairs, spans)) && continue
        lo(cols, x, y) = [min(ω[x, k], ω[y, k]) for k in cols]
        hi(cols, x, y) = [max(ω[x, k], ω[y, k]) for k in cols]
        seg(cols, e) = (q=qs[cols], ω=e, I=[peak(k, v) for (k, v) in zip(cols, e)])
        seeds = (seg(left, lo(left, a, b)), seg(left, hi(left, a, b)), seg(right, lo(right, a, b)), seg(right, hi(right, a, b)))
        E = [lo(left, a, b); hi(left, a, b); lo(right, a, b); hi(right, a, b)]
        f = fit_junction(I, σ, qs, ωs, first(left):last(right), (minimum(E) - window_pad, maximum(E) + window_pad), seeds;
                         resolution=(Γq, Γω), valid=measured, kw...)
        lower_left = ω[a, last(left)] < ω[b, last(left)] ? a : b
        upper_left = lower_left == a ? b : a
        push!(fits, f); push!(pairs, (lower_left, upper_left)); push!(spans, (qs[first(J)], qs[last(J)]))
        exchange = labelling === :adiabatic ? false :
                   labelling === :diabatic ? f.verdict in (:cross, :avoid) : f.verdict === :cross
        # linking of step 1: does the left-lower branch continue as the right-lower one?
        linked_exchange = ω[lower_left, first(right)] > ω[upper_left, first(right)]
        if exchange != linked_exchange
            tail = first(J):nq
            ω[[a, b], tail] .= ω[[b, a], tail]
            from_model[[a, b], tail] .= from_model[[b, a], tail]
        end
        e1L, e2L, _ = diabatic_energies(f, qs[first(left)])
        for k in J
            e1, e2, g = diabatic_energies(f, qs[k])
            if exchange
                x, y = e1L < e2L ? (e1, e2) : (e2, e1)
            else
                m, d = (e1 + e2) / 2, (e1 - e2) / 2
                x, y = m - sqrt(d^2 + g^2), m + sqrt(d^2 + g^2)
            end
            ω[lower_left, k], ω[upper_left, k] = x, y
            from_model[lower_left, k] = from_model[upper_left, k] = true
        end
    end
    return BranchResolution(qs, ω, from_model, fits, pairs, spans)
end

"""
    branch_supports(b::BranchResolution, ωs; halfwidth, valid=nothing) -> Vector{BitMatrix}

One support per branch, `|ω − ω_branch(q)| ≤ halfwidth`, restricted to measured bins. Supports of crossing
branches overlap at the crossing.
"""
function branch_supports(b::BranchResolution, ωs; halfwidth::Real, valid=nothing)
    halfwidth > 0 || throw(ArgumentError("halfwidth must be positive"))
    measured = valid === nothing ? trues(length(b.q), length(ωs)) : BitMatrix(valid)
    return [BitMatrix([isfinite(b.ω[k, i]) && abs(ωs[j] - b.ω[k, i]) <= halfwidth && measured[i, j]
                       for i in eachindex(b.q), j in eachindex(ωs)]) for k in axes(b.ω, 1)]
end

"""
    branch_trajectories(I, b::BranchResolution, ωs; halfwidth, valid=nothing, centre=:weighted)
        -> (; q, ω, imputed, supports)

Trajectories of resolved branches: [`mask_dispersions`](@ref) inside [`branch_supports`](@ref), with the
junction columns taken from the fit and flagged in `imputed`. The result plots with
[`plot_ridge_workflow`](@ref) (`trajectories=`, and `supports=` for the supports).
"""
function branch_trajectories(I::AbstractMatrix{<:Real}, b::BranchResolution, ωs; halfwidth::Real, valid=nothing,
                             centre::Symbol=:weighted)
    supports = branch_supports(b, ωs; halfwidth, valid)
    t = mask_dispersions(I, supports, b.q, ωs; centre, valid)
    ω = copy(t.ω); ω[b.from_model] .= b.ω[b.from_model]
    return (q=t.q, ω=ω, imputed=t.imputed .| b.from_model, supports=supports)
end
