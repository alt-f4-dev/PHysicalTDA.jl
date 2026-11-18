function _ridge_kernel(axis, fwhm)
    isfinite(fwhm) && fwhm >= 0 || throw(ArgumentError("resolution widths must be finite and nonnegative"))
    σ = fwhm / sqrt(8log(2))
    return map(eachindex(axis)) do i
        iszero(σ) && return (i:i, [1.0])
        bins = searchsortedfirst(axis, axis[i] - 3σ):searchsortedlast(axis, axis[i] + 3σ)
        (bins, [exp(-0.5 * ((axis[j] - axis[i]) / σ)^2) for j in bins])
    end
end

# Carry the numerator, valid weight and variance through both separable passes.
function _smooth_ridge_measurements(I, σ, valid, qkernel, ωkernel)
    nq, nω = size(I)
    signal = zeros(nq, nω)
    weight = zeros(nq, nω)
    variance = zeros(nq, nω)
    for j in 1:nω, i in 1:nq
        bins, weights = qkernel[i]
        d = w = v = 0.0
        for (k, a) in zip(bins, weights)
            valid[k, j] || continue
            d += a * I[k, j]
            w += a
            v += abs2(a * σ[k, j])
        end
        signal[i, j], weight[i, j], variance[i, j] = d, w, v
    end
    S, error = fill(NaN, nq, nω), fill(NaN, nq, nω)
    ds, ws, vs = zeros(nq), zeros(nq), zeros(nq)
    for j in 1:nω
        fill!(ds, 0.0); fill!(ws, 0.0); fill!(vs, 0.0)
        bins, weights = ωkernel[j]
        for (k, a) in zip(bins, weights), i in 1:nq
            ds[i] += a * signal[i, k]
            ws[i] += a * weight[i, k]
            vs[i] += a^2 * variance[i, k]
        end
        for i in 1:nq
            valid[i, j] && ws[i] > 0 || continue
            S[i, j] = ds[i] / ws[i]
            error[i, j] = sqrt(vs[i]) / ws[i]
        end
    end
    return S, error
end

function _ridge_eligible(S, error, z, fraction)
    maxima = zeros(size(S, 1))
    for j in axes(S, 2), i in axes(S, 1)
        x = S[i, j]
        isfinite(x) && (maxima[i] = max(maxima[i], x))
    end
    return BitMatrix([S[i, j] > 0 && S[i, j] >= max(z * error[i, j], fraction * maxima[i])
                      for i in axes(S, 1), j in axes(S, 2)])
end

@inline function _visit_ridge_neighbors(f::F, i, j, nq, nω, adjacent) where {F}
    j > 1 && f(i, j-1)
    j < nω && f(i, j+1)
    for k in adjacent[j]
        i > 1 && f(i-1, k)
        i < nq && f(i+1, k)
    end
end

function _assign_ridge_supports(eligible, seeds, adjacent, minpix)
    nq, nω = size(eligible)
    owner = zeros(Int, nq, nω)
    queue = CartesianIndex{2}[]
    for (parent, mask) in pairs(seeds), i in CartesianIndices(eligible)
        eligible[i] && mask[i] && owner[i] == 0 || continue
        owner[i] = parent
        push!(queue, i)
    end
    first = 1
    while first <= length(queue)
        point = queue[first]
        i, j = Tuple(point)
        parent = owner[point]
        first += 1
        _visit_ridge_neighbors(i, j, nq, nω, adjacent) do k, l
            if eligible[k, l] && owner[k, l] == 0
                owner[k, l] = parent
                push!(queue, CartesianIndex(k, l))
            end
        end
    end
    minpix <= 1 && return owner
    seen = falses(nq, nω)
    for point in CartesianIndices(owner)
        parent = owner[point]
        parent > 0 && !seen[point] || continue
        empty!(queue)
        push!(queue, point)
        seen[point] = true
        first = 1
        while first <= length(queue)
            i, j = Tuple(queue[first])
            first += 1
            _visit_ridge_neighbors(i, j, nq, nω, adjacent) do k, l
                if owner[k, l] == parent && !seen[k, l]
                    seen[k, l] = true
                    push!(queue, CartesianIndex(k, l))
                end
            end
        end
        if length(queue) < minpix
            for i in queue
                owner[i] = 0
            end
        end
    end
    return owner
end

function _close_ridge_supports!(owner, valid, width)
    width == 0 && return owner
    for i in axes(owner, 1)
        left = 0
        for j in axes(owner, 2)
            if !valid[i, j]
                left = 0
            elseif owner[i, j] > 0
                if left > 0 && j - left <= width + 1 && owner[i, left] == owner[i, j]
                    @views owner[i, left+1:j-1] .= owner[i, j]
                end
                left = j
            end
        end
    end
    return owner
end

"""
    refine_supports(I, σ, valid, qs, ωs, coarse_masks; resolution, z=3.0, fraction=0.15,
                    bridge=nothing, close=2, minpix=10, rebin=nothing)

Refine existing parent supports using Gaussian FWHM widths `resolution = (Γq, Γω)` in axis units. Smooth only finite valid measurements and propagate independent input errors. Retain positive bins above `z` standard deviations and a fraction of the row maximum. Grow from the supplied parent seeds along adjacent momentum rows with at most `bridge` energy displacement per row (default the larger of the mean energy spacing and `Γω / 2`). Equidistant bins go to the earlier parent. Drop connected supports smaller than `minpix`, then close gaps of at most `close` valid energy bins inside each parent. Return one fine-grid mask per parent, in input order. Supply a `RebinResult` as `rebin` to map coarse masks to the fine grid first.
"""
function refine_supports(I::AbstractMatrix{<:Real}, σ::AbstractMatrix{<:Real}, valid::AbstractMatrix{Bool}, qs, ωs, coarse_masks;
                         resolution, z::Real=3.0, fraction::Real=0.15, bridge=nothing, close::Integer=2, minpix::Integer=10,
                         rebin=nothing)
    size(I) == size(σ) == size(valid) == (length(qs), length(ωs)) || throw(DimensionMismatch("I, σ, valid and axes must agree"))
    _ridge_axis(qs, "qs")
    _ridge_axis(ωs, "ωs")
    isfinite(z) && z >= 0 && 0 <= fraction <= 1 && close >= 0 && minpix >= 1 ||
        throw(ArgumentError("invalid significance, fraction, gap width or minimum support size"))
    Γq, Γω = resolution
    gap = bridge === nothing ? max(mean(diff(ωs)), Γω / 2) : bridge
    isfinite(gap) && gap >= 0 || throw(ArgumentError("bridge must be finite and nonnegative"))
    coarse_shape = rebin === nothing ? size(I) : size(rebin.I)
    masks = _ridge_parent_masks(coarse_masks, coarse_shape)
    seeds = [_checked_ridge_mask(rebin === nothing ? mask : upsample_mask(mask, rebin), size(I)) for mask in masks]
    isempty(seeds) && return BitMatrix[]
    measured = valid .& isfinite.(I) .& isfinite.(σ)
    any(i -> measured[i] && σ[i] < 0, eachindex(σ)) && throw(ArgumentError("measured uncertainties must be nonnegative"))
    S, error = _smooth_ridge_measurements(I, σ, measured, _ridge_kernel(qs, Γq), _ridge_kernel(ωs, Γω))
    eligible = _ridge_eligible(S, error, z, fraction)
    adjacent = [searchsortedfirst(ωs, ω - gap):searchsortedlast(ωs, ω + gap) for ω in ωs]
    owner = _assign_ridge_supports(eligible, seeds, adjacent, minpix)
    _close_ridge_supports!(owner, measured, close)
    return [owner .== parent for parent in eachindex(seeds)]
end

"""
    mask_dispersions(I, masks, qs, ωs; centre=:weighted, valid=nothing)

Return `(q, ω, imputed)` for arbitrary parent masks. Rows of `ω` correspond to masks. Use positive measured intensity as weights, or `centre=:bin` for the strongest measured bin. With no positive weight, use the strongest measured bin and flag the estimate. Interpolate missing internal rows in physical momentum and flag them. Leading and trailing gaps remain `NaN`.
"""
function mask_dispersions(I::AbstractMatrix{<:Real}, masks, qs, ωs; centre::Symbol=:weighted, valid=nothing)
    size(I) == (length(qs), length(ωs)) || throw(DimensionMismatch("axes must match I"))
    _ridge_axis(qs, "qs")
    _ridge_axis(ωs, "ωs")
    centre in (:weighted, :bin) || throw(ArgumentError("centre must be :weighted or :bin"))
    supports = _ridge_parent_masks(masks, size(I))
    measured = _ridge_valid(I, valid)
    ω = fill(NaN, length(supports), length(qs))
    imputed = falses(size(ω))
    for (parent, mask) in pairs(supports)
        _checked_ridge_mask(mask, size(I))
        _mask_trajectory!(@view(ω[parent, :]), @view(imputed[parent, :]), I, measured, mask, qs, ωs, centre)
    end
    return (q=Float64.(qs), ω=ω, imputed=imputed)
end

function _mask_trajectory!(ω, imputed, I, measured, mask, qs, ωs, centre)
    energy, weights = zeros(length(qs)), zeros(length(qs))
    peak, best = zeros(Int, length(qs)), fill(-Inf, length(qs))
    for j in eachindex(ωs), i in eachindex(qs)
        mask[i, j] && measured[i, j] || continue
        x = I[i, j]
        if x > best[i]
            best[i], peak[i] = x, j
        end
        w = max(x, 0)
        energy[i] += w * ωs[j]
        weights[i] += w
    end
    left = 0
    for i in eachindex(qs)
        if peak[i] == 0
            imputed[i] = true
            continue
        end
        weighted = centre === :weighted && weights[i] > 0
        ω[i] = weighted ? energy[i] / weights[i] : ωs[peak[i]]
        imputed[i] = centre === :weighted && !weighted
        if left > 0 && i > left + 1
            for k in left+1:i-1
                t = (qs[k] - qs[left]) / (qs[i] - qs[left])
                ω[k] = (1-t) * ω[left] + t * ω[i]
            end
        end
        left = i
    end
end
