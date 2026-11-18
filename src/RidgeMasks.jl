function _checked_ridge_mask(mask::AbstractMatrix{Bool}, shape)
    size(mask) == shape || throw(DimensionMismatch("mask must have size $shape"))
    return mask
end

function _ridge_valid(I, valid)
    valid === nothing || _checked_ridge_mask(valid, size(I))
    return (valid === nothing ? isfinite.(I) : valid .& isfinite.(I))
end

_ridge_parent_masks(mask::AbstractMatrix{Bool}, shape) = (mask,)
_ridge_parent_masks(masks::Union{AbstractVector,Tuple}, shape) = masks
_ridge_parent_masks(result, shape) = band_supports(result, shape...)

function _check_disjoint_ridge_masks(masks, shape)
    occupied = falses(shape)
    for mask in masks, i in eachindex(mask)
        mask[i] || continue
        occupied[i] && throw(ArgumentError("physical parent masks must be disjoint on measured bins"))
        occupied[i] = true
    end
end

function _nearest_ridge_bin(axis, x)
    j = searchsortedfirst(axis, x)
    j <= 1 && return 1
    j > length(axis) && return length(axis)
    return x - axis[j-1] <= axis[j] - x ? j-1 : j
end

function _check_masked_tracking_params(p)
    scales = (p.prediction_scale, p.max_slope)
    weights = (p.curvature_weight, p.score_weight, p.fallback_shape_weight)
    all(x -> isfinite(x) && x > 0, scales) && all(isfinite, weights) &&
        p.max_children >= 1 && p.extra_assignment_candidates >= 0 && p.max_peaks_per_slice >= 1 ||
        throw(ArgumentError("invalid tracking parameters"))
end

# Scan measured runs separately so that an invalid bin never supplies a peak's neighbours or baseline.
function _masked_ridge_candidates(I, q, ω, p, valid, domains; σ=nothing)
    _ridge_axis(q, "qs")
    _ridge_axis(ω, "ωs")
    p.peak_smooth_radius >= 0 && isfinite(p.fallback_window) && p.fallback_window > 0 ||
        throw(ArgumentError("masked recovery requires a nonnegative smoothing radius and a positive finite fallback window"))
    _check_masked_tracking_params(p)
    nq, nω = size(I)
    A = zeros(Float32, nq, nω)
    nodes, groups = RidgeNode[], RidgeSupportGroup[]
    byq = [Int32[] for _ in q]
    groups_by_q = [Int32[] for _ in q]
    parent_groups = [Int32[] for _ in domains]
    prefix, work = zeros(nω + 1), zeros(nω)
    scores = Tuple{Int,Float64,Bool}[]
    kept = Int[]
    for (parent, domain) in pairs(domains), iq in 1:nq
        empty!(scores)
        j = 1
        while j <= nω
            if !domain[iq, j]
                j += 1
                continue
            end
            lo = j
            while j < nω && domain[iq, j+1]
                j += 1
            end
            _masked_peak_run!(scores, work, prefix, A, I, valid, iq, lo, j, p, σ)
            j += 1
        end
        isempty(scores) && continue
        sort!(scores; by=x -> (-x[2], x[1]))
        floor = p.relative_score_floor * first(scores)[2]
        empty!(kept)
        for (bin, score, significant) in scores
            length(kept) < p.max_peaks_per_slice || break
            score >= floor || significant || continue
            any(k -> abs(k - bin) < max(2, p.peak_smooth_radius) &&
                     all(@view domain[iq, min(k, bin):max(k, bin)]), kept) && continue
            push!(kept, bin)
        end
        sort!(kept)
        for bin in kept
            push!(nodes, RidgeNode(iq, bin, Float64(q[iq]), Float64(ω[bin]), Float64(A[iq, bin])))
            id = Int32(length(nodes))
            push!(byq[iq], id)
            threshold = max(p.support_min_level, p.support_fraction * A[iq, bin])
            lo = hi = bin
            while lo > 1 && domain[iq, lo-1] && ω[bin] - ω[lo-1] <= p.support_max_half_width && A[iq, lo-1] >= threshold
                lo -= 1
            end
            while hi < nω && domain[iq, hi+1] && ω[hi+1] - ω[bin] <= p.support_max_half_width && A[iq, hi+1] >= threshold
                hi += 1
            end
            gs = parent_groups[parent]
            if !isempty(gs) && groups[last(gs)].iq == iq && lo <= groups[last(gs)].hi
                g = groups[last(gs)]
                push!(g.candidates, id)
                groups[last(gs)] = RidgeSupportGroup(g.iq, min(g.lo, Int32(lo)), max(g.hi, Int32(hi)), g.candidates)
            else
                push!(groups, RidgeSupportGroup(Int32(iq), Int32(lo), Int32(hi), Int32[id]))
                push!(gs, Int32(length(groups)))
                push!(groups_by_q[iq], Int32(length(groups)))
            end
        end
    end
    return RidgeCandidates(I, A, q, ω, nodes, byq, groups, groups_by_q, p), parent_groups
end

function _masked_peak_run!(scores, work, prefix, A, I, valid, iq, lo, hi, p, σ=nothing)
    radius = p.peak_smooth_radius
    prefix[lo] = 0.0
    for j in lo:hi
        prefix[j+1] = prefix[j] + I[iq, j]
    end
    for j in lo:hi
        a, b = max(lo, j-radius), min(hi, j+radius)
        work[j] = (prefix[b+1] - prefix[a]) / (b-a+1)
    end
    baseline, top = extrema(@view work[lo:hi])
    nω = size(I, 2)
    # At a parent boundary, a measured neighbouring bin can establish a maximum without contributing to smoothing.
    left = lo > 1 && valid[iq, lo-1] ? Float64(I[iq, lo-1]) : Inf
    right = hi < nω && valid[iq, hi+1] ? Float64(I[iq, hi+1]) : Inf
    lo == 1 && p.include_energy_edges && (left = -Inf)
    hi == nω && p.include_energy_edges && (right = -Inf)
    isfinite(left) && (baseline = min(baseline, left))
    isfinite(right) && (baseline = min(baseline, right))
    top > baseline || return
    for j in lo:hi
        x = work[j]
        a, b = j == lo ? left : work[j-1], j == hi ? right : work[j+1]
        score = (x - baseline) / (top - baseline)
        A[iq, j] = score
        x > a && x >= b || continue
        lmin = minimum(@view work[max(lo, j-4radius-2):j])
        rmin = minimum(@view work[j:min(hi, j+4radius+2)])
        j == lo && (lmin = isfinite(left) ? left : baseline)
        j == hi && (rmin = isfinite(right) ? right : baseline)
        prominence = (x - max(lmin, rmin)) / (top - baseline)
        significant = false
        if σ !== nothing && isfinite(p.min_significance)
            a2, b2 = max(lo, j-radius), min(hi, j+radius)
            σs = sqrt(sum(abs2, @view σ[iq, a2:b2])) / (b2 - a2 + 1)
            significant = σs > 0 && (x - max(lmin, rmin)) / σs >= p.min_significance
        end
        prominence > 0 && (prominence >= p.min_prominence || significant) && push!(scores, (j, score, significant))
    end
end
