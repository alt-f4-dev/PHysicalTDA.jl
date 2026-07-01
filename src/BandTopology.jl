#--------------------------------#
# Domain-Level Branch Topology   #
#--------------------------------#

"""
    restrict_to_domain(Iqω, domains, domain; pad_ω=0)

Return a copy of `Iqω` restricted to one energy-domain index.

`domains` is usually `result.qfirst.domains` from `adaptive_autosplitbands2d`.
The returned matrix has the same size as `Iqω`; pixels outside the requested
domain are set to zero.

This is a geometry operation only. It does not assume that the domain is
acoustic, optical, magnetic, phononic, or otherwise physically labeled.
"""
function restrict_to_domain(Iqω::AbstractMatrix{<:Real},
                            domains,
                            domain::Integer;
                            pad_ω::Int=0)
    domain < 1 && throw(ArgumentError("domain index must be positive"))

    nq, nω = size(Iqω)
    A = zeros(Float64, nq, nω)

    @inbounds for iq in 1:nq
        domain <= length(domains[iq]) || continue

        r = domains[iq][domain]
        isempty(r) && continue

        lo = max(1, first(r) - pad_ω)
        hi = min(nω, last(r) + pad_ω)

        A[iq, lo:hi] .= Iqω[iq, lo:hi]
    end

    return A
end


"""
    robust_branch_count(peak_counts; quantile_level=0.85, max_branches=nothing)

Estimate the resolved ridge/branch count from local peak multiplicities.

Unlike the split-merge slide helper, this function does not cap the count at 2
unless `max_branches=2` is explicitly passed. This makes it usable for more
general cases, including multiple branches, crossings, avoided crossings,
and complicated parent domains.

Returns `0` if no peaks are detected.
"""
function robust_branch_count(peak_counts::AbstractVector{<:Integer};
                             quantile_level::Float64=0.85,
                             max_branches::Union{Nothing,Int}=nothing)
    nonzero = [Int(c) for c in peak_counts if c > 0]
    isempty(nonzero) && return 0

    sort!(nonzero)
    idx = clamp(ceil(Int, quantile_level * length(nonzero)), 1, length(nonzero))
    count_est = nonzero[idx]

    if max_branches !== nothing
        count_est = min(count_est, max_branches)
    end

    return max(count_est, 0)
end


"""
    has_split_merge_pattern(peak_counts; low_count=1, high_count=2)

Return `true` when the peak-count sequence has a low → high → low pattern.

The default detects the usual one-branch → two-branch → one-branch signature,
but the function is written generically enough to support other multiplicity
thresholds.

This is a multiplicity diagnostic, not a complete physical classifier.
A true crossing, avoided crossing, diffuse continuum, or intensity exchange
may require additional domain-interaction diagnostics.
"""
function has_split_merge_pattern(peak_counts::AbstractVector{<:Integer};
                                 low_count::Int=1,
                                 high_count::Int=2)
    isempty(peak_counts) && return false

    high_idxs = findall(c -> c >= high_count, peak_counts)
    isempty(high_idxs) && return false

    qlo = first(high_idxs)
    qhi = last(high_idxs)

    left_has_low = qlo > 1 && any(c -> c <= low_count, @view peak_counts[1:qlo-1])
    right_has_low = qhi < length(peak_counts) &&
                    any(c -> c <= low_count, @view peak_counts[qhi+1:end])

    return left_has_low && right_has_low
end


"""
    analyze_domain_topologies(Iqω, qs, ωs, result; kwargs...)

Run ridge-topology analysis separately inside each inferred parent/domain.

This generalizes the family-aware slide-script logic by removing acoustic/
optical labels. The output is domain-indexed, not family-labeled.

Returns a named tuple:

    (
        topologies = Vector{RidgeTopologyResult},
        domain_images = Vector{Matrix{Float64}},
        summary = summarize_domain_branches(topologies)
    )

Important distinction:

    adaptive_autosplitbands2d(...; ridge_topology=true)

runs one global ridge topology analysis on the full image, while

    analyze_domain_topologies(...)

runs one topology analysis per inferred domain.
"""
function analyze_domain_topologies(Iqω::AbstractMatrix{<:Real},
                                   qs,
                                   ωs,
                                   result;
                                   domains = result.qfirst.domains,
                                   domain_count::Union{Nothing,Int}=nothing,
                                   pad_ω::Union{Nothing,Int}=nothing,
                                   smooth_radius_ω::Union{Nothing,Int}=nothing,
                                   include_energy_edges::Bool=false,
                                   min_prominence::Union{Nothing,Real}=nothing,
                                   max_peaks_per_slice::Int=3,
                                   max_jump_bins=nothing,
                                   tube_radius_ω::Int=2,
                                   periodic_q::Bool=false,
                                   pixel_adjacency::Symbol=:full,
                                   relative_score_floor::Float64=0.10,
                                   min_branch_separation_bins=nothing,
                                   max_branch_separation_frac::Float64=0.32,
                                   min_cycle_q_fraction::Float64=0.08,
                                   min_two_peak_fraction::Float64=0.15,
                                   min_cycle_confidence::Float64=0.55,
                                   branch_quantile_level::Float64=0.85,
                                   max_branches::Union{Nothing,Int}=nothing)

    ndomains = domain_count === nothing ? result.qfirst.nparents : Int(domain_count)

    pad = pad_ω === nothing ? max(1, Int(result.radius_ω)) : max(0, Int(pad_ω))
    rω = smooth_radius_ω === nothing ? max(1, Int(result.peak_radius_ω)) :
                                      max(1, Int(smooth_radius_ω))

    topologies = Vector{RidgeTopologyResult}(undef, ndomains)
    domain_images = Vector{Matrix{Float64}}(undef, ndomains)

    for b in 1:ndomains
        Idomain = restrict_to_domain(Iqω, domains, b; pad_ω=pad)
        domain_images[b] = Idomain

        topologies[b] = build_ridge_topology(
            Idomain, qs, ωs;
            smooth_radius_ω=rω,
            include_energy_edges=include_energy_edges,
            min_prominence=min_prominence,
            max_peaks_per_slice=max_peaks_per_slice,
            max_jump_bins=max_jump_bins,
            tube_radius_ω=tube_radius_ω,
            periodic_q=periodic_q,
            pixel_adjacency=pixel_adjacency,
            relative_score_floor=relative_score_floor,
            min_branch_separation_bins=min_branch_separation_bins,
            max_branch_separation_frac=max_branch_separation_frac,
            min_cycle_q_fraction=min_cycle_q_fraction,
            min_two_peak_fraction=min_two_peak_fraction,
            min_cycle_confidence=min_cycle_confidence,
        )
    end

    summary = summarize_domain_branches(
        topologies;
        quantile_level=branch_quantile_level,
        max_branches=max_branches,
    )

    return (
        topologies = topologies,
        domain_images = domain_images,
        summary = summary,
    )
end


"""
    summarize_domain_branches(topologies; quantile_level=0.85, max_branches=nothing)

Summarize one `RidgeTopologyResult` per domain.

Returns a named tuple:

    (
        branch_counts,
        beta0,
        beta1,
        split_merge,
        confidence,
        status,
        total_branches
    )

Statuses are heuristic diagnostic labels:

    :empty        no peak-count data
    :ridge_poor   too few q-slices contain resolved ridges
    :single       one resolved branch
    :multiple     multiple resolved branches
    :branching    split/merge or cycle-like topology detected

The status should be interpreted as a topology diagnostic, not as a final
physical label.
"""
function summarize_domain_branches(topologies::AbstractVector;
                                   quantile_level::Float64=0.85,
                                   max_branches::Union{Nothing,Int}=nothing,
                                   min_resolved_fraction::Float64=0.20)

    ndomains = length(topologies)

    branch_counts = zeros(Int, ndomains)
    beta0 = zeros(Int, ndomains)
    beta1 = zeros(Int, ndomains)
    split_merge = falses(ndomains)
    confidence = zeros(Float64, ndomains)
    status = Vector{Symbol}(undef, ndomains)

    for b in 1:ndomains
        topo = topologies[b]
        peak_counts = topo.peak_counts

        branch_counts[b] = robust_branch_count(
            peak_counts;
            quantile_level=quantile_level,
            max_branches=max_branches,
        )

        beta0[b] = topo.beta0
        beta1[b] = topo.beta1

        split_merge[b] = topo.beta1 > 0 || has_split_merge_pattern(peak_counts)

        if isempty(peak_counts)
            confidence[b] = 0.0
            status[b] = :empty
            continue
        end

        resolved_fraction = count(c -> c > 0, peak_counts) / length(peak_counts)
        modal_support = branch_counts[b] == 0 ? 0.0 :
            count(c -> c >= branch_counts[b], peak_counts) / length(peak_counts)

        confidence[b] = clamp(0.5 * resolved_fraction + 0.5 * modal_support, 0.0, 1.0)

        status[b] =
            resolved_fraction == 0.0 ? :empty :
            resolved_fraction < min_resolved_fraction ? :ridge_poor :
            split_merge[b] ? :branching :
            branch_counts[b] <= 1 ? :single :
            :multiple
    end

    return (
        branch_counts = branch_counts,
        beta0 = beta0,
        beta1 = beta1,
        split_merge = split_merge,
        confidence = confidence,
        status = status,
        total_branches = sum(branch_counts),
    )
end


"""
    domain_marginals(Iqω, result, qs, ωs; domain_count=nothing)

Compute domain-level marginal spectra.

This differs from `components_to_band_marginals(result)`: component marginals
are one row per retained connected support component, while domain marginals
are one row per inferred parent/domain.

Returns:

    (
        domainIq = Matrix{Float64},  # ndomains × nq
        domainIω = Matrix{Float64},  # ndomains × nω
    )
"""
function domain_marginals(Iqω::AbstractMatrix{<:Real},
                          result,
                          qs,
                          ωs;
                          domain_count::Union{Nothing,Int}=nothing)

    ndomains = domain_count === nothing ? result.qfirst.nparents : Int(domain_count)

    return domain_marginals(
        Iqω,
        result.qfirst.domains,
        qs,
        ωs;
        domain_count=ndomains,
    )
end


"""
    domain_marginals(Iqω, domains, qs, ωs; domain_count=nothing)

Compute marginal spectra directly from a supplied domain list.

Use this method when the domains do not come directly from an
`adaptive_autosplitbands2d` result.
"""
function domain_marginals(Iqω::AbstractMatrix{<:Real},
                          domains::AbstractVector,
                          qs,
                          ωs;
                          domain_count::Union{Nothing,Int}=nothing)

    nq, nω = size(Iqω)

    ndomains = if domain_count === nothing
        maximum(length(d) for d in domains)
    else
        Int(domain_count)
    end

    Δq = mean(diff(qs))
    Δω = mean(diff(ωs))

    domainIq = zeros(Float64, ndomains, nq)
    domainIω = zeros(Float64, ndomains, nω)

    @inbounds for b in 1:ndomains
        for iq in 1:nq
            b <= length(domains[iq]) || continue

            r = domains[iq][b]
            isempty(r) && continue

            for iω in r
                val = Float64(Iqω[iq, iω])
                domainIq[b, iq] += val * Δω
                domainIω[b, iω] += val * Δq
            end
        end
    end

    return (
        domainIq = domainIq,
        domainIω = domainIω,
    )
end
