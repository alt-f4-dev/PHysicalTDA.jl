#----------------------------#
# Physical Band Segmentation #
#----------------------------#----------------------------------------------#
# Resolution-aware 2D band segmentation for I(q,ω) matrices.                #
# Returns parent-band supports, component marginals, and ridge diagnostics. #
#                                                                           #
# Dependencies: Statistics and Base.Threads.                                #
#---------------------------------------------------------------------------#

using Statistics
using Base.Threads

#---------------------------------#
# Robust Statistics and Smoothing #
#---------------------------------#
@inline periodic_index(i::Integer, n::Integer) = mod1(i, n)
function robust_mad(x)
    μ = median(x)
    return 1.4826 * median(abs.(x .- μ))
end
@inline robust_sigma(x) = max(robust_mad(x), eps(Float64))

function robust_background_sigma(A::AbstractMatrix{<:Real})
    vals = vec(Float64.(A))
    μ = median(vals)
    σ = robust_sigma(vals)
    return μ, σ
end

function normalize01!(y::AbstractVector{Float64})
    ymin, ymax = extrema(y)
    if ymax == ymin
        fill!(y, 0.0)
    else
        invrng = 1.0 / (ymax - ymin)
        @inbounds @simd for i in eachindex(y)
            y[i] = (y[i] - ymin) * invrng
        end
    end
    return y
end

function moving_average_1d!(dest::Vector{Float64}, src::AbstractVector{<:Real}, radius::Int,
                            prefix::Union{Nothing,Vector{Float64}}=nothing)
    # O(n) box smoothing by prefix sums. Energy is not periodic.
    n = length(src)
    if radius <= 0
        @inbounds @simd for i in 1:n
            dest[i] = Float64(src[i])
        end
        return dest
    end
    p = prefix === nothing ? Vector{Float64}(undef, n + 1) : prefix
    @assert length(p) >= n + 1 "prefix buffer too short"
    p[1] = 0.0
    @inbounds for i in 1:n
        p[i + 1] = p[i] + Float64(src[i])
    end
    @inbounds @simd for i in 1:n
        lo = max(1, i - radius)
        hi = min(n, i + radius)
        dest[i] = (p[hi + 1] - p[lo]) / (hi - lo + 1)
    end
    return dest
end

function separable_smooth2d!(B::Matrix{Float64}, A::AbstractMatrix{<:Real};
                             radius_q::Int=1, radius_ω::Int=2,
                             periodic_q::Bool=true)
    # Separable smoothing. Thread-local prefix buffers avoid per-row allocation.
    nq, nω = size(A)
    tmp = Matrix{Float64}(undef, nq, nω)

    @threads for iω in 1:nω
        if radius_q >= nq ÷ 2
            acc = 0.0
            @inbounds for iq in 1:nq
                acc += A[iq, iω]
            end
            μ = acc / nq
            @inbounds @simd for iq in 1:nq
                tmp[iq, iω] = μ
            end
        else
            @inbounds for iq in 1:nq
                acc = 0.0
                cnt = 0
                for dq in -radius_q:radius_q
                    qj_raw = iq + dq
                    if periodic_q
                        qj = periodic_index(qj_raw, nq)
                        acc += A[qj, iω]
                        cnt += 1
                    elseif 1 <= qj_raw <= nq
                        acc += A[qj_raw, iω]
                        cnt += 1
                    end
                end
                tmp[iq, iω] = acc / cnt
            end
        end
    end

    prefixes = [Vector{Float64}(undef, nω + 1) for _ in 1:Threads.maxthreadid()]
    @threads for iq in 1:nq
        prefix = prefixes[Threads.threadid()]
        prefix[1] = 0.0
        @inbounds for iω in 1:nω
            prefix[iω + 1] = prefix[iω] + tmp[iq, iω]
        end
        @inbounds @simd for iω in 1:nω
            lo = max(1, iω - radius_ω)
            hi = min(nω, iω + radius_ω)
            B[iq, iω] = (prefix[hi + 1] - prefix[lo]) / (hi - lo + 1)
        end
    end
    return B
end


#-------------------------#
# Adaptive Peak Detection #
#-------------------------#
function estimate_energy_resolution(Iqω::AbstractMatrix{<:Real}, ωs;
                                    nslices::Int=48,
                                    frac::Float64=0.5)
    # Estimate a characteristic energy width from sparse half-height widths.
    # The scale is approximately invariant under grid refinement.
    nq, nω = size(Iqω)
    idxs = round.(Int, range(1, nq; length=min(nq, nslices)))
    widths = Float64[]
    @inbounds for iq in idxs
        spec = @view Iqω[iq, :]
        imax = argmax(spec)
        peak = Float64(spec[imax])
        base = median(Float64.(spec))
        amp = peak - base
        amp <= 0 && continue
        level = base + frac * amp
        lo = imax
        while lo > 1 && spec[lo] > level
            lo -= 1
        end
        hi = imax
        while hi < nω && spec[hi] > level
            hi += 1
        end
        hi > lo && push!(widths, abs(ωs[hi] - ωs[lo]))
    end
    isempty(widths) && return 2 * mean(diff(ωs))
    return max(median(widths), 2 * mean(diff(ωs)))
end

function adaptive_resolution_radii(qs, ωs, Iqω=nothing;
                                    resolution_q=nothing,
                                    resolution_ω=nothing)
    Δq = mean(diff(qs))
    Δω = mean(diff(ωs))

    # Supplied instrumental resolution is authoritative.
    # Otherwise infer an energy width and use a conservative q radius.
    Γω = isnothing(resolution_ω) ?
        (Iqω === nothing ? 2 * Δω : estimate_energy_resolution(Iqω, ωs)) :
        Float64(resolution_ω)
    # Default q resolution is physical-scale, not purely pixel-scale.
    # Pass instrument-derived resolution when it is available.
    qspan = maximum(qs) - minimum(qs)
    Γq = isnothing(resolution_q) ? max(2 * Δq, 0.01 * qspan) : Float64(resolution_q)

    rq = max(1, ceil(Int, 0.5 * Γq / Δq))
    rω = max(1, ceil(Int, 0.5 * Γω / Δω))
    return rq, rω, Γq, Γω
end

function adaptive_prominence(s::AbstractVector{Float64}; min_floor::Float64=0.025,
                             kσ::Float64=3.5)
    # Resolution-invariant default: estimate pointwise noise from first
    # differences of the already resolution-smoothed, normalized slice.
    n = length(s)
    n < 3 && return min_floor
    diffs = Vector{Float64}(undef, n-1)
    @inbounds for i in 1:n-1
        diffs[i] = s[i+1] - s[i]
    end
    σ = robust_sigma(diffs) / sqrt(2)
    return max(min_floor, kσ*σ)
end

function peak_candidates!(peaks::Vector{Int}, scores::Vector{Float64},
                          work::Vector{Float64}, spec::AbstractVector{<:Real}, ωs;
                          smooth_radius::Int=2,
                          include_energy_edges::Bool=true,
                          min_prominence::Union{Nothing,Real}=nothing,
                          min_separation_bins::Int=2,
                          prefix::Union{Nothing,Vector{Float64}}=nothing)
    empty!(peaks); empty!(scores)
    moving_average_1d!(work, spec, smooth_radius, prefix)
    normalize01!(work)
    n = length(work)
    prom_cut = isnothing(min_prominence) ? adaptive_prominence(work) : Float64(min_prominence)

    # Energy is not periodic, but edge peaks are valid measured-window extrema.
    if include_energy_edges && n >= 2 && work[1] >= work[2]
        prom = work[1] - minimum(work)
        prom >= prom_cut && push!(peaks, 1)
    end
    edge_window = min(max(2*smooth_radius+2, min_separation_bins +2), max(2, n ÷ 5))
    @inbounds for j in 2:n-1
        if work[j] > work[j-1] && work[j] >= work[j+1]
            # Local prominence using nearest valley estimate. This is cheaper
            # and less overzealous than full left/right global minima scans.
            lo = max(1, j-4*smooth_radius-2)
            hi = min(n, j+4*smooth_radius+2)
            left_min = minimum(@view work[lo:j])
            right_min = minimum(@view work[j:hi])
            prom = if include_energy_edges && j <= edge_window
                work[j] - right_min
            elseif include_energy_edges && j >= n - edge_window + 1
                work[j] - left_min
            else
                work[j] - max(left_min, right_min)
            end
            prom >= prom_cut && push!(peaks, j)
        end
    end

    if include_energy_edges && n >= 2 && work[n] >= work[n-1]
        prom = work[n] - minimum(work)
        prom >= prom_cut && push!(peaks, n)
    end

    isempty(peaks) && return peaks, scores, prom_cut

    # Non-maximum suppression by score.
    sort!(peaks; by=j -> -work[j])
    kept = Int[]
    @inbounds for p in peaks
        ok = true
        for k in kept
            if abs(p-k) < min_separation_bins
                ok = false
                break
            end
        end
        ok && push!(kept, p)
    end
    sort!(kept)
    empty!(peaks)
    append!(peaks, kept)
    for p in peaks
        push!(scores, work[p])
    end
    return peaks, scores, prom_cut
end

function infer_parent_count(peak_counts::Vector{Int}; method::Symbol=:quantile,
                            max_bands::Union{Nothing,Int}=nothing,
                            quantile_level::Float64=0.85)
    nonzero = [c for c in peak_counts if c > 0]
    isempty(nonzero) && return 1
    N = if method === :max
        maximum(nonzero)
    elseif method === :median
        max(1, round(Int, median(nonzero)))
    elseif method === :mode
        counts = Dict{Int,Int}()
        for c in nonzero
            counts[c] = get(counts, c, 0) + 1
        end
        best_c, best_n = first(counts)
        for (c,n) in counts
            if n > best_n || (n == best_n && c > best_c)
                best_c, best_n = c, n
            end
        end
        best_c
    else
        # High-quantile count is robust to missed weak branches.
        # It is safer than the maximum when rare noise creates extra peaks.
        sort!(nonzero)
        idx = clamp(ceil(Int, quantile_level * length(nonzero)), 1, length(nonzero))
        nonzero[idx]
    end
    max_bands !== nothing && (N = min(N, max_bands))
    return max(N, 1)
end

function fill_missing_centers!(centers::Matrix{Int}, imputed::BitMatrix;
                               periodic_q::Bool=true)
    # O(nparents*nq) nearest-neighbor imputation along q.
    # q may be periodic; energy is not.
    nparents, nq = size(centers)
    @inbounds for b in 1:nparents
        valid = Int[]
        for iq in 1:nq
            centers[b, iq] != 0 && push!(valid, iq)
        end
        if isempty(valid)
            for iq in 1:nq
                centers[b, iq] = 1
                imputed[b, iq] = true
            end
            continue
        end
        if length(valid) == nq
            continue
        end

        if periodic_q
            extpos = Vector{Int}(undef, 3*length(valid))
            extval = Vector{Int}(undef, 3*length(valid))
            k = 1
            for offset in (-nq, 0, nq)
                for v in valid
                    extpos[k] = v + offset
                    extval[k] = centers[b, v]
                    k += 1
                end
            end
            p = 1
            for iq in 1:nq
                centers[b, iq] != 0 && continue
                while p < length(extpos) && abs(extpos[p+1] - iq) <= abs(extpos[p] - iq)
                    p += 1
                end
                centers[b, iq] = extval[p]
                imputed[b, iq] = true
            end
        else
            left_pos = fill(0, nq); left_val = fill(0, nq)
            right_pos = fill(0, nq); right_val = fill(0, nq)
            lastp = 0; lastv = 0
            for iq in 1:nq
                if centers[b, iq] != 0
                    lastp = iq; lastv = centers[b, iq]
                end
                left_pos[iq] = lastp; left_val[iq] = lastv
            end
            lastp = 0; lastv = 0
            for iq in nq:-1:1
                if centers[b, iq] != 0
                    lastp = iq; lastv = centers[b, iq]
                end
                right_pos[iq] = lastp; right_val[iq] = lastv
            end
            for iq in 1:nq
                centers[b, iq] != 0 && continue
                lp = left_pos[iq]; rp = right_pos[iq]
                if lp == 0
                    centers[b, iq] = right_val[iq]
                elseif rp == 0
                    centers[b, iq] = left_val[iq]
                else
                    centers[b, iq] = (iq-lp <= rp-iq) ? left_val[iq] : right_val[iq]
                end
                imputed[b, iq] = true
            end
        end
    end
    return centers, imputed
end

function candidate_reliability(score::Float64, prom_cut::Float64)
    # Dimensionless candidate reliability.
    # Approximates persistence-weighted confidence without per-slice PH.
    return score / max(prom_cut, eps(Float64))
end

function nearest_full_center_column(centers::Matrix{Int}, iq::Int; periodic_q::Bool=true)
    nparents, nq = size(centers)
    best_j = 0
    best_d = typemax(Int)
    @inbounds for j in 1:nq
        j == iq && continue
        full = true
        for b in 1:nparents
            if centers[b,j] == 0
                full = false
                break
            end
        end
        full || continue
        d = abs(j - iq)
        periodic_q && (d = min(d, nq-d))
        if d < best_d
            best_d = d
            best_j = j
        end
    end
    return best_j
end

function assign_peaks_to_reference!(centers::Matrix{Int},iq::Int, ps::Vector{Int}, refs::AbstractVector{<:Integer})
    nparents = size(centers,1)
    used_parent = falses(nparents)
    used_peak = falses(length(ps))
    pairs = Tuple{Int,Int,Int}[]
    sizehint!(pairs, nparents*length(ps))
    @inbounds for b in 1:nparents
        refs[b] == 0 && continue
        for k in eachindex(ps)
            push!(pairs, (abs(ps[k] - Int(refs[b])), b, k))
        end
    end
    sort!(pairs; by = x -> x[1])
    @inbounds for pair in pairs
        b = pair[2]
        k = pair[3]
        if !used_parent[b] && !used_peak[k]
            centers[b,iq] = ps[k]
            used_parent[b] = true
            used_peak[k] = true
        end
    end
    return centers
end


function qfirst_parent_domains(Iqω::AbstractMatrix{<:Real}, qs, ωs;
                               smooth_radius_ω::Int=2,
                               include_energy_edges::Bool=true,
                               periodic_q::Bool=true,
                               min_prominence::Union{Nothing,Real}=nothing,
                               parent_count_method::Symbol=:quantile,
                               max_bands::Union{Nothing,Int}=nothing,
                               verbose::Bool=true)
    nq, nω = size(Iqω)
    works = [Vector{Float64}(undef, nω) for _ in 1:Threads.maxthreadid()]
    prefixes = [Vector{Float64}(undef, nω + 1) for _ in 1:Threads.maxthreadid()]
    peaks_local = [Int[] for _ in 1:nq]
    scores_local = [Float64[] for _ in 1:nq]
    prom_cuts = zeros(Float64, nq)

    @threads for iq in 1:nq
        tid = Threads.threadid()
        work = works[tid]
        prefix = prefixes[tid]
        p = Int[]; sc = Float64[]
        peak_candidates!(p, sc, work, @view(Iqω[iq, :]), ωs;
                         smooth_radius=smooth_radius_ω,
                         include_energy_edges=include_energy_edges,
                         min_prominence=min_prominence,
                         min_separation_bins=max(2, smooth_radius_ω),
                         prefix=prefix)
        peaks_local[iq] = copy(p)
        scores_local[iq] = copy(sc)
        prom_cuts[iq] = adaptive_prominence(work)
    end

    peak_counts = length.(peaks_local)
    Nparents = infer_parent_count(peak_counts; method=parent_count_method, max_bands=max_bands)

    centers = zeros(Int, Nparents, nq)
    imputed = falses(Nparents, nq)

    # Parent labels are energy ordered.
    # Extra candidates are filtered by reliability before sorting.
    selected_peaks_by_q = [Int[] for _ in 1:nq]
    @inbounds for iq in 1:nq
        ps = peaks_local[iq]
        isempty(ps) && continue
        if length(ps) > Nparents
            ss = scores_local[iq]
            pc = prom_cuts[iq]
            rel = [candidate_reliability(Float64(x), pc) for x in ss]
            order = sortperm(rel; rev=true)[1:Nparents]
            ps = sort(ps[order])
        else
            ps = sort(ps)
        end
        selected_peaks_by_q[iq] = ps
        if length(ps) == Nparents
            for b in 1:Nparents
                centers[b,iq] = ps[b]
            end
        end
    end

    # Partial rows: match observed peaks to the nearest full row.
    @inbounds for iq in 1:nq
        ps = selected_peaks_by_q[iq]
        isempty(ps) && continue
        length(ps) == Nparents && continue
        ref_iq = nearest_full_center_column(centers, iq; periodic_q=periodic_q)
        if ref_iq == 0
            m = min(length(ps),Nparents)
            for b in 1:m
                centers[b,iq] = ps[b]
            end
        else
            refs = @view centers[:, ref_iq]
            assign_peaks_to_reference!(centers, iq, ps, refs)
        end
    end

    fill_missing_centers!(centers, imputed; periodic_q=periodic_q)

    # Voronoi energy domains from tracked centers.
    domains = Vector{Vector{UnitRange{Int}}}(undef, nq)
    @inbounds for iq in 1:nq
        c = sort(collect(@view centers[:, iq]))
        ranges = Vector{UnitRange{Int}}(undef, Nparents)
        start = 1
        for b in 1:Nparents-1
            bd = clamp((c[b] + c[b+1]) ÷ 2, 1, nω-1)
            ranges[b] = start:bd
            start = bd+1
        end
        ranges[Nparents] = start:nω
        domains[iq] = ranges
    end

    if verbose
        @info "q-first adaptive detection: boundary-aware peak counts=$(sort(unique(peak_counts)))"
        @info "q-first adaptive detection: inferred global parent-band count=$Nparents using method=$parent_count_method"
        @info "q-first adaptive detection: prominence cutoff median=$(round(median(prom_cuts); digits=5)), range=($(round(minimum(prom_cuts); digits=5)), $(round(maximum(prom_cuts); digits=5)))"
        @info "q-first adaptive detection: imputed centers total=$(count(imputed)), by parent=$([count(@view imputed[b, :]) for b in 1:Nparents])"
        @info "q-first adaptive detection: endpoint centers=$([(b, centers[b,1], centers[b,end]) for b in 1:Nparents])"
    end

    return (domains=domains, centers=centers, imputed=imputed,
            peaks=peaks_local, scores=scores_local, nparents=Nparents,
            peak_counts=peak_counts)
end

#--------------#
# Result Types #
#--------------#
struct ComponentResult
    parent_band::Int
    component::Int
    qbox::Tuple{Float64,Float64}
    ωbox::Tuple{Float64,Float64}
    qinds::UnitRange{Int}       # bounding-box q indices, not full membership
    ωinds::UnitRange{Int}       # bounding-box energy indices, not full membership
    npixels::Int
    weight::Float64
    bandIq::Vector{Float64}
    bandIω::Vector{Float64}
    linear_indices::Vector{Int} # exact component pixel membership, column-major
end



#----------------------#
# Ridge Topology Types #
#----------------------#
struct RidgeNode
    iq::Int
    iω::Int
    q::Float64
    ω::Float64
    score::Float64
end

struct RidgeCycleCandidate
    component::Int
    beta1::Int
    qspan::Tuple{Int,Int}
    q_fraction::Float64
    two_peak_fraction::Float64
    median_separation_bins::Float64
    median_score::Float64
    confidence::Float64
    accepted::Bool
    reason::String
end

struct RidgeTopologyResult
    nodes::Vector{RidgeNode}
    edges::Vector{Tuple{Int,Int}}
    beta0::Int
    beta1::Int
    peak_counts::Vector{Int}
    ridge_mask::BitMatrix
    ridge_labels::Matrix{Int32}
    ncomponents::Int
    cycle_candidates::Vector{RidgeCycleCandidate}
    accepted_cycle_count::Int
    status::Symbol
    diagnostics::Dict{Symbol,Any}
end

#-----------------#
# Graph Utilities #
#-----------------#
function _uf_find!(parent::Vector{Int}, x::Int)
    y = x
    @inbounds while parent[y] != y
        y = parent[y]
    end
    root = y
    y = x
    @inbounds while parent[y] != y
        nxt = parent[y]
        parent[y] = root
        y = nxt
    end
    return root
end

function _uf_union!(parent::Vector{Int}, rank::Vector{Int}, a::Int, b::Int)
    ra = _uf_find!(parent, a)
    rb = _uf_find!(parent, b)
    ra == rb && return false
    if rank[ra] < rank[rb]
        parent[ra] = rb
    elseif rank[ra] > rank[rb]
        parent[rb] = ra
    else
        parent[rb] = ra
        rank[ra] += 1
    end
    return true
end

function graph_betti_numbers(nnodes::Int, edges::Vector{Tuple{Int,Int}})
    nnodes == 0 && return 0, 0
    parent = collect(1:nnodes)
    rank = zeros(Int, nnodes)
    for (a, b) in edges
        1 <= a <= nnodes && 1 <= b <= nnodes || continue
        _uf_union!(parent, rank, a, b)
    end
    roots = Set{Int}()
    for i in 1:nnodes
        push!(roots, _uf_find!(parent, i))
    end
    beta0 = length(roots)
    beta1 = max(length(edges) - nnodes + beta0, 0)
    return beta0, beta1
end

function graph_component_labels(nnodes::Int, edges::Vector{Tuple{Int,Int}})
    nnodes == 0 && return Int[], 0
    parent = collect(1:nnodes)
    rank = zeros(Int, nnodes)
    for (a, b) in edges
        1 <= a <= nnodes && 1 <= b <= nnodes || continue
        _uf_union!(parent, rank, a, b)
    end
    root_to_label = Dict{Int,Int}()
    labels = zeros(Int, nnodes)
    next_label = 0
    for i in 1:nnodes
        r = _uf_find!(parent, i)
        if !haskey(root_to_label, r)
            next_label += 1
            root_to_label[r] = next_label
        end
        labels[i] = root_to_label[r]
    end
    return labels, next_label
end

function _mark_tube!(mask::BitMatrix, iq::Int, iω::Int, radius_ω::Int)
    nq, nω = size(mask)
    1 <= iq <= nq || return nothing
    lo = max(1, iω - radius_ω)
    hi = min(nω, iω + radius_ω)
    @inbounds for j in lo:hi
        mask[iq, j] = true
    end
    return nothing
end

function rasterize_ridge_graph(nodes::Vector{RidgeNode}, edges::Vector{Tuple{Int,Int}},
                               nq::Int, nω::Int; tube_radius_ω::Int=1,
                               periodic_q::Bool=false)
    mask = falses(nq, nω)
    for node in nodes
        _mark_tube!(mask, node.iq, node.iω, tube_radius_ω)
    end
    for (a, b) in edges
        na = nodes[a]
        nb = nodes[b]
        dq = nb.iq - na.iq
        if periodic_q && abs(dq) > 1
            continue
        end
        steps = max(abs(dq), 1)
        for t in 0:steps
            u = steps == 0 ? 0.0 : t / steps
            iq = round(Int, (1-u)*na.iq + u*nb.iq)
            iω = round(Int, (1-u)*na.iω + u*nb.iω)
            1 <= iq <= nq || continue
            _mark_tube!(mask, iq, iω, tube_radius_ω)
        end
    end
    return mask
end

#--------------------------#
# Ridge Graph Construction #
#--------------------------#
function ridge_peak_slices(Iqω::AbstractMatrix{<:Real}, ωs;
                           smooth_radius_ω::Int=2,
                           include_energy_edges::Bool=false,
                           min_prominence::Union{Nothing,Real}=nothing,
                           max_peaks_per_slice::Int=2,
                           relative_score_floor::Float64=0.1)
    nq, nω = size(Iqω)
    works = [Vector{Float64}(undef, nω) for _ in 1:Threads.maxthreadid()]
    prefixes = [Vector{Float64}(undef, nω + 1) for _ in 1:Threads.maxthreadid()]
    peaks_local = [Int[] for _ in 1:nq]
    scores_local = [Float64[] for _ in 1:nq]

    @threads for iq in 1:nq
        tid = Threads.threadid()
        p = Int[]
        sc = Float64[]
        peak_candidates!(p, sc, works[tid], @view(Iqω[iq, :]), ωs;
                         smooth_radius=smooth_radius_ω,
                         include_energy_edges=include_energy_edges,
                         min_prominence=min_prominence,
                         min_separation_bins=max(2, smooth_radius_ω),
                         prefix=prefixes[tid])
        if !isempty(p)
            smax = maximum(sc)
            keep = findall(x -> x >= relative_score_floor * smax, sc)
            p = p[keep]
            sc = sc[keep]
        end
        if length(p) > max_peaks_per_slice
            order = sortperm(sc; rev=true)[1:max_peaks_per_slice]
            p = p[order]
            sc = sc[order]
        end
        ordω = sortperm(p)
        peaks_local[iq] = copy(p[ordω])
        scores_local[iq] = copy(sc[ordω])
    end
    return peaks_local, scores_local
end

function _nearest_node_id(nodes::Vector{RidgeNode}, ids::Vector{Int}, target_iω::Int)
    @assert !isempty(ids)
    best = ids[1]
    bestdist = abs(nodes[best].iω - target_iω)
    @inbounds for k in 2:length(ids)
        id = ids[k]
        dist = abs(nodes[id].iω - target_iω)
        if dist < bestdist
            best = id
            bestdist = dist
        end
    end
    return best
end

@inline function _valid_splitmerge_geometry(single_iω::Int, low_iω::Int, high_iω::Int;
                                            min_sep::Int,
                                            max_sep::Int,
                                            midpoint_tol::Float64)
    low_iω < high_iω || return false
    sep = high_iω - low_iω
    min_sep <= sep <= max_sep || return false
    low_iω <= single_iω <= high_iω || return false
    mid = 0.5 * (low_iω + high_iω)
    return abs(single_iω - mid) <= max(midpoint_tol, 0.35 * sep)
end

function _add_ordered_matching_edges!(edges::Vector{Tuple{Int,Int}},
                                      nodes::Vector{RidgeNode},
                                      a_sorted::Vector{Int}, b_sorted::Vector{Int},
                                      jump::Int)
    m = min(length(a_sorted), length(b_sorted))
    @inbounds for k in 1:m
        a = a_sorted[k]
        b = b_sorted[k]
        abs(nodes[a].iω - nodes[b].iω) <= jump && push!(edges, (a, b))
    end
    return edges
end

function _component_edges(edges::Vector{Tuple{Int,Int}}, labels::Vector{Int}, comp::Int)
    return [(a,b) for (a,b) in edges if labels[a] == comp && labels[b] == comp]
end

function ridge_cycle_candidates(nodes::Vector{RidgeNode}, edges::Vector{Tuple{Int,Int}},
                                peak_counts::Vector{Int}, ncomponents::Int, labels::Vector{Int};
                                min_cycle_q_fraction::Float64=0.08,
                                min_two_peak_fraction::Float64=0.15,
                                min_median_score::Float64=0.10,
                                min_confidence::Float64=0.55)
    cands = RidgeCycleCandidate[]
    isempty(nodes) && return cands
    nq = length(peak_counts)
    for comp in 1:ncomponents
        ids = findall(==(comp), labels)
        isempty(ids) && continue
        comp_edges = _component_edges(edges, labels, comp)
        idmap = Dict{Int,Int}(id => k for (k, id) in pairs(ids))
        local_edges = Tuple{Int,Int}[]
        for (a, b) in comp_edges
            if haskey(idmap, a) && haskey(idmap, b)
                push!(local_edges, (idmap[a], idmap[b]))
            end
        end
        β0c, β1c = graph_betti_numbers(length(ids), local_edges)
        β1c == 0 && continue
        qvals = [nodes[id].iq for id in ids]
        qlo, qhi = minimum(qvals), maximum(qvals)
        qspan = max(qhi - qlo + 1, 1)
        qfrac = qspan / max(nq, 1)
        two_frac = count(i -> peak_counts[i] >= 2, qlo:qhi) / qspan
        byq = Dict{Int,Vector{Int}}()
        for id in ids
            push!(get!(byq, nodes[id].iq, Int[]), id)
        end
        seps = Float64[]
        scores = Float64[]
        for (_, local_ids) in byq
            append!(scores, [nodes[id].score for id in local_ids])
            if length(local_ids) >= 2
                sorted = sort(local_ids; by=id -> nodes[id].iω)
                push!(seps, Float64(nodes[sorted[end]].iω - nodes[sorted[1]].iω))
            end
        end
        medsep = isempty(seps) ? 0.0 : median(seps)
        medscore = isempty(scores) ? 0.0 : median(scores)
        conf = clamp(0.35 * min(qfrac / min_cycle_q_fraction, 1.0) +
                     0.35 * min(two_frac / min_two_peak_fraction, 1.0) +
                     0.30 * min(medscore / min_median_score, 1.0), 0.0, 1.0)
        accepted = qfrac >= min_cycle_q_fraction && two_frac >= min_two_peak_fraction && medscore >= min_median_score && conf >= min_confidence
        reason = accepted ? "accepted" : "rejected: insufficient q-span, two-peak support, or ridge score"
        push!(cands, RidgeCycleCandidate(comp, β1c, (qlo, qhi), qfrac, two_frac, medsep, medscore, conf, accepted, reason))
    end
    return cands
end


#----------------------------#
# Split-Merge Pattern Rescue #
#----------------------------#
function _nearest_singleton_slice(peak_counts::Vector{Int}, start::Int, stop::Int, step::Int)
    i = start
    while step > 0 ? i <= stop : i >= stop
        peak_counts[i] == 1 && return i
        i += step
    end
    return 0
end

function _edge_median(xs::Vector{Float64}; frac::Float64=0.12)
    isempty(xs) && return 0.0
    n = length(xs)
    k = clamp(ceil(Int, frac*n), 1, n)
    return median(vcat(xs[1:k], xs[max(1, n-k+1):n]))
end

function _center_median(xs::Vector{Float64}; frac::Float64=0.35)
    isempty(xs) && return 0.0
    n = length(xs)
    k = clamp(floor(Int, frac*n), 0, max(0, (n-1) ÷ 2))
    lo = 1 + k
    hi = n - k
    lo > hi && return median(xs)
    return median(@view xs[lo:hi])
end

function _splitmerge_pattern_run(peaks_local::Vector{Vector{Int}}, scores_local::Vector{Vector{Float64}}, nω::Int;
                                 min_cycle_q_fraction::Float64=0.08,
                                 min_two_peak_fraction::Float64=0.15,
                                 min_median_score::Float64=0.10,
                                 min_confidence::Float64=0.55,
                                 min_branch_separation_bins::Int=4,
                                 max_branch_separation_frac::Float64=0.28)
    # Detect resolved split/rejoin intervals from ridge multiplicity.
    # The lens profile rejects ordinary weak independent branches.

    nq = length(peaks_local)
    peak_counts = length.(peaks_local)
    max_sep = max(min_branch_separation_bins, floor(Int, max_branch_separation_frac * nω))

    best = nothing
    i = 1
    while i <= nq
        if peak_counts[i] < 2
            i += 1
            continue
        end

        qlo = i
        while i <= nq && peak_counts[i] >= 2
            i += 1
        end
        qhi = i - 1
        qspan = qhi - qlo + 1
        qfrac = qspan / max(nq, 1)
        two_frac = count(j -> peak_counts[j] >= 2, qlo:qhi) / qspan

        seps = Float64[]
        scores = Float64[]
        for q in qlo:qhi
            ps = peaks_local[q]
            sc = scores_local[q]
            length(ps) >= 2 || continue
            push!(seps, Float64(maximum(ps) - minimum(ps)))
            append!(scores, sc)
        end
        isempty(seps) && continue

        medsep = median(seps)
        medscore = isempty(scores) ? 0.0 : median(scores)
        min_branch_separation_bins <= medsep <= max_sep || continue

        edge_sep = _edge_median(seps)
        center_sep = _center_median(seps)
        lens_gain = center_sep - edge_sep
        lens_ratio = center_sep / max(edge_sep, 1.0)

        # Broadened split/merge runs may not have near-zero endpoint separation.
        # Use a softer lens test while preserving weak-optical rejection.
        min_lens_gain = max(1.0, 0.35 * min_branch_separation_bins)
        lens_ok = (lens_gain >= min_lens_gain) && (lens_ratio >= 1.05)

        left = _nearest_singleton_slice(peak_counts, qlo - 1, max(1, qlo - max(6, qspan)), -1)
        right = _nearest_singleton_slice(peak_counts, qhi + 1, min(nq, qhi + max(6, qspan)), 1)

        midpoint_votes = 0
        midpoint_possible = 0
        if left != 0 && !isempty(peaks_local[left])
            midpoint_possible += 1
            lo_left, hi_left = minimum(peaks_local[qlo]), maximum(peaks_local[qlo])
            left_single = peaks_local[left][1]
            if _valid_splitmerge_geometry(left_single, lo_left, hi_left;
                                          min_sep=min_branch_separation_bins, max_sep=max_sep,
                                          midpoint_tol=Float64(max(min_branch_separation_bins, 
                                                                   2*min_branch_separation_bins)))
                midpoint_votes += 1
            end
        end
        if right != 0 && !isempty(peaks_local[right])
            midpoint_possible += 1
            lo_right, hi_right = minimum(peaks_local[qhi]), maximum(peaks_local[qhi])
            right_single = peaks_local[right][1]
            if _valid_splitmerge_geometry(right_single, lo_right, hi_right;
                                          min_sep=min_branch_separation_bins, max_sep=max_sep,
                                          midpoint_tol=Float64(max(min_branch_separation_bins, 
                                                                   2*min_branch_separation_bins)))
                midpoint_votes += 1
            end
        end
        midpoint_ok = midpoint_possible == 0 ? lens_ok : (midpoint_votes >= 1 || lens_ok)
        (lens_ok && midpoint_ok) || continue

        conf = clamp(0.30 * min(qfrac / min_cycle_q_fraction, 1.0) +
                     0.25 * min(two_frac / min_two_peak_fraction, 1.0) +
                     0.20 * min(medscore / min_median_score, 1.0) +
                     0.25 * min(lens_ratio / 1.50, 1.0), 0.0, 1.0)
        accepted = qfrac >= min_cycle_q_fraction && two_frac >= min_two_peak_fraction && medscore >= min_median_score && conf >= min_confidence
        accepted || continue

        score = conf + 0.05*qfrac + 0.05*min(lens_ratio/2, 1.0)
        if best === nothing || score > best[:score]
            best = Dict{Symbol,Any}(:qlo=>qlo, :qhi=>qhi, :left=>left, :right=>right,
                                    :qfrac=>qfrac, :two_frac=>two_frac, :medsep=>medsep,
                                    :edge_sep=>edge_sep, :center_sep=>center_sep,
                                    :lens_ratio=>lens_ratio, :medscore=>medscore,
                                    :confidence=>conf, :score=>score)
        end
    end
    return best
end

function _node_at_slice_energy(nodes::Vector{RidgeNode}, node_ids::Vector{Vector{Int}}, iq::Int, target_iω::Int)
    ids = node_ids[iq]
    isempty(ids) && return 0
    return _nearest_node_id(nodes, ids, target_iω)
end

function add_splitmerge_pattern_edges!(edges::Vector{Tuple{Int,Int}}, nodes::Vector{RidgeNode},
                                       node_ids::Vector{Vector{Int}}, peaks_local::Vector{Vector{Int}},
                                       scores_local::Vector{Vector{Float64}}, nω::Int;
                                       min_cycle_q_fraction::Float64=0.08,
                                       min_two_peak_fraction::Float64=0.15,
                                       min_median_score::Float64=0.10,
                                       min_confidence::Float64=0.55,
                                       min_branch_separation_bins::Int=4,
                                       max_branch_separation_frac::Float64=0.28)
    run = _splitmerge_pattern_run(peaks_local, scores_local, nω;
                                  min_cycle_q_fraction=min_cycle_q_fraction,
                                  min_two_peak_fraction=min_two_peak_fraction,
                                  min_median_score=min_median_score,
                                  min_confidence=min_confidence,
                                  min_branch_separation_bins=min_branch_separation_bins,
                                  max_branch_separation_frac=max_branch_separation_frac)
    run === nothing && return 0, nothing

    qlo = run[:qlo]; qhi = run[:qhi]
    left = run[:left]; right = run[:right]

    lo_left, hi_left = minimum(peaks_local[qlo]), maximum(peaks_local[qlo])
    lo_right, hi_right = minimum(peaks_local[qhi]), maximum(peaks_local[qhi])

    b1 = _node_at_slice_energy(nodes, node_ids, qlo, lo_left)
    b2 = _node_at_slice_energy(nodes, node_ids, qlo, hi_left)
    c1 = _node_at_slice_energy(nodes, node_ids, qhi, lo_right)
    c2 = _node_at_slice_energy(nodes, node_ids, qhi, hi_right)

    added = 0
    candidate_edges = Tuple{Int,Int}[]

    # Virtual closure edges encode inferred branch points.
    # This helps when broad one-peak slices defeat the local linker.
    push!(candidate_edges, (b1, b2))
    push!(candidate_edges, (c1, c2))

    if left != 0 && !isempty(peaks_local[left])
        a = _node_at_slice_energy(nodes, node_ids, left, peaks_local[left][1])
        push!(candidate_edges, (a, b1)); push!(candidate_edges, (a, b2))
    end
    if right != 0 && !isempty(peaks_local[right])
        d = _node_at_slice_energy(nodes, node_ids, right, peaks_local[right][1])
        push!(candidate_edges, (c1, d)); push!(candidate_edges, (c2, d))
    end

    for e in candidate_edges
        if e[1] != 0 && e[2] != 0 && e[1] != e[2] && !(e in edges)
            push!(edges, e)
            added += 1
        end
    end
    return added, run
end

function build_ridge_topology(Iqω::AbstractMatrix{<:Real}, qs, ωs;
                              smooth_radius_ω::Int=2,
                              include_energy_edges::Bool=false,
                              min_prominence::Union{Nothing,Real}=nothing,
                              max_peaks_per_slice::Int=2,
                              max_jump_bins::Union{Nothing,Int}=nothing,
                              tube_radius_ω::Int=1,
                              periodic_q::Bool=false,
                              pixel_adjacency::Symbol=:full,
                              relative_score_floor::Float64=0.10,
                              min_branch_separation_bins::Union{Nothing,Int}=nothing,
                              max_branch_separation_frac::Float64=0.32,
                              min_cycle_q_fraction::Float64=0.08,
                              min_two_peak_fraction::Float64=0.15,
                              min_cycle_confidence::Float64=0.45)
    nq, nω = size(Iqω)
    peaks_local, scores_local = ridge_peak_slices(Iqω, ωs; smooth_radius_ω=smooth_radius_ω,
                                                  include_energy_edges=include_energy_edges,
                                                  min_prominence=min_prominence,
                                                  max_peaks_per_slice=max_peaks_per_slice,
                                                  relative_score_floor=relative_score_floor)

    nodes = RidgeNode[]
    node_ids = [Int[] for _ in 1:nq]
    for iq in 1:nq
        for (k, iω) in enumerate(peaks_local[iq])
            push!(nodes, RidgeNode(iq, iω, Float64(qs[iq]), Float64(ωs[iω]), scores_local[iq][k]))
            push!(node_ids[iq], length(nodes))
        end
    end

    jump = isnothing(max_jump_bins) ? max(4, 4 * smooth_radius_ω) : Int(max_jump_bins)
    min_sep = isnothing(min_branch_separation_bins) ? max(4, 2 * smooth_radius_ω) : Int(min_branch_separation_bins)
    max_sep = max(min_sep, floor(Int, max_branch_separation_frac * nω))
    midpoint_tol = Float64(max(jump, min_sep))
    edges = Tuple{Int,Int}[]
    nlinks = periodic_q ? nq : nq - 1

    for iq in 1:nlinks
        jq = iq == nq ? 1 : iq + 1
        a_ids = node_ids[iq]
        b_ids = node_ids[jq]
        if isempty(a_ids) || isempty(b_ids)
            continue
        end
        a_sorted = sort(a_ids; by=id -> nodes[id].iω)
        b_sorted = sort(b_ids; by=id -> nodes[id].iω)
        na = length(a_sorted)
        nb = length(b_sorted)

        if na == nb
            _add_ordered_matching_edges!(edges, nodes, a_sorted, b_sorted, jump)
        elseif na == 1 && nb == 2
            s = a_sorted[1]
            lo, hi = b_sorted[1], b_sorted[2]
            if _valid_splitmerge_geometry(nodes[s].iω, nodes[lo].iω, nodes[hi].iω; min_sep=min_sep, max_sep=max_sep, midpoint_tol=midpoint_tol)
                push!(edges, (s, lo)); push!(edges, (s, hi))
            else
                b = _nearest_node_id(nodes, b_sorted, nodes[s].iω)
                abs(nodes[s].iω - nodes[b].iω) <= max(jump, min_sep) && push!(edges, (s, b))
            end
        elseif na == 2 && nb == 1
            s = b_sorted[1]
            lo, hi = a_sorted[1], a_sorted[2]
            if _valid_splitmerge_geometry(nodes[s].iω, nodes[lo].iω, nodes[hi].iω; min_sep=min_sep, max_sep=max_sep, midpoint_tol=midpoint_tol)
                push!(edges, (lo, s)); push!(edges, (hi, s))
            else
                a = _nearest_node_id(nodes, a_sorted, nodes[s].iω)
                abs(nodes[a].iω - nodes[s].iω) <= max(jump, min_sep) && push!(edges, (a, s))
            end
        else
            _add_ordered_matching_edges!(edges, nodes, a_sorted, b_sorted, jump)
        end
    end

    pattern_edges_added, splitmerge_pattern = add_splitmerge_pattern_edges!(edges, nodes, node_ids, peaks_local, scores_local, nω;
                                                                            min_cycle_q_fraction=min_cycle_q_fraction,
                                                                            min_two_peak_fraction=min_two_peak_fraction,
                                                                            min_median_score=0.10,
                                                                            min_confidence=min_cycle_confidence,
                                                                            min_branch_separation_bins=min_sep,
                                                                            max_branch_separation_frac=max_branch_separation_frac)

    beta0, beta1_raw = graph_betti_numbers(length(nodes), edges)
    comp_labels, ngraph_components = graph_component_labels(length(nodes), edges)
    candidates = ridge_cycle_candidates(nodes, edges, length.(peaks_local), ngraph_components, comp_labels;
                                        min_cycle_q_fraction=min_cycle_q_fraction,
                                        min_two_peak_fraction=min_two_peak_fraction,
                                        min_confidence=min_cycle_confidence)
    # Count accepted split/merge events, not raw graph cycle rank.
    # Rescue edges can create multiple cycles around one physical lens.
    accepted_cycles = isempty(candidates) ? 0 : count(c -> c.accepted, candidates)
    beta1 = accepted_cycles
    status = beta1 > 0 ? :splitmerge_detected : :no_splitmerge_cycle

    ridge_mask = rasterize_ridge_graph(nodes, edges, nq, nω;
                                       tube_radius_ω=tube_radius_ω,
                                       periodic_q=periodic_q)
    ridge_labels, ncomponents = label_support_components(ridge_mask;
                                                         pixel_adjacency=pixel_adjacency,
                                                         periodic_q=periodic_q)
    diagnostics = Dict{Symbol,Any}(
        :raw_beta1 => beta1_raw,
        :graph_components => ngraph_components,
        :jump_bins => jump,
        :min_branch_separation_bins => min_sep,
        :max_branch_separation_bins => max_sep,
        :relative_score_floor => relative_score_floor,
        :min_cycle_q_fraction => min_cycle_q_fraction,
        :min_two_peak_fraction => min_two_peak_fraction,
        :min_cycle_confidence => min_cycle_confidence,
        :splitmerge_pattern_edges_added => pattern_edges_added,
        :splitmerge_pattern => splitmerge_pattern,
    )
    return RidgeTopologyResult(nodes, edges, beta0, beta1,
                               length.(peaks_local), ridge_mask,
                               ridge_labels, ncomponents,
                               candidates, accepted_cycles,
                               status, diagnostics)
end

#----------------------#
# Parent Support Masks #
#----------------------#
@inline function logistic_confidence(x::Float64)
    x >= 40 && return 1.0
    x <= -40 && return 0.0
    return 1.0 / (1.0 + exp(-x))
end

function build_parent_mask!(mask::BitMatrix, S::Matrix{Float64}, domains,
                            parent::Int, global_bg::Float64, global_σ::Float64;
                            kσ_support::Float64=3.5,
                            local_fraction_floor::Float64=0.80,
                            confidence_min::Float64=0.75)
    # Soft support uses confidence above adaptive background/local thresholds.
    # Binarization is kept only for connected-component labeling.
    fill!(mask, false)
    nq, nω = size(S)
    attempted = 0
    active = 0
    σeff = max(global_σ, eps(Float64))
    @inbounds for iq in 1:nq
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue
        attempted += length(r)
        local_peak = maximum(@view S[iq, r])
        τ_abs = global_bg + kσ_support * global_σ
        τ_rel = local_fraction_floor * local_peak
        τ = max(τ_abs, τ_rel)
        for iω in r
            conf = logistic_confidence((S[iq, iω] - τ) / σeff)
            if conf >= confidence_min
                mask[iq, iω] = true
                active += 1
            end
        end
    end
    return attempted, active
end

function row_activity(mask::BitMatrix, domains, parent::Int)
    nq, _ = size(mask)
    active = falses(nq)
    @inbounds for iq in 1:nq
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue
        for iω in r
            if mask[iq, iω]
                active[iq] = true
                break
            end
        end
    end
    return active
end

function close_short_gaps(active::BitVector, max_gap::Int; periodic_q::Bool=true)
    n = length(active)
    closed = copy(active)
    n == 0 && return closed
    nactive = count(active)
    (nactive == 0 || nactive == n || max_gap <= 0) && return closed

    if periodic_q
        start = findfirst(active)
        i = mod1(start + 1, n)
        run = Int[]
        while i != start
            if active[i]
                if !isempty(run) && length(run) <= max_gap
                    for j in run
                        closed[j] = true
                    end
                end
                empty!(run)
            else
                push!(run, i)
            end
            i = mod1(i + 1, n)
        end
        # Final run, if present, is the wrap-around bounded gap.
        if !isempty(run) && length(run) <= max_gap
            for j in run
                closed[j] = true
            end
        end
    else
        i = 1
        while i <= n
            if active[i]
                i += 1
                continue
            end
            lo = i
            while i <= n && !active[i]
                i += 1
            end
            hi = i - 1
            bounded = lo > 1 && hi < n && active[lo - 1] && active[hi + 1]
            if bounded && (hi - lo + 1) <= max_gap
                for j in lo:hi
                    closed[j] = true
                end
            end
        end
    end
    return closed
end

function repair_parent_mask_continuity!(mask::BitMatrix, S::Matrix{Float64}, domains,
                                        centers::Matrix{Int}, parent::Int;
                                        radius_q::Int=1,
                                        radius_ω::Int=2,
                                        periodic_q::Bool=true)
    # Resolution-scale continuity repair, not blind dilation.
    # Only short q-gaps near a tracked centerline are bridged.
    nq, nω = size(mask)
    parent <= size(centers, 1) || return 0, 0, 0

    before = count(mask)
    active = row_activity(mask, domains, parent)
    max_gap_q_bins = max(2, 2 * radius_q + 1)
    bridge_halfwidth_ω = max(1, cld(radius_ω, 2))
    closed = close_short_gaps(active, max_gap_q_bins; periodic_q=periodic_q)
    repaired_rows = 0

    @inbounds for iq in 1:nq
        closed[iq] || continue
        active[iq] && continue
        repaired_rows += 1
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue
        c = clamp(centers[parent, iq], first(r), last(r))
        lo = max(first(r), c - bridge_halfwidth_ω)
        hi = min(last(r), c + bridge_halfwidth_ω)
        for iω in lo:hi
            mask[iq, iω] = true
        end
    end

    after = count(mask)
    return max(after - before, 0), repaired_rows, max_gap_q_bins
end

function connect_parent_row_fragments!(mask::BitMatrix, S::Matrix{Float64}, domains,
                                       parent::Int,
                                       global_bg::Float64, global_σ::Float64;
                                       radius_ω::Int=2,
                                       max_gap_factor::Float64=1.25,
                                       kσ_gap::Float64=0.25)
    # Connect short vertical holes within one parent domain on a q-row.
    # This handles resolvable overlaps with no inactive q-row.
    nq, nω = size(mask)
    parent <= length(domains[1]) || return 0, 0

    max_gap_ω = max(1, ceil(Int, max_gap_factor * radius_ω))
    τ_gap = global_bg + kσ_gap * max(global_σ, eps(Float64))

    before = count(mask)
    touched_rows = 0

    @inbounds for iq in 1:nq
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue

        last_hi = 0
        have_run = false
        row_touched = false

        iω = first(r)
        while iω <= last(r)
            if !mask[iq, iω]
                iω += 1
                continue
            end

            lo = iω
            while iω <= last(r) && mask[iq, iω]
                iω += 1
            end
            hi = iω - 1

            if have_run
                gap_lo = last_hi + 1
                gap_hi = lo - 1
                gap = gap_hi - gap_lo + 1

                if 1 <= gap <= max_gap_ω
                    gap_peak = maximum(@view S[iq, gap_lo:gap_hi])
                    left_edge = S[iq, last_hi]
                    right_edge = S[iq, lo]

                    # Bridge if the gap or both adjacent supports carry evidence.
                    if gap_peak >= τ_gap || min(left_edge, right_edge) >= τ_gap
                        for j in gap_lo:gap_hi
                            mask[iq, j] = true
                        end
                        row_touched = true
                    end
                end
            end

            last_hi = hi
            have_run = true
        end

        row_touched && (touched_rows += 1)
    end

    after = count(mask)
    return max(after - before, 0), touched_rows
end

function enforce_parent_centerline_support!(mask::BitMatrix, S::Matrix{Float64}, domains,
                                            centers::Matrix{Int}, imputed::BitMatrix,
                                            parent::Int,
                                            global_bg::Float64, global_σ::Float64;
                                            radius_q::Int=1,
                                            radius_ω::Int=2,
                                            periodic_q::Bool=true,
                                            kσ_center::Float64=2.0,
                                            parent_fraction_floor::Float64=0.025,
                                            imputed_kσ_center::Float64=1.25,
                                            imputed_fraction_floor::Float64=0.030,
                                            imputed_search_radius_factor::Int=2)
    nq, nω = size(mask)
    parent <= size(centers, 1) || return 0, 0

    parent_peak = 0.0
    @inbounds for iq in 1:nq
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue
        parent_peak = max(parent_peak, maximum(@view S[iq, r]))
    end

    σeff = max(global_σ, eps(Float64))

    τ_detected = max(global_bg + kσ_center * σeff, parent_fraction_floor * parent_peak)

    # Imputed rows often occur where branches collapse into one maximum.
    # Require evidence, but use a weaker threshold than detected rows.
    τ_imputed = max(global_bg + imputed_kσ_center * σeff, imputed_fraction_floor * parent_peak)

    search_radius_detected = max(1, radius_ω)
    search_radius_imputed = max(search_radius_detected, imputed_search_radius_factor * radius_ω)

    tube_halfwidth_ω = max(2, cld(radius_ω, 2))
    max_gap_q_bins = max(2, 2 * radius_q + 1)

    supported = falses(nq)
    center_idx = zeros(Int, nq)

    # First pass: detected rows follow tracked centers.
    # Imputed rows do not snap to the global parent-domain maximum.
    @inbounds for iq in 1:nq
        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue

        c0 = clamp(centers[parent, iq], first(r), last(r))
        sr = imputed[parent, iq] ? search_radius_imputed : search_radius_detected

        slo = max(first(r), c0 - sr)
        shi = min(last(r), c0 + sr)

        c = c0
        best = S[iq, c0]
        for iω in slo:shi
            if S[iq, iω] > best
                best = S[iq, iω]
                c = iω
            end
        end

        threshold = imputed[parent, iq] ? τ_imputed : τ_detected
        if best >= threshold
            supported[iq] = true
            center_idx[iq] = c
        end
    end

    # Second pass: interpolate imputed runs between non-imputed centers.
    # Evidence tests prevent bridging truly empty q-gaps.
    @inbounds for iq in 1:nq
        imputed[parent, iq] || continue
        supported[iq] && continue

        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue

        left = 0
        right = 0

        for d in 1:nq
            il = periodic_q ? periodic_index(iq - d, nq) : iq - d
            1 <= il <= nq || break
            if !imputed[parent, il] && centers[parent, il] != 0
                left = il
                break
            end
        end

        for d in 1:nq
            ir = periodic_q ? periodic_index(iq + d, nq) : iq + d
            1 <= ir <= nq || break
            if !imputed[parent, ir] && centers[parent, ir] != 0
                right = ir
                break
            end
        end

        left == 0 && right == 0 && continue

        cpath = if left != 0 && right != 0
            dl = periodic_q ? min(abs(iq - left), nq - abs(iq - left)) : abs(iq - left)
            dr = periodic_q ? min(abs(right - iq), nq - abs(right - iq)) : abs(right - iq)
            denom = max(dl + dr, 1)
            round(Int, (dr * centers[parent, left] + dl * centers[parent, right]) / denom)
        elseif left != 0
            centers[parent, left]
        else
            centers[parent, right]
        end

        cpath = clamp(cpath, first(r), last(r))

        sr = search_radius_imputed
        slo = max(first(r), cpath - sr)
        shi = min(last(r), cpath + sr)

        c = cpath
        best = S[iq, cpath]
        for iω in slo:shi
            if S[iq, iω] > best
                best = S[iq, iω]
                c = iω
            end
        end

        if best >= τ_imputed
            supported[iq] = true
            center_idx[iq] = c
        end
    end

    closed = close_short_gaps(supported, max_gap_q_bins; periodic_q=periodic_q)

    @inbounds for iq in 1:nq
        if closed[iq] && center_idx[iq] == 0
            ranges = domains[iq]
            parent <= length(ranges) || continue
            r = ranges[parent]
            isempty(r) && continue

            left = 0
            right = 0

            for d in 1:max_gap_q_bins+1
                il = periodic_q ? periodic_index(iq - d, nq) : iq - d
                if 1 <= il <= nq && supported[il] && center_idx[il] != 0
                    left = il
                    break
                end
            end

            for d in 1:max_gap_q_bins+1
                ir = periodic_q ? periodic_index(iq + d, nq) : iq + d
                if 1 <= ir <= nq && supported[ir] && center_idx[ir] != 0
                    right = ir
                    break
                end
            end

            if left != 0 && right != 0
                center_idx[iq] = clamp(round(Int, 0.5 * (center_idx[left] + center_idx[right])),
                                       first(r), last(r))
            elseif left != 0
                center_idx[iq] = clamp(center_idx[left], first(r), last(r))
            elseif right != 0
                center_idx[iq] = clamp(center_idx[right], first(r), last(r))
            else
                center_idx[iq] = clamp(centers[parent, iq], first(r), last(r))
            end
        end
    end

    before = count(mask)
    anchored_rows = 0

    @inbounds for iq in 1:nq
        closed[iq] || continue
        c = center_idx[iq]
        c == 0 && continue

        ranges = domains[iq]
        parent <= length(ranges) || continue
        r = ranges[parent]
        isempty(r) && continue

        anchored_rows += 1

        lo = max(first(r), c - tube_halfwidth_ω)
        hi = min(last(r), c + tube_halfwidth_ω)

        for iω in lo:hi
            mask[iq, iω] = true
        end
    end

    after = count(mask)
    return max(after - before, 0), anchored_rows
end
#-----------------#
# Pixel Adjacency #
#-----------------#
const EDGE_PIXEL_NEIGHBORS = ((-1,0), (1,0), (0,-1), (0,1))
const FULL_PIXEL_NEIGHBORS = ((-1,-1), (-1,0), (-1,1),
                              (0,-1),           (0,1),
                              (1,-1),  (1,0),   (1,1))

@inline function pixel_neighbors(pixel_adjacency::Symbol)
    pixel_adjacency === :edge && return EDGE_PIXEL_NEIGHBORS
    pixel_adjacency === :full && return FULL_PIXEL_NEIGHBORS
    throw(ArgumentError("pixel_adjacency must be :edge or :full"))
end

#--------------------#
# Component Labeling #
#--------------------#
function label_support_components(mask::BitMatrix;
                                  pixel_adjacency::Symbol=:full,
                                  periodic_q::Bool=true)
    nq, nω = size(mask)
    labels = zeros(Int32, nq, nω)
    ncomp = 0

    neighbors = pixel_neighbors(pixel_adjacency)

    qqueue = Vector{Int}(undef, nq*nω)
    wqueue = Vector{Int}(undef, nq*nω)

    @inbounds for iq0 in 1:nq, iw0 in 1:nω
        if mask[iq0, iw0] && labels[iq0, iw0] == 0
            ncomp += 1
            head = 1
            tail = 1
            qqueue[tail] = iq0
            wqueue[tail] = iw0
            labels[iq0, iw0] = ncomp

            while head <= tail
                q0 = qqueue[head]
                w0 = wqueue[head]
                head += 1

                for (dq, dw) in neighbors
                    qraw = q0 + dq
                    q1 = periodic_q ? periodic_index(qraw, nq) : qraw
                    w1 = w0 + dw

                    if 1 <= q1 <= nq && 1 <= w1 <= nω
                        if mask[q1, w1] && labels[q1, w1] == 0
                            labels[q1, w1] = ncomp
                            tail += 1
                            qqueue[tail] = q1
                            wqueue[tail] = w1
                        end
                    end
                end
            end
        end
    end

    return labels, ncomp
end

function label_components_and_summarize!(components::Vector{ComponentResult},
                                         mask::BitMatrix, Iqω::Matrix{Float64},
                                         qs::Vector{Float64}, ωs::Vector{Float64},
                                         parent::Int;
                                         pixel_adjacency::Symbol=:full,
                                         periodic_q::Bool=true,
                                         min_resolution_pixels::Int=12)
    nq, nω = size(mask)
    visited = falses(nq, nω)
    ncomp = 0
    kept_for_parent = 0
    Δq = mean(diff(qs))
    Δω = mean(diff(ωs))
    lin = LinearIndices(mask)

    neighbors = pixel_neighbors(pixel_adjacency)

    qqueue = Vector{Int}(undef, nq*nω)
    wqueue = Vector{Int}(undef, nq*nω)
    pixbuf = Vector{Int}(undef, nq*nω)
    scratch_Iq = zeros(Float64, nq)
    scratch_Iω = zeros(Float64, nω)

    @inbounds for iq0 in 1:nq, iw0 in 1:nω
        if mask[iq0, iw0] && !visited[iq0, iw0]
            ncomp += 1
            fill!(scratch_Iq, 0.0)
            fill!(scratch_Iω, 0.0)

            head = 1
            tail = 1
            qqueue[tail] = iq0
            wqueue[tail] = iw0
            visited[iq0, iw0] = true

            qmin = iq0
            qmax = iq0
            wmin = iw0
            wmax = iw0
            npix = 0
            weight = 0.0

            while head <= tail
                q0 = qqueue[head]
                w0 = wqueue[head]
                head += 1

                npix += 1
                pixbuf[npix] = lin[q0, w0]
                val = Iqω[q0, w0]
                weight += val
                scratch_Iq[q0] += val * Δω
                scratch_Iω[w0] += val * Δq
                qmin = min(qmin, q0)
                qmax = max(qmax, q0)
                wmin = min(wmin, w0)
                wmax = max(wmax, w0)

                for (dq, dw) in neighbors
                    qraw = q0 + dq
                    q1 = periodic_q ? periodic_index(qraw, nq) : qraw
                    w1 = w0 + dw
                    if 1 <= q1 <= nq && 1 <= w1 <= nω
                        if mask[q1, w1] && !visited[q1, w1]
                            visited[q1, w1] = true
                            tail += 1
                            qqueue[tail] = q1
                            wqueue[tail] = w1
                        end
                    end
                end
            end

            if npix >= min_resolution_pixels
                kept_for_parent += 1
                push!(components, ComponentResult(parent, kept_for_parent,
                                                  (qs[qmin], qs[qmax]), (ωs[wmin], ωs[wmax]),
                                                  qmin:qmax, wmin:wmax, npix, weight,
                                                  copy(scratch_Iq), copy(scratch_Iω),
                                                  copy(@view pixbuf[1:npix])))
            end
        end
    end

    return ncomp, kept_for_parent
end

#------------------------------#
# Adaptive Segmentation Driver #
#------------------------------#
"""
    adaptive_autosplitbands2d(Iqω, qs, ωs; kwargs...)

Segment a 2D intensity matrix `Iqω` into physically tracked parent bands and
connected support components. The first dimension is momentum `q`, the second
is energy `ω`.

Returns a named tuple containing parent tracking data, final components,
component counts, smoothed support images, resolution estimates, bridge
statistics, and optional ridge-topology split/merge diagnostics.

Key options include `periodic_q`, `resolution_q`, `resolution_ω`,
`parent_count_method`, `max_bands`, `pixel_adjacency`, and `ridge_*` controls.
Use `pixel_adjacency=:edge` for edge-sharing support and `:full` for the
full 3×3 neighborhood excluding the center pixel.
"""
function adaptive_autosplitbands2d(Iqω_in::AbstractMatrix{<:Real}, qs, ωs;
                                   periodic_q::Bool=true,
                                   include_energy_edges::Bool=true,
                                   resolution_q=nothing,
                                   resolution_ω=nothing,
                                   support_smooth_radius_ω=nothing,
                                   peak_smooth_radius_ω=nothing,
                                   max_support_smooth_radius_ω::Int=9,
                                   max_peak_smooth_radius_ω::Int=6,
                                   parent_count_method::Symbol=:quantile,
                                   max_bands::Union{Nothing,Int}=nothing,
                                   pixel_adjacency::Symbol=:full,
                                   ridge_topology::Bool=true,
                                   ridge_smooth_radius_ω=nothing,
                                   ridge_max_peaks_per_slice::Int=2,
                                   ridge_max_jump_bins=nothing,
                                   ridge_tube_radius_ω::Int=2,
                                   ridge_periodic_q::Bool=false,
                                   ridge_relative_score_floor::Float64=0.10,
                                   ridge_min_branch_separation_bins=nothing,
                                   ridge_max_branch_separation_frac::Float64=0.28,
                                   ridge_min_cycle_q_fraction::Float64=0.08,
                                   ridge_min_two_peak_fraction::Float64=0.15,
                                   ridge_min_cycle_confidence::Float64=0.55,
                                   verbose::Bool=true)
    Iqω = Matrix{Float64}(Iqω_in)
    nq, nω = size(Iqω)
    radius_q, radius_ω_est, Γq, Γω = adaptive_resolution_radii(qs, ωs, Iqω; resolution_q=resolution_q, resolution_ω=resolution_ω)

    # Half-height estimates can inflate in overlap regions.
    # Do not use them directly for masks or peak tracking.
    support_radius_ω = isnothing(support_smooth_radius_ω) ? min(radius_ω_est, max_support_smooth_radius_ω) : max(1, Int(support_smooth_radius_ω))

    peak_radius_ω = isnothing(peak_smooth_radius_ω) ? min(support_radius_ω, max_peak_smooth_radius_ω) : max(1, Int(peak_smooth_radius_ω))

    # Support image: moderately smoothed, used for binary masks/components.
    S = Matrix{Float64}(undef, nq, nω)
    separable_smooth2d!(S, Iqω; radius_q=radius_q, radius_ω=support_radius_ω, periodic_q=periodic_q)

    bg, σ = robust_background_sigma(S)
    minpix = max(4, (2*radius_q + 1) * (2*support_radius_ω + 1))

    # Peak-tracking image is q-smoothed only.
    # Energy smoothing occurs inside qfirst_parent_domains.
    S_peak = Matrix{Float64}(undef, nq, nω)
    separable_smooth2d!(S_peak, Iqω; radius_q=radius_q, radius_ω=0, periodic_q=periodic_q)

    qfirst = qfirst_parent_domains(S_peak, qs, ωs;
                                   smooth_radius_ω=peak_radius_ω,
                                   include_energy_edges=include_energy_edges,
                                   periodic_q=periodic_q,
                                   min_prominence=nothing,
                                   parent_count_method=parent_count_method,
                                   max_bands=max_bands,
                                   verbose=verbose)



    mask = falses(nq, nω)
    components = ComponentResult[]
    attempted = zeros(Int, qfirst.nparents)
    active = zeros(Int, qfirst.nparents)
    raw_counts = zeros(Int, qfirst.nparents)
    centerline_pixels = zeros(Int, qfirst.nparents)
    centerline_rows = zeros(Int, qfirst.nparents)
    bridge_pixels = zeros(Int, qfirst.nparents)
    bridge_rows = zeros(Int, qfirst.nparents)
    bridge_max_gap = zeros(Int, qfirst.nparents)
    for b in 1:qfirst.nparents
        attempted[b], active[b] = build_parent_mask!(mask, S, qfirst.domains, b, bg, σ; kσ_support=3.5, local_fraction_floor=0.45, confidence_min=0.17)
        centerline_pixels[b], centerline_rows[b] = enforce_parent_centerline_support!(mask, S, qfirst.domains, qfirst.centers, qfirst.imputed, b, bg, σ;
                                                                                      radius_q=radius_q,
                                                                                      radius_ω=support_radius_ω,
                                                                                      periodic_q=periodic_q,
                                                                                      kσ_center=2.0,
                                                                                      parent_fraction_floor=0.025)

        fragment_pixels, fragment_rows = connect_parent_row_fragments!(mask, S, qfirst.domains, b, bg, σ; 
                                                                       radius_ω=support_radius_ω, 
                                                                       max_gap_factor=1.25,
                                                                       kσ_gap=0.25)

        centerline_pixels[b] += fragment_pixels
        centerline_rows[b] = max(centerline_rows[b], fragment_rows)

        bridge_pixels[b], bridge_rows[b], bridge_max_gap[b] = repair_parent_mask_continuity!(mask, S, qfirst.domains, qfirst.centers, b;
                                                                                             radius_q=radius_q,
                                                                                             radius_ω=support_radius_ω,
                                                                                             periodic_q=periodic_q)

        raw_counts[b], kept = label_components_and_summarize!(components, mask, Iqω, qs, ωs, b;
                                                              pixel_adjacency=pixel_adjacency,
                                                              periodic_q=periodic_q,
                                                              min_resolution_pixels=minpix)

        verbose && @info "parent band $b: attempted=$(attempted[b]), active=$(active[b]), active fraction=$(round(active[b]/max(attempted[b],1); digits=4)), centerline_pixels=$(centerline_pixels[b]), centerline_rows=$(centerline_rows[b]), bridge_pixels=$(bridge_pixels[b]), bridge_rows=$(bridge_rows[b]), max_gap_q=$(bridge_max_gap[b]), raw components=$(raw_counts[b]), kept components=$kept"
    end

    component_counts = [count(c -> c.parent_band == b, components) for b in 1:qfirst.nparents]
    r_smooth = isnothing(ridge_smooth_radius_ω) ? max(1, min(6, cld(radius_ω_est, 6))) : Int(ridge_smooth_radius_ω)
    effective_ridge_max_peaks = max(ridge_max_peaks_per_slice, qfirst.nparents)
    ridge_topology_result = ridge_topology ? build_ridge_topology(Iqω, qs, ωs; 
                                                                  smooth_radius_ω=r_smooth, 
                                                                  include_energy_edges=false,
                                                                  min_prominence=nothing,
                                                                  max_peaks_per_slice=effective_ridge_max_peaks,
                                                                  max_jump_bins=ridge_max_jump_bins,
                                                                  tube_radius_ω=ridge_tube_radius_ω,
                                                                  periodic_q=ridge_periodic_q,
                                                                  pixel_adjacency=pixel_adjacency,
                                                                  relative_score_floor=ridge_relative_score_floor,
                                                                  min_branch_separation_bins=ridge_min_branch_separation_bins,
                                                                  max_branch_separation_frac=ridge_max_branch_separation_frac,
                                                                  min_cycle_q_fraction=ridge_min_cycle_q_fraction,
                                                                  min_two_peak_fraction=ridge_min_two_peak_fraction,
                                                                  min_cycle_confidence=ridge_min_cycle_confidence) :
        RidgeTopologyResult(RidgeNode[], Tuple{Int,Int}[], 0, 0, Int[], falses(nq, nω), zeros(Int32, nq, nω), 0, RidgeCycleCandidate[], 0, :disabled, Dict{Symbol,Any}())
    if verbose
        @info "resolution-aware smoothing: radius_q=$radius_q, radius_ω=$support_radius_ω, Γq=$(round(Γq; digits=5)), Γω=$(round(Γω; digits=5))"
        @info "soft adaptive support: robust background=$(round(bg; digits=5)), robust σ=$(round(σ; digits=5)), min component pixels=$minpix"
        @info "centerline support: added pixels=$centerline_pixels, anchored q-rows=$centerline_rows"
        @info "continuity repair: bridge pixels=$bridge_pixels, repaired q-rows=$bridge_rows"
        @info "validation: parent bands=$(qfirst.nparents)"
        @info "validation: peak counts before tracking=$(sort(unique(qfirst.peak_counts)))"
        @info "validation: kept components per parent=$component_counts"
        @info "validation: total final components=$(length(components))"
        @info "ridge topology: nodes=$(length(ridge_topology_result.nodes)), edges=$(length(ridge_topology_result.edges)), beta0=$(ridge_topology_result.beta0), beta1=$(ridge_topology_result.beta1), raw_beta1=$(get(ridge_topology_result.diagnostics, :raw_beta1, 0)), accepted_cycles=$(ridge_topology_result.accepted_cycle_count), status=$(ridge_topology_result.status), peak_counts=$(sort(unique(ridge_topology_result.peak_counts))), ridge_smooth_radius_ω=$r_smooth, effective_max_peaks=$effective_ridge_max_peaks"
        for c in components
            @info "component parent=$(c.parent_band), local=$(c.component): qbox=$(c.qbox), ωbox=$(c.ωbox), pixels=$(c.npixels), weight=$(round(c.weight; digits=4))"
        end
    end

    return (qfirst=qfirst, components=components, component_counts=component_counts,
            raw_component_counts=raw_counts, smooth=S,
            ridge_topology=ridge_topology_result,
            ridge_beta0=ridge_topology_result.beta0,
            ridge_beta1=ridge_topology_result.beta1,
            ridge_peak_counts=ridge_topology_result.peak_counts,
            ridge_mask=ridge_topology_result.ridge_mask,
            ridge_labels=ridge_topology_result.ridge_labels,
            ridge_cycle_candidates=ridge_topology_result.cycle_candidates,
            ridge_cycle_count=ridge_topology_result.accepted_cycle_count,
            ridge_status=ridge_topology_result.status,
            ridge_diagnostics=ridge_topology_result.diagnostics,
            centerline_pixels=centerline_pixels,
            centerline_rows=centerline_rows,
            bridge_pixels=bridge_pixels,
            bridge_rows=bridge_rows,
            background=bg, sigma=σ,
            radius_q=radius_q, radius_ω=support_radius_ω,
            estimated_radius_ω=radius_ω_est,
            peak_radius_ω=peak_radius_ω,
            resolution_q=Γq, resolution_ω=Γω)
end

#-----------------#
# User-Facing API #
#-----------------#
"""
    components_to_band_marginals(result)

Convert an `adaptive_autosplitbands2d` result into band marginal arrays.

Returns `(bandIq, bandIω)`, where each row corresponds to one retained
connected component. `bandIq` is integrated over energy, and `bandIω` is
integrated over momentum.
"""
function components_to_band_marginals(result)
    if isempty(result.components)
        @warn "splitbands: No components found during band-segmentation!"
        _, nω = size(result.smooth)
        nq = size(result.smooth, 1)
        return (bandIq = zeros(Float64, 0, nq), bandIω = zeros(Float64, 0, nω))
    end
    ncomp = length(result.components)
    nq = length(result.components[1].bandIq)
    nω = length(result.components[1].bandIω)

    bandIq = zeros(Float64, ncomp, nq)
    bandIω = zeros(Float64, ncomp, nω)

    for (i, c) in enumerate(result.components)
        bandIq[i, :] .= c.bandIq
        bandIω[i, :] .= c.bandIω
    end

    return (bandIq = bandIq, bandIω = bandIω)
end
"""
    splitbands(Iqω, qs, ωs; kwargs...)

User-facing band-segmentation API. Runs `adaptive_autosplitbands2d` and returns
only the component marginal arrays `(bandIq, bandIω)`.

Use this when the downstream workflow expects marginal spectra rather than the
full segmentation and ridge-topology diagnostics.
"""
function splitbands(Iqω::AbstractMatrix{<:Real}, qs, ωs; kwargs...)
    result = adaptive_autosplitbands2d(Iqω, qs, ωs; kwargs...)
    return components_to_band_marginals(result)
end
