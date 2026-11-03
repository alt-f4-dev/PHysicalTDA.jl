#-------------------------------------#
# Child-Ridge Recovery                #
#-------------------------------------#
# Recovers the individual ridges (child dispersions) inside connected parent supports: local
# peak candidates and their energy supports, parents as connected supports across momentum,
# a persistent child count per parent, tangent-aware tracking of each child, and the resulting
# ridge graph with split, merge and crossing points.
#
# Ported from the manuscript script ridge-dispersion-recovery.jl. With the default
# RidgeRecoveryParams() the results are identical to that script.

"""
    RidgeRecoveryParams(; kwargs...)
    RidgeRecoveryParams(ωs; max_slope, kwargs...)

Settings for [`find_ridge_candidates`](@ref) and [`recover_child_ridges`](@ref). Energies are in the
units of `ωs`, momenta in the units of `qs`, and slopes in energy per momentum unit.

The keyword constructor's defaults are the manuscript values (powder SEQUOIA data, 0.1 meV
energy bins).

The second constructor is a preset for coarse, resolution-matched grids, written in units of the
energy bin `Δω = mean(diff(ωs))`: half-maximum supports (`support_fraction = 0.5`),
`support_max_half_width = 2Δω`, `support_touch = Δω`, `prediction_scale = Δω`,
`fallback_window = Δω`, and `min_child_q_fraction = 0.15`. `max_slope` (the largest expected
`|dω/dq|`, with some margin) is required. Any keyword overrides the preset. This preset has been
checked only on simulated single-crystal data (Rb₂MnF₄ along (H, H, 0), 0.4 meV bins).

Fields:

| Field | Default | Meaning |
|---|---|---|
| `peak_smooth_radius` | 1 | energy smoothing radius (bins) for peak candidates |
| `min_prominence` | 0.005 | absolute prominence of a candidate peak |
| `relative_score_floor` | 0.0035 | candidate score floor relative to the strongest peak |
| `max_peaks_per_slice` | 16 | candidates kept per momentum slice |
| `support_fraction` | 0.95 | a support extends while the normalized intensity is ≥ this fraction of the peak |
| `support_min_level` | 0.001 | absolute floor of that threshold |
| `support_max_half_width` | 0.22 | largest support half-width (energy) |
| `support_touch` | 0.020 | supports closer than this (energy) are joined, within and across slices |
| `min_parent_q_fraction` | 0.035 | parents spanning a smaller fraction of the momentum axis are dropped |
| `min_child_q_fraction` | 0.025 | fraction of a parent's slices with ≥ k supports needed to count k children |
| `max_children` | 6 | largest number of children per parent |
| `extra_assignment_candidates` | 2 | candidates beyond the child count considered per slice |
| `max_slope` | 16 | slope scale (energy per momentum unit); larger slopes are penalized |
| `prediction_scale` | 0.05 | energy scale of the deviation from the predicted position |
| `curvature_weight` | 0.65 | weight of changes in slope |
| `score_weight` | 0.085 | weight of candidate intensity |
| `resolved_velocity_memory` | 0.35 | slope memory when all children are resolved |
| `unresolved_velocity_memory` | 0.97 | slope memory when children overlap |
| `fallback_window` | 0.05 | energy window searched when a child has no candidate |
| `fallback_shape_weight` | 0.50 | penalty on departing from the prediction inside that window |
"""
Base.@kwdef struct RidgeRecoveryParams
    peak_smooth_radius::Int = 1
    min_prominence::Float64 = 0.005
    relative_score_floor::Float64 = 0.0035
    max_peaks_per_slice::Int = 16
    support_fraction::Float32 = 0.95f0
    support_min_level::Float32 = 0.001f0
    support_max_half_width::Float32 = 0.22f0
    support_touch::Float32 = 0.020f0
    min_parent_q_fraction::Float32 = 0.035f0
    min_child_q_fraction::Float32 = 0.025f0
    max_children::Int = 6
    extra_assignment_candidates::Int = 2
    max_slope::Float32 = 16.0f0
    prediction_scale::Float32 = 0.05f0
    curvature_weight::Float32 = 0.65f0
    score_weight::Float32 = 0.085f0
    resolved_velocity_memory::Float32 = 0.35f0
    unresolved_velocity_memory::Float32 = 0.97f0
    fallback_window::Float32 = 0.05f0
    fallback_shape_weight::Float32 = 0.50f0
end

function RidgeRecoveryParams(ωs::AbstractVector{<:Real}; max_slope::Real, kwargs...)
    Δω = Float32(mean(abs.(diff(Float64.(ωs)))))
    preset = (support_fraction=0.5f0, support_max_half_width=2Δω, support_touch=Δω,
              prediction_scale=Δω, fallback_window=Δω, min_child_q_fraction=0.15f0,
              max_slope=Float32(max_slope))
    return RidgeRecoveryParams(; preset..., kwargs...)
end

"""
    RidgeSupportGroup

Energy support `lo:hi` (bin indices) at momentum row `iq`, and the candidate node ids it contains.
"""
struct RidgeSupportGroup
    iq::Int32
    lo::Int32
    hi::Int32
    candidates::Vector{Int32}
end

"""
    RidgeCandidates

Result of [`find_ridge_candidates`](@ref).

Fields: `I` (input, `Float32`), `A` (energy-smoothed map normalized to [0, 1] per momentum row),
`q`, `ω`, `nodes` (candidate peaks as `RidgeNode`s), `byq` (node ids per row), `groups`
([`RidgeSupportGroup`](@ref)s), `groups_by_q` (group ids per row), and `params`.
"""
struct RidgeCandidates
    I::Matrix{Float32}
    A::Matrix{Float32}
    q::Vector{Float32}
    ω::Vector{Float32}
    nodes::Vector{RidgeNode}
    byq::Vector{Vector{Int32}}
    groups::Vector{RidgeSupportGroup}
    groups_by_q::Vector{Vector{Int32}}
    params::RidgeRecoveryParams
end

"""
    ChildRidgeResult

Result of [`recover_child_ridges`](@ref).

Fields:

- `candidates` :: the [`RidgeCandidates`](@ref) used;
- `group_parent` :: parent index of each support group (0 when its parent was dropped);
- `parent_groups` :: support-group ids of each parent, ordered by lowest energy;
- `local_multiplicity` :: `nparents × nq`, number of support groups of each parent per row;
- `global_children` :: number of children tracked in each parent;
- `tracks` :: per parent, a `children × nq` matrix of child energies (`NaN` outside the parent);
- `support_mask` :: union of all parent supports, on the candidate grid;
- `graph_nodes`, `graph_edges` :: ridge graph built from the tracks;
- `split_nodes`, `merge_nodes` :: graph nodes with out-degree > 1 (splits) or in-degree > 1 (merges);
- `crossing_points` :: `(q, ω)` where two children of a parent exchange order.
"""
struct ChildRidgeResult
    candidates::RidgeCandidates
    group_parent::Vector{Int32}
    parent_groups::Vector{Vector{Int32}}
    local_multiplicity::Matrix{Int16}
    global_children::Vector{Int16}
    tracks::Vector{Matrix{Float32}}
    support_mask::BitMatrix
    graph_nodes::Vector{RidgeNode}
    graph_edges::Vector{Tuple{Int,Int}}
    split_nodes::Vector{Int}
    merge_nodes::Vector{Int}
    crossing_points::Vector{NTuple{2,Float32}}
end

function Base.show(io::IO, r::ChildRidgeResult)
    print(io, "ChildRidgeResult($(length(r.parent_groups)) parents, children=$(Int.(r.global_children)), ",
          "$(length(r.split_nodes)) splits, $(length(r.merge_nodes)) merges, $(length(r.crossing_points)) crossings)")
end

# Energy smoothing (box filter) and per-row normalization to [0, 1].
function _ridge_smooth_normalize(I, r)
    nq, nω = size(I)
    A = Matrix{Float32}(undef, nq, nω)
    Threads.@threads for iq in 1:nq
        row = Vector{Float64}(undef, nω)
        moving_average_1d!(row, @view(I[iq, :]), r)
        normalize01!(row)
        @inbounds for j in 1:nω
            A[iq, j] = Float32(row[j])
        end
    end
    return A
end

"""
    find_ridge_candidates(I, qs, ωs; params=RidgeRecoveryParams(), periodic_q=false) -> RidgeCandidates

Local peak candidates of each momentum row of `I` (`nq × nω`), found with
[`build_ridge_topology`](@ref), and their energy supports. A support extends from a peak while
the normalized intensity stays above `max(support_min_level, support_fraction × peak)`, up to
`support_max_half_width`; supports within `support_touch` of each other are joined into one
[`RidgeSupportGroup`](@ref).
"""
function find_ridge_candidates(Iqω::AbstractMatrix{<:Real}, qs, ωs;
                               params::RidgeRecoveryParams=RidgeRecoveryParams(), periodic_q::Bool=false)
    p = params
    size(Iqω) == (length(qs), length(ωs)) ||
        throw(DimensionMismatch("I is $(size(Iqω)); qs and ωs give $((length(qs), length(ωs)))."))
    I = Matrix{Float32}(Iqω); q = Float32.(qs); ω = Float32.(ωs)
    A = _ridge_smooth_normalize(I, p.peak_smooth_radius)
    seed = build_ridge_topology(I, q, ω; smooth_radius_ω=p.peak_smooth_radius,
                                min_prominence=p.min_prominence, max_peaks_per_slice=p.max_peaks_per_slice,
                                relative_score_floor=p.relative_score_floor, max_jump_bins=nothing, tube_radius_ω=1,
                                periodic_q=periodic_q, min_cycle_confidence=1.0)
    nodes = seed.nodes; nq, nω = size(I)
    byq = [Int32[] for _ in 1:nq]
    @inbounds for (id, n) in pairs(nodes)
        push!(byq[n.iq], Int32(id))
    end
    groups = RidgeSupportGroup[]; groups_by_q = [Int32[] for _ in 1:nq]
    Δω = mean(diff(ω))
    maxr = max(1, round(Int, p.support_max_half_width / Δω))
    touch = max(0, round(Int, p.support_touch / Δω))
    for iq in 1:nq
        ids = sort(byq[iq]; by=id -> nodes[id].iω)
        isempty(ids) && continue
        intervals = Tuple{Int,Int,Int32}[]
        @inbounds for id in ids
            j = nodes[id].iω
            τ = max(p.support_min_level, p.support_fraction * A[iq, j])
            lo = j; hi = j
            while lo > 1 && j - lo < maxr && A[iq, lo-1] >= τ; lo -= 1; end
            while hi < nω && hi - j < maxr && A[iq, hi+1] >= τ; hi += 1; end
            push!(intervals, (lo, hi, id))
        end
        lo, hi, id = intervals[1]; members = Int32[id]
        for k in 2:length(intervals)
            a, b, c = intervals[k]
            if a <= hi + touch
                hi = max(hi, b); push!(members, c)
            else
                push!(groups, RidgeSupportGroup(Int32(iq), Int32(lo), Int32(hi), members))
                push!(groups_by_q[iq], Int32(length(groups)))
                lo, hi, members = a, b, Int32[c]
            end
        end
        push!(groups, RidgeSupportGroup(Int32(iq), Int32(lo), Int32(hi), members))
        push!(groups_by_q[iq], Int32(length(groups)))
    end
    return RidgeCandidates(I, A, q, ω, nodes, byq, groups, groups_by_q, p)
end

# Best support-guided continuation near a predicted energy, for a child with no candidate.
function _ridge_support_candidate(c, parent_of_group, parent, iq, pred)
    p = c.params
    Δω = mean(diff(c.ω))
    j0 = clamp(round(Int, 1 + (pred - c.ω[1]) / Δω), 1, length(c.ω))
    r = max(1, round(Int, p.fallback_window / Δω))
    best = -Inf32; bestj = j0
    @inbounds for gid in c.groups_by_q[iq]
        parent_of_group[gid] == parent || continue
        g = c.groups[gid]
        lo = max(Int(g.lo), j0 - r); hi = min(Int(g.hi), j0 + r)
        lo <= hi || continue
        for j in lo:hi
            z = (c.ω[j] - pred) / p.fallback_window
            value = c.A[iq, j] - p.fallback_shape_weight * z * z
            value > best && (best = value; bestj = j)
        end
    end
    return c.ω[bestj], max(best, 0f0)
end

# Exhaustive tangent-aware assignment of m children to n candidates.
function _ridge_best_assignment(p, pred, prev, vel, E, S, dq, unique_candidates)
    m = length(pred); n = length(E)
    best = fill(1, m); trial = similar(best); bestcost = Inf32; total = n^m
    @inbounds for code in 0:total-1
        x = code; used = UInt64(0); cost = 0f0; valid = true
        for b in 1:m
            k = x % n + 1; x ÷= n; trial[b] = k
            if unique_candidates
                bit = UInt64(1) << (k - 1)
                (used & bit) != 0 && (valid = false; break)
                used |= bit
            end
            slope = (E[k] - prev[b]) / dq
            de = (E[k] - pred[b]) / p.prediction_scale
            dv = (slope - vel[b]) / p.max_slope
            cost += de * de + p.curvature_weight * dv * dv - p.score_weight * S[k]
            abs(slope) > p.max_slope && (cost += 4f0 * (abs(slope) / p.max_slope - 1f0)^2)
        end
        valid && cost < bestcost && (bestcost = cost; copyto!(best, trial))
    end
    return best
end

# Tracks m children of one parent across momentum.
function _ridge_track_parent(c, parent_of_group, parent, groups, m)
    p = c.params
    nq = length(c.q); ids_by_q = [Int32[] for _ in 1:nq]
    for gid in groups
        append!(ids_by_q[c.groups[gid].iq], c.groups[gid].candidates)
    end
    qinds = findall(x -> !isempty(x), ids_by_q); qlo, qhi = first(qinds), last(qinds)
    T = fill(Float32(NaN), m, nq); V = zeros(Float32, m)
    ids = ids_by_q[qlo]
    order = sort(ids; by=id -> c.nodes[id].score, rev=true)[1:min(length(ids), m)]
    E = sort(Float32[c.nodes[id].ω for id in order])
    @inbounds for b in 1:m
        T[b, qlo] = E[clamp(round(Int, 1 + (b - 1) * (length(E) - 1) / max(m - 1, 1)), 1, length(E))]
    end
    for iq in qlo+1:qhi
        dq = max(c.q[iq] - c.q[iq-1], eps(Float32)); prev = @view T[:, iq-1]
        pred = Float32[prev[b] + V[b] * dq for b in 1:m]
        ids = ids_by_q[iq]; nobs = length(ids)
        if nobs > m + p.extra_assignment_candidates
            ids = sort(ids; by=id -> c.nodes[id].score, rev=true)[1:m+p.extra_assignment_candidates]
        end
        E = Float32[c.nodes[id].ω for id in ids]; S = Float32[c.nodes[id].score for id in ids]
        if nobs < m || isempty(E)
            for b in 1:m
                e, s = _ridge_support_candidate(c, parent_of_group, parent, iq, pred[b])
                push!(E, e); push!(S, 0.5f0 * s)
            end
        end
        assign = _ridge_best_assignment(p, pred, prev, V, E, S, dq, length(E) >= m)
        memory = nobs < m ? p.unresolved_velocity_memory : p.resolved_velocity_memory
        @inbounds for b in 1:m
            e = E[assign[b]]; slope = (e - prev[b]) / dq
            T[b, iq] = e; V[b] = memory * V[b] + (1f0 - memory) * slope
        end
    end
    return T
end

"""
    recover_child_ridges(candidates::RidgeCandidates) -> ChildRidgeResult
    recover_child_ridges(I, qs, ωs; params=RidgeRecoveryParams(), periodic_q=false) -> ChildRidgeResult

Child ridges inside connected parent supports.

1. Support groups that overlap (within `support_touch`) in neighbouring momentum rows are joined
   into parents; parents spanning less than `min_parent_q_fraction` of the momentum axis are dropped.
2. A parent has `k` children when at least `min_child_q_fraction` of its rows contain ≥ `k`
   support groups (at most `max_children`).
3. The children are tracked across momentum by assigning candidates with a cost on the deviation
   from the predicted energy, the change in slope, slopes above `max_slope`, and candidate
   intensity. Where fewer candidates than children are found, the best point inside the parent's
   support near the prediction is used.
4. The tracks form a ridge graph whose split, merge and crossing points are reported.

Use [`child_ridge_dispersions`](@ref) for the tracks as variable-length segments. The second form
runs [`find_ridge_candidates`](@ref) first.
"""
function recover_child_ridges(Iqω::AbstractMatrix{<:Real}, qs, ωs;
                              params::RidgeRecoveryParams=RidgeRecoveryParams(), periodic_q::Bool=false)
    return recover_child_ridges(find_ridge_candidates(Iqω, qs, ωs; params=params, periodic_q=periodic_q))
end

function recover_child_ridges(c::RidgeCandidates)
    p = c.params
    ng = length(c.groups); par = collect(1:ng); rank = zeros(Int, ng)
    touch = max(0, round(Int, p.support_touch / mean(diff(c.ω))))
    for iq in 1:length(c.q)-1
        A = c.groups_by_q[iq]; B = c.groups_by_q[iq+1]; ia = 1; ib = 1
        while ia <= length(A) && ib <= length(B)
            ga = c.groups[A[ia]]; gb = c.groups[B[ib]]
            if Int(ga.hi) + touch < Int(gb.lo)
                ia += 1
            elseif Int(gb.hi) + touch < Int(ga.lo)
                ib += 1
            else
                _uf_union!(par, rank, Int(A[ia]), Int(B[ib]))
                Int(ga.hi) <= Int(gb.hi) ? (ia += 1) : (ib += 1)
            end
        end
    end
    roots = Dict{Int,Vector{Int32}}()
    for gid in 1:ng
        push!(get!(roots, _uf_find!(par, gid), Int32[]), Int32(gid))
    end
    parent_groups = Vector{Vector{Int32}}()
    for gs in values(roots)
        iqs = Int[c.groups[g].iq for g in gs]
        (maximum(iqs) - minimum(iqs) + 1) / length(c.q) >= p.min_parent_q_fraction && push!(parent_groups, gs)
    end
    sort!(parent_groups; by=gs -> minimum(c.groups[g].lo for g in gs))
    np = length(parent_groups); group_parent = zeros(Int32, ng)
    for (ip, gs) in pairs(parent_groups), g in gs
        group_parent[g] = Int32(ip)
    end
    local_multiplicity = zeros(Int16, np, length(c.q)); support = falses(size(c.I))
    children = zeros(Int16, np); tracks = Vector{Matrix{Float32}}(undef, np)
    for (ip, gs) in pairs(parent_groups)
        for gid in gs
            g = c.groups[gid]
            local_multiplicity[ip, g.iq] += Int16(1)
            @views support[g.iq, g.lo:g.hi] .= true
        end
        qactive = findall(>(0), @view local_multiplicity[ip, :]); span = length(qactive)
        m = 1
        for k in 2:p.max_children
            count(x -> x >= k, @view local_multiplicity[ip, qactive]) / span >= p.min_child_q_fraction && (m = k)
        end
        children[ip] = Int16(m)
        tracks[ip] = _ridge_track_parent(c, group_parent, Int32(ip), gs, m)
    end
    graph_nodes = RidgeNode[]; edges = Set{Tuple{Int,Int}}(); key = Dict{Tuple{Int,Int},Int}()
    Δω = mean(diff(c.ω))
    for ip in 1:np, b in axes(tracks[ip], 1)
        lastid = 0
        for iq in axes(tracks[ip], 2)
            e = tracks[ip][b, iq]; isfinite(e) || continue
            j = clamp(round(Int, 1 + (e - c.ω[1]) / Δω), 1, length(c.ω))
            id = get!(key, (iq, j)) do
                push!(graph_nodes, RidgeNode(iq, j, Float64(c.q[iq]), Float64(e), 1.0))
                length(graph_nodes)
            end
            lastid > 0 && push!(edges, (lastid, id)); lastid = id
        end
    end
    edgevec = sort!(collect(edges))
    indeg = zeros(Int, length(graph_nodes)); outdeg = zeros(Int, length(graph_nodes))
    for (a, b) in edgevec
        outdeg[a] += 1; indeg[b] += 1
    end
    splits = findall(i -> outdeg[i] > 1 && indeg[i] <= 1, eachindex(graph_nodes))
    merges = findall(i -> indeg[i] > 1 && outdeg[i] <= 1, eachindex(graph_nodes))
    crossings = NTuple{2,Float32}[]
    for ip in 1:np, a in 1:(size(tracks[ip], 1)-1), b in (a+1):size(tracks[ip], 1), iq in 1:(length(c.q)-1)
        d1 = tracks[ip][a, iq] - tracks[ip][b, iq]; d2 = tracks[ip][a, iq+1] - tracks[ip][b, iq+1]
        if isfinite(d1) && isfinite(d2) && d1 * d2 < 0
            t = abs(d1) / (abs(d1) + abs(d2))
            push!(crossings, ((1f0 - t) * c.q[iq] + t * c.q[iq+1],
                              0.5f0 * ((1f0 - t) * (tracks[ip][a, iq] + tracks[ip][b, iq]) +
                                       t * (tracks[ip][a, iq+1] + tracks[ip][b, iq+1]))))
        end
    end
    return ChildRidgeResult(c, group_parent, parent_groups, local_multiplicity, children, tracks, support,
                            graph_nodes, edgevec, splits, merges, crossings)
end

"""
    child_ridge_dispersions(result::ChildRidgeResult)

The child tracks of [`recover_child_ridges`](@ref) as variable-length segments.

Returns `(segments, global_children, local_multiplicity, split_nodes, merge_nodes, crossing_points)`,
where each segment is `(parent, child, q, ω, local_multiplicity)` over the rows the child spans.
"""
function child_ridge_dispersions(t::ChildRidgeResult)
    segs = NamedTuple[]
    for ip in eachindex(t.tracks), b in axes(t.tracks[ip], 1)
        inds = findall(isfinite, @view(t.tracks[ip][b, :]))
        isempty(inds) && continue
        push!(segs, (parent=ip, child=b, q=t.candidates.q[inds], ω=copy(@view(t.tracks[ip][b, inds])),
                     local_multiplicity=copy(@view(t.local_multiplicity[ip, inds]))))
    end
    return (segments=segs, global_children=t.global_children, local_multiplicity=t.local_multiplicity,
            split_nodes=t.split_nodes, merge_nodes=t.merge_nodes, crossing_points=t.crossing_points)
end
