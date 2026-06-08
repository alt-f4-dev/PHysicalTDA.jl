#--------------------------------------------------------------------
# Complex.jl
#
# Periodic cubical persistent homology via Ripserer.Custom.
#
# PeriodicCubical(data; periodic_dims, threshold) builds an explicit
# CW-complex with correct periodic identification and returns it as a
# Ripserer.Custom filtration.  Standard cohomology then gives correct
# essential H₁, H₂, etc. for cylinders, tori, etc.
#
# Public API:
#   PeriodicCubical(data; periodic_dims, threshold) -> Custom filtration
#   pd_array_intensities(Z; periodic_dims=nothing, ...)
#
# Dependencies: Ripserer (loaded by parent module)
#--------------------------------------------------------------------

#------------------------------------------------#
#  Section 1: parse periodic_dims specification  #
#------------------------------------------------#

"""
    _parse_periodic_dims(K, periodic_dims) -> NTuple{K,Bool}

Accept either an NTuple{K,Bool} mask or an iterable of 1-based axis indices.
"""
function _parse_periodic_dims(K::Int, periodic_dims::NTuple{<:Any,Bool})
    return periodic_dims
end

function _parse_periodic_dims(K::Int, periodic_dims)
    return ntuple(i -> i in periodic_dims, K)
end

#-----------------------------------------------------------------#
#  Section 2: _tile_periodic — kept for PHysicalTDA test access   #
#-----------------------------------------------------------------#

"""
    _tile_periodic(data, periodic) -> Array

Return a copy of `data` with each periodic axis extended by one slice
(first slice repeated at the end).
"""
function _tile_periodic(
    data::AbstractArray{T,K},
    periodic::NTuple{K,Bool},
) where {T,K}
    new_size = ntuple(i -> periodic[i] ? size(data, i) + 1 : size(data, i), Val(K))
    result   = Array{T,K}(undef, new_size)
    orig_ranges = ntuple(i -> 1:size(data, i), Val(K))
    result[orig_ranges...] = data
    for ax in 1:K
        periodic[ax] || continue
        n   = size(data, ax)
        src = ntuple(i -> i == ax ? (1:1)     : (1:size(result, i)), Val(K))
        dst = ntuple(i -> i == ax ? (n+1:n+1) : (1:size(result, i)), Val(K))
        result[dst...] = result[src...]
    end
    return result
end

#--------------------------------------------------------------#
#  Section 3: _periodic_cubemap — kept for test inspection     #
#--------------------------------------------------------------#

@inline function _periodic_cubemap_size(
    sz::NTuple{K,Int}, periodic::NTuple{K,Bool}
) where K
    ntuple(i -> periodic[i] ? 2 * sz[i] : 2 * sz[i] - 1, Val(K))
end

"""
    _periodic_cubemap(data, periodic) -> Array

Build the periodic cubemap.  Retained for test inspection only.
"""
function _periodic_cubemap(
    data::AbstractArray{T,K},
    periodic::NTuple{K,Bool},
) where {T,K}
    csz    = _periodic_cubemap_size(size(data), periodic)
    result = fill(typemin(T), csz)
    for ci in CartesianIndices(data)
        _update_cubemap_rec!(result, data[ci], Tuple(ci), periodic, csz, 1,
                             Vector{Int}(undef, K))
    end
    return result
end

function _update_cubemap_rec!(result, v, t, periodic, csz, dim, idx)
    K    = length(csz)
    mid  = 2 * t[dim] - 1
    lo   = mid - 1
    hi   = mid + 1
    lo_c = periodic[dim] ? mod1(lo == 0 ? csz[dim] : lo, csz[dim]) : max(lo, 1)
    hi_c = periodic[dim] ? mod1(hi, csz[dim])                       : min(hi, csz[dim])
    for c in (lo_c, mid, hi_c)
        idx[dim] = c
        if dim == K
            ci = CartesianIndex(NTuple{length(idx),Int}(idx))
            checkbounds(Bool, result, ci) && (result[ci] = max(result[ci], v))
        else
            _update_cubemap_rec!(result, v, t, periodic, csz, dim + 1, idx)
        end
    end
end

#------------------------------------------------------------#
#  Section 4: PeriodicCubical constructor → Custom          #
#------------------------------------------------------------#

"""
    PeriodicCubical(data; periodic_dims, threshold=maximum(data))

Build a periodic cubical filtration and return it as a `Ripserer.Custom`
filtration with explicit periodic identification.

Periodic axes are closed by wrapping vertex n+1 ≡ vertex 1, producing
the correct topology (cylinder, torus, etc.) which Ripserer's standard
cohomology algorithm reduces correctly.

- `data`: `AbstractArray{T,K}` of filtration values.
- `periodic_dims`: `NTuple{K,Bool}` mask or iterable of 1-based axis
  indices that are periodic.
- `threshold`: cells born above this value are excluded.

# Example

```julia
img  = rand(Float64, 50, 200)
cf   = PeriodicCubical(img; periodic_dims=(true, false))
dgms = ripserer(cf; dim_max=1)
```
"""
function PeriodicCubical(
    data::AbstractArray{T,K};
    periodic_dims,
    threshold = maximum(data),
) where {T,K}
    periodic = _parse_periodic_dims(K, periodic_dims)
    τ        = T(threshold)
    sz       = size(data)
    lin      = LinearIndices(data)

    # Map a CartesianIndex (possibly exceeding bounds on periodic axes) to
    # the canonical linear vertex index in the original data array.
    function canonical(ci::CartesianIndex{K})
        t = ntuple(Val(K)) do i
            periodic[i] ? mod1(Tuple(ci)[i], sz[i]) : Tuple(ci)[i]
        end
        return lin[CartesianIndex{K}(t)]
    end

    simplices = Pair{Vector{Int}, T}[]

    # --- vertices (0-cubes) ---
    for ci in CartesianIndices(data)
        data[ci] > τ && continue
        push!(simplices, [lin[ci]] => data[ci])
    end

    # --- d-cubes for d = 1 … K ---
    # Enumerate all size-d subsets of axes, then for each corner in data,
    # build the 2^d-vertex cube extending +1 along those axes.
    axis_subsets = _all_axis_subsets(K)   # pre-computed for all d

    for d in 1:K
        for axes_combo in axis_subsets[d]
            for ci in CartesianIndices(data)
                t = Tuple(ci)

                # Build 2^d corners; return nothing if aperiodic axis
                # would go out of range.
                corners = _cube_corners(t, axes_combo, sz, periodic, K)
                isnothing(corners) && continue

                birth_val = maximum(data[c] for c in corners)
                birth_val > τ && continue

                verts = sort!([canonical(c) for c in corners])
                push!(simplices, verts => birth_val)
            end
        end
    end

    return Custom(simplices; threshold = τ)
end

# Pre-compute all non-empty subsets of {1,…,K} grouped by size.
# Returns a Vector of length K where element d contains all size-d subsets.
function _all_axis_subsets(K::Int)
    result = [Vector{Int}[] for _ in 1:K]
    # iterate over all 2^K - 1 non-empty bitmasks
    for mask in 1:(2^K - 1)
        axes = [i for i in 1:K if (mask >> (i-1)) & 1 == 1]
        d    = length(axes)
        push!(result[d], axes)
    end
    return result
end

# Build the 2^d corners of a d-cube anchored at t, extended +1 along
# each axis in axes_combo.  Returns nothing if any corner is out of
# range on an aperiodic axis.
function _cube_corners(
    t::NTuple{K,Int},
    axes_combo::Vector{Int},
    sz::NTuple{K,Int},
    periodic::NTuple{K,Bool},
    ::Int,
) where K
    d = length(axes_combo)
    corners = CartesianIndex{K}[]
    for bits in 0:(2^d - 1)
        nt    = collect(t)
        valid = true
        for (k, ax) in enumerate(axes_combo)
            if (bits >> (k - 1)) & 1 == 1
                nt[ax] += 1
                if !periodic[ax] && nt[ax] > sz[ax]
                    valid = false
                    break
                end
            end
        end
        valid || return nothing
        push!(corners, CartesianIndex{K}(Tuple(nt)))
    end
    return corners
end

#----------------------------------------------#
#  Section 5: PHysicalTDA.jl integration       #
#----------------------------------------------#

"""
    pd_array_intensities(Z::AbstractMatrix{<:Real}; kwargs...)

Compute persistent homology of the 2D intensity array `Z`.

# Keywords
- `maxdim::Int = 1`: maximum homological dimension.
- `threshold::Real`: filtration cutoff (default: maximum of `Z`).
- `superlevel::Bool = false`: negate `Z` for superlevel filtration.
- `normalize::Bool = false`: scale `Z` to [0,1] before filtration.
- `periodic_dims`: `nothing` (default, standard `Cubical`), `NTuple{K,Bool}`
  mask, or iterable of periodic axis indices.

# Returns
`Vector{PersistenceDiagram}` of length `maxdim + 1`.
"""
function pd_array_intensities(
    Z::AbstractMatrix{<:Real};
    maxdim::Int      = 1,
    threshold        = nothing,
    superlevel::Bool = false,
    normalize::Bool  = false,
    periodic_dims    = nothing,
)
    arr = Float64.(Z)

    if normalize
        lo, hi = extrema(arr)
        span   = hi - lo
        arr    = span > 0 ? (arr .- lo) ./ span : arr .- lo
    end

    superlevel && (arr = -arr)

    τ  = isnothing(threshold) ? maximum(arr) : Float64(threshold)
    cf = isnothing(periodic_dims) ? Cubical(arr; threshold = τ) :
                                    PeriodicCubical(arr; periodic_dims = periodic_dims,
                                                        threshold = τ)

    return collect(ripserer(cf; dim_max = maxdim))
end

# N-dimensional overload
function pd_array_intensities(
    Z::AbstractArray{<:Real,N};
    maxdim::Int      = 1,
    threshold        = nothing,
    superlevel::Bool = false,
    normalize::Bool  = false,
    periodic_dims    = nothing,
) where N
    arr = Float64.(Z)

    if normalize
        lo, hi = extrema(arr)
        span   = hi - lo
        arr    = span > 0 ? (arr .- lo) ./ span : arr .- lo
    end

    superlevel && (arr = -arr)

    τ  = isnothing(threshold) ? maximum(arr) : Float64(threshold)
    cf = isnothing(periodic_dims) ? Cubical(arr; threshold = τ) :
                                    PeriodicCubical(arr; periodic_dims = periodic_dims,
                                                        threshold = τ)

    return collect(ripserer(cf; dim_max = maxdim))
end
