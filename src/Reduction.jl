const _CovarianceEntry = Tuple{Int,Int,Float64}
_numerator_covariance(h) = get(h.meta, :numerator_covariance, ())

function _record_covariance!(cov, variance, i, j, value)
    iszero(value) && return
    if i == j
        variance[i] += 2value
    else
        key = minmax(i, j)
        cov[key] = get(cov, key, 0.0) + value
    end
end

_covariance_entries(cov) = sort!(_CovarianceEntry[(i,j,v) for ((i,j),v) in cov if !iszero(v)])

function _reduced_histogram(h, D, D2, N, edges, centres, cov; names=h.names, units=h.units, frames=h.frames)
    return NormalizedHistogram(D, D2, N, edges, centres, names, units, frames,
                               merge(h.meta, (; numerator_covariance=cov)))
end

function _reduction_dims(dims, rank)
    ds = dims isa Integer ? (Int(dims),) : Tuple(Int(d) for d in dims)
    !isempty(ds) && allunique(ds) && all(d -> 1 <= d <= rank, ds) ||
        throw(ArgumentError("dims must specify distinct histogram dimensions"))
    return ds
end

_bin_value(x::Number, i) = x
_bin_value(x::AbstractArray, i) = x[i]
_bin_value(::Nothing, i) = 1.0
_bin_value(f, i) = f(i)

function _check_bin_array(x, h, name)
    x isa AbstractArray && axes(x) != axes(h.D) && throw(DimensionMismatch("$name must match histogram axes"))
end

function _accumulate_groups(D, D2::V, N, groups::NTuple{R}, weight::W, source_cov) where {V,R,W}
    shape = map(length, groups)
    sums, exposure = zeros(shape), zeros(shape)
    variance = D2 === nothing ? nothing : zeros(shape)
    scatter = isempty(source_cov) ? nothing : zeros(Int, size(D))
    weights = scatter === nothing ? nothing : zeros(size(D))
    target_ids = LinearIndices(shape)
    for k in CartesianIndices(sums)
        block = ntuple(d -> groups[d][k[d]], Val(R))
        numerator = denominator = var = 0.0
        for i in CartesianIndices(block)
            n = N[i]
            n > 0 || continue
            a = _bin_value(weight, i)
            isfinite(a) && a >= 0 || throw(ArgumentError("weights must be finite and nonnegative"))
            iszero(a) && continue
            numerator += a * D[i]
            denominator += a * n
            D2 === nothing || (var += abs2(a) * D2[i])
            if scatter !== nothing
                scatter[i] = target_ids[k]
                weights[i] = a
            end
        end
        sums[k], exposure[k] = numerator, denominator
        variance === nothing || (variance[k] = var)
    end
    scatter === nothing && return sums, variance, exposure, ()
    cov = Dict{Tuple{Int,Int},Float64}()
    for (i,j,v) in source_cov
        k,l = scatter[i], scatter[j]
        k > 0 && l > 0 && _record_covariance!(cov, variance, k, l, weights[i]*weights[j]*v)
    end
    return sums, variance, exposure, _covariance_entries(cov)
end

function _reduce_groups(h, groups, weight)
    _check_bin_array(weight, h, "weight")
    edges = [[h.edges[d][first(first(gs))]; [h.edges[d][last(g)+1] for g in gs]] for (d,gs) in enumerate(groups)]
    centres = [[length(g) == 1 ? h.centres[d][first(g)] : (h.edges[d][first(g)]+h.edges[d][last(g)+1])/2
                for g in gs] for (d,gs) in enumerate(groups)]
    D,D2,N,cov = _accumulate_groups(h.D, h.D2, h.N, groups, weight, _numerator_covariance(h))
    return _reduced_histogram(h, D, D2, N, edges, centres, cov)
end

"""
    rebin(h::NormalizedHistogram; factors, remainder=:drop)

Sum measured numerator, variance and exposure over contiguous blocks. Give one positive integer factor per axis or a common factor. Trailing bins can be dropped, merged into the last block, or kept as a short block using `:drop`, `:merge`, or `:keep`.
"""
function rebin(h::NormalizedHistogram; factors, remainder::Symbol=:drop)
    remainder in (:drop,:merge,:keep) || throw(ArgumentError("invalid remainder"))
    fs = factors isa Integer ? ntuple(_ -> factors, ndims(h.D)) : Tuple(factors)
    length(fs) == ndims(h.D) || throw(DimensionMismatch("supply one factor per axis"))
    all(f -> f isa Integer && f > 0, fs) || throw(ArgumentError("factors must be positive integers"))
    groups = ntuple(ndims(h.D)) do d
        n = size(h.D,d)
        f = Int(min(fs[d], n))
        m = remainder === :keep ? cld(n,f) : fld(n,f)
        [(i-1)*f+1:(remainder === :merge && i == m ? n : min(i*f,n)) for i in 1:m]
    end
    return _reduce_groups(h, groups, nothing)
end

function _integration_range(h,d,spec)
    spec isa Colon && return 1:size(h.D,d)
    if spec isa AbstractUnitRange{<:Integer}
        isempty(spec) && throw(ArgumentError("empty integration range"))
        checkbounds(h.centres[d], spec)
        return Int(first(spec)):Int(last(spec))
    end
    spec isa Tuple{Real,Real} && spec[1] <= spec[2] || throw(ArgumentError("expected a bin range or (lo, hi)"))
    c = h.centres[d]
    bins = searchsortedfirst(c,spec[1]):searchsortedlast(c,spec[2])
    isempty(bins) && throw(ArgumentError("integration interval has no bin centres"))
    return bins
end

"""
    integrate(h::NormalizedHistogram; dims, range=:, weight=nothing)

Integrate into singleton dimensions. A range can be contiguous bin indices or an inclusive physical interval `(lo, hi)` tested against bin centres. For several dimensions supply one selection per axis. Weights may be scalar, histogram-shaped, or callable on the Cartesian bin index. Variance uses squared weights.
"""
function integrate(h::NormalizedHistogram; dims, range=Colon(), weight=nothing)
    ds = _reduction_dims(dims,ndims(h.D))
    specs = if range isa Colon
        ntuple(_ -> range,length(ds))
    elseif length(ds) == 1 && !(range isa Tuple && length(range) == 1)
        (range,)
    else
        range
    end
    length(specs) == length(ds) || throw(DimensionMismatch("supply one range per integrated axis"))
    groups = ntuple(ndims(h.D)) do d
        j = findfirst(==(d),ds)
        j === nothing ? [i:i for i in axes(h.D,d)] : [_integration_range(h,d,specs[j])]
    end
    return _reduce_groups(h,groups,weight)
end

"""
    weighted_average(h::NormalizedHistogram, f; dims)

Average a scalar, histogram-shaped field, or Cartesian-index callable with MDNorm exposure weights. Keep reduced dimensions as singleton axes. Exclude nonfinite field values and bins without exposure.
"""
function weighted_average(h::NormalizedHistogram,f;dims)
    _check_bin_array(f,h,"f")
    ds = _reduction_dims(dims,ndims(h.D))
    return _weighted_average(h.N,f,ntuple(d -> d in ds,ndims(h.D)))
end

function _weighted_average(N::Array{Float64,R},f::F,collapsed::NTuple{R,Bool}) where {R,F}
    shape = ntuple(d -> collapsed[d] ? 1 : size(N,d),Val(R))
    sums, weights = zeros(shape), zeros(shape)
    for i in CartesianIndices(N)
        n = N[i]
        n > 0 || continue
        x = _bin_value(f,i)
        isfinite(x) || continue
        k = CartesianIndex(ntuple(d -> collapsed[d] ? 1 : i[d],Val(R)))
        sums[k] += n*x
        weights[k] += n
    end
    sums ./= weights
    return sums
end

"""
    rebin_to_resolution(h::NormalizedHistogram; factors=nothing, resolution_q=nothing, resolution_ω=nothing, remainder=:drop)

Use explicit block factors or determine them from the resolution of a histogram with one nonsingleton momentum axis and energy. Reduce the measured numerator and exposure before forming intensity.
"""
function rebin_to_resolution(h::NormalizedHistogram; factors=nothing,resolution_q=nothing,resolution_ω=nothing,
                             remainder::Symbol=:drop)
    factors === nothing || return rebin(h; factors,remainder)
    e = energy_axis(h)
    spatial = [d for d in 1:ndims(h.D) if d != e && size(h.D,d) > 1]
    length(spatial) == 1 && size(h.D,e) > 1 || throw(ArgumentError("supply factors for a multidimensional grid"))
    all(x -> x === nothing || (isfinite(x) && x>0), (resolution_q,resolution_ω)) ||
        throw(ArgumentError("invalid resolution width"))
    q = only(spatial)
    order = (q,e,(d for d in 1:ndims(h.D) if d != q && d != e)...)
    signal = reshape(permutedims(intensity(h),order),size(h.D,q),size(h.D,e))
    pos = minimum((x for x in signal if isfinite(x) && x>0);init=Inf)
    signal[.!isfinite.(signal)] .= isfinite(pos) ? pos : 0.0
    qs, es = h.centres[q], h.centres[e]
    _,_,Γq,Γe = adaptive_resolution_radii(qs,es,signal;resolution_q,resolution_ω)
    factor(Γ,axis) = clamp(floor(Int,Γ/mean(diff(axis))+0.5),1,length(axis))
    fq,fe = factor(Γq,qs),factor(Γe,es)
    return rebin(h;factors=ntuple(d -> d == q ? fq : d == e ? fe : 1,ndims(h.D)),remainder)
end

function _nearest_bin(centres,edges,x)
    first(edges) <= x <= last(edges) || return 0
    j = searchsortedfirst(centres,x)
    j == 1 && return 1
    j > length(centres) && return length(centres)
    return x-centres[j-1] <= centres[j]-x ? j-1 : j
end

"""
    map_nearest_bins(values, source, target; dims=nothing)

Match source and target axes by name, unit and frame, and map a field by nearest bin centre. Omitted singleton field axes broadcast. Use `dims` to select matching axes explicitly. Ties pick the lower bin and out-of-coverage points return `NaN` with a false validity flag.
"""
function map_nearest_bins(values::AbstractArray{<:Real,R},source::NormalizedHistogram{R},target::NormalizedHistogram;
                          dims=nothing) where {R}
    keys = collect(zip(source.names,source.units,source.frames))
    targets = collect(zip(target.names,target.units,target.frames))
    ds = dims === nothing ? Tuple(d for d in 1:R if size(values,d)>1 || (size(source.D,d)==1 && keys[d] in targets)) :
         _reduction_dims(dims,R)
    maps = ntuple(Val(R)) do d
        if d in ds
            size(values,d)==size(source.D,d) || throw(DimensionMismatch("field axis must match source"))
            matches = findall(==(keys[d]),targets)
            length(matches)==1 || throw(ArgumentError("expected one matching target axis for $(source.names[d])"))
            k = only(matches)
            return (k,map(x -> _nearest_bin(source.centres[d],source.edges[d],x),target.centres[k]))
        end
        size(values,d)==1 || throw(DimensionMismatch("unmapped axes must have one bin"))
        return (0,Int[])
    end
    return _map_nearest_field(values,size(target.D),maps)
end

function _map_nearest_field(values::AbstractArray{<:Real,R},shape,maps) where {R}
    mapped,valid = fill(NaN,shape),falses(shape)
    for i in CartesianIndices(mapped)
        j = ntuple(d -> maps[d][1]==0 ? 1 : maps[d][2][i[maps[d][1]]],Val(R))
        all(>(0),j) || continue
        x = values[j...]
        isfinite(x) || continue
        mapped[i],valid[i] = x,true
    end
    return (; values=mapped,valid)
end
