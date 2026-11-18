function _point_operations(group,::Val{P}) where {P}
    ops = SMatrix{P,P,Float64,P*P}[]
    for matrix in (group isa AbstractMatrix ? (group,) : group)
        size(matrix)==(P,P) && all(isfinite,matrix) || throw(ArgumentError("expected finite $P×$P point operations"))
        op = SMatrix{P,P,Float64}(matrix)
        isapprox(abs(det(op)),1.0) || throw(ArgumentError("point operations must preserve reciprocal volume"))
        push!(ops,op)
    end
    isempty(ops) && throw(ArgumentError("supply at least one point operation"))
    return ops
end

function _zone_origin(origin,::Val{P}) where {P}
    origin === nothing && return SVector{P,Float64}(ntuple(_ -> -0.5,Val(P)))
    length(origin)==P && all(isfinite,origin) || throw(ArgumentError("invalid zone origin"))
    return SVector{P,Float64}(origin)
end

function _zone_coordinate(x,lo)
    lo<=x<lo+1 && return x
    y = lo+mod(x-lo,1.0)
    return abs(y-(lo+1)) <= 8eps(max(1.0,abs(y))) ? lo : y
end

_zone_point(q,origin,translations) = translations ? _zone_coordinate.(q,origin) : q

"""
    symmetry_images(hkl, group; translations=true, origin=nothing)

Apply reciprocal-space point group matrices. With translations, wrap into the half-open unit cell starting at `origin` (default -0.5 on each axis). Return distinct images.
"""
function symmetry_images(hkl,group;translations::Bool=true,origin=nothing)
    q = _hkl_point(hkl)
    return _symmetry_images(q,_point_operations(group,Val(length(q))),_zone_origin(origin,Val(length(q))),translations)
end

function _symmetry_images(q::SVector{P},ops,origin,translations) where {P}
    images = NTuple{P,Float64}[]
    for op in ops
        image = Tuple(_zone_point(op*q,origin,translations))
        any(p -> all(isapprox.(image,p;atol=1e-12,rtol=0)),images) || push!(images,image)
    end
    return images
end

_point_at(::Nothing,g,i,dim) = _index_components(_stored_hkl(g,i),dim)
function _point_at(indices::AbstractArray,g,i,dim::Val{P}) where {P}
    return _index_components(indices[i],dim)
end

function _momentum_targets(target::F,g,indices,ops,origin::SVector{P},translations) where {F,P}
    if indices !== nothing
        axes(indices)==map(Base.OneTo,g.shape) || throw(DimensionMismatch("indices must match the momentum grid"))
        all(q -> P<=length(q)<=3 && all(isfinite,q),indices) || throw(ArgumentError("invalid corrected HKL coordinates"))
    end
    starts = Vector{Int}(undef,prod(g.shape)+1)
    destinations, distinct = Int[],Int[]
    for (j,i) in enumerate(CartesianIndices(g.shape))
        starts[j] = length(destinations)+1
        q = _point_at(indices,g,i,Val(P))
        if ops === nothing
            k = target(q)
            k>0 && push!(destinations,k)
        else
            empty!(distinct)
            for op in ops
                k = target(_zone_point(op*q,origin,translations))
                k>0 && !(k in distinct) && push!(distinct,k)
            end
            append!(destinations,distinct)
        end
    end
    starts[end] = length(destinations)+1
    return (;starts,destinations)
end

function _target_entries(targets,momentum_ids,g,i,n)
    k = momentum_ids[_momentum_index(g,i)]
    return targets.starts[k]:targets.starts[k+1]-1,n*(i[g.energy]-1)
end

function _scatter_reciprocal(h,g,targets,shape,conv,background,norms,mask)
    return _scatter_reciprocal(h.D,h.D2,h.N,h.centres[g.energy],g,targets,shape,
                               conv,background,norms,mask,_numerator_covariance(h))
end

function _scatter_reciprocal(D,D2::V,N,energies,g,targets,shape,conv::C,background::B,norms,mask,source_cov) where {V,C,B}
    numerator,weight = zeros(shape),zeros(shape)
    variance = D2 === nothing ? nothing : zeros(shape)
    cov = Dict{Tuple{Int,Int},Float64}()
    active = isempty(source_cov) ? nothing : falses(size(D))
    n = prod(shape) ÷ last(shape)
    momentum_ids = LinearIndices(g.shape)
    dest = targets.destinations
    for i in CartesianIndices(D)
        mask !== nothing && mask[i] && continue
        entries,offset = _target_entries(targets,momentum_ids,g,i,n)
        isempty(entries) && continue
        q = norms === nothing ? 0.0 : norms[_momentum_index(g,i)]
        bin = _conversion_bin(D,N,conv,background,i,q,energies[i[g.energy]])
        bin === nothing && continue
        d,w,var = bin
        active === nothing || (active[i] = true)
        D2 === nothing || (var += D2[i])
        for j in entries
            k = dest[j]+offset
            numerator[k] += d
            weight[k] += w
            if variance !== nothing
                variance[k] += var
                for l in first(entries):j-1
                    _record_covariance!(cov,variance,k,dest[l]+offset,var)
                end
            end
        end
    end
    if active !== nothing
        source_ids = CartesianIndices(D)
        for (a,b,v) in source_cov
            active[a] && active[b] || continue
            as,oa = _target_entries(targets,momentum_ids,g,source_ids[a],n)
            bs,ob = _target_entries(targets,momentum_ids,g,source_ids[b],n)
            for j in as,k in bs
                _record_covariance!(cov,variance,dest[j]+oa,dest[k]+ob,v)
            end
        end
    end
    return numerator,variance,weight,_covariance_entries(cov)
end

function _zone_bin(q,edges,linear,translations)
    ids = map(edges,Tuple(q)) do edge,x
        k = searchsortedlast(edge,x)
        if k < length(edge) && edge[k+1]-x <= 8eps(max(1.0,abs(x)))
            k += 1
        end
        k == length(edge) && !translations && x<=last(edge) ? k-1 : k
    end
    all(d -> 1<=ids[d]<length(edges[d]),eachindex(ids)) || return 0
    return linear[ids...]
end

"""
    fold_to_zone(h, group; grid, indices=nothing, translations=true)

Accumulate distinct images into a grid described by two or three HKL edge vectors. Energy becomes the last output axis. Each source contributes once per destination. Shared source contributions are recorded as off-diagonal linear-index numerator covariance in `meta.numerator_covariance`. Use `indices` for corrected HKL.
"""
function fold_to_zone(h::NormalizedHistogram,group;grid,indices=nothing,translations::Bool=true)
    edges = map(x -> collect(Float64,x),Tuple(grid))
    length(edges) in (2,3) || throw(ArgumentError("grid must specify two or three HKL axes"))
    for edge in edges
        length(edge)>=2 && all(isfinite,edge) && all(i -> edge[i]>edge[i-1],2:length(edge)) ||
            throw(ArgumentError("zone edges must be finite and strictly increasing"))
        translations && last(edge)-first(edge)>1+8eps(Float64) && throw(ArgumentError("zone spans more than one period"))
    end
    return _fold_to_zone(h,group,edges,indices,translations)
end

function _fold_to_zone(h,group,edges::NTuple{P},indices,translations) where {P}
    shape = map(edge -> length(edge)-1,edges)
    g = _momentum_geometry(h)
    ops = _point_operations(group,Val(P))
    linear = LinearIndices(shape)
    target = q -> _zone_bin(q,edges,linear,translations)
    targets = _momentum_targets(target,g,indices,ops,SVector{P,Float64}(map(first,edges)),translations)
    e = g.energy
    D,D2,N,cov = _scatter_reciprocal(h,g,targets,(shape...,size(h.D,e)),nothing,nothing,nothing,nothing)
    centres = map(edge -> [(edge[i]+edge[i+1])/2 for i in 1:length(edge)-1],edges)
    names = [collect(("H","K","L")[1:P]);h.names[e]]
    units,frames = [fill("r.l.u.",P);h.units[e]],[fill("HKL",P);h.frames[e]]
    return _reduced_histogram(h,D,D2,N,[collect(edges);[h.edges[e]]],[collect(centres);[h.centres[e]]],cov;
                              names,units,frames)
end

function _path_geometry(path,width,step)
    vertices = map(_hkl_point,collect(path))
    length(vertices)>=2 || throw(ArgumentError("path needs at least two vertices"))
    p = length(first(vertices))
    all(q -> length(q)==p,vertices) || throw(DimensionMismatch("inconsistent path dimensions"))
    points = SVector{p,Float64}.(vertices)
    segments = diff(points)
    lengths = norm.(segments)
    all(>(0),lengths) || throw(ArgumentError("consecutive path vertices must differ"))
    isfinite(width) && width>=0 && isfinite(step) && step>0 || throw(ArgumentError("invalid path width or step"))
    cumulative = [0.0;cumsum(lengths)]
    edges = collect(0.0:step:last(cumulative))
    if length(edges)>1 && isapprox(last(edges),last(cumulative);atol=0,rtol=8eps(Float64))
        edges[end] = last(cumulative)
    else
        push!(edges,last(cumulative))
    end
    return (;points,directions=segments ./ lengths,lengths,cumulative,edges,width2=Float64(width)^2)
end

function _path_bin(q,path)
    distance = Inf
    s = 0.0
    for k in eachindex(path.lengths)
        Δ = q-path.points[k]
        t = clamp(dot(Δ,path.directions[k]),0,path.lengths[k])
        d = sum(abs2,Δ-t*path.directions[k])
        if d<distance
            distance,s = d,path.cumulative[k]+t
        end
    end
    distance <= path.width2+(16eps(max(1.0,maximum(abs,q))))^2 || return 0
    return min(searchsortedlast(path.edges,s),length(path.edges)-1)
end

"""
    cut_path(h, path; width, step, fold=nothing, background=nothing, conv=nothing, exclude=nothing, indices=nothing,
             translations=true, origin=nothing, reciprocal_basis=nothing)

Cut along a 2D or 3D HKL polyline using a Euclidean half-width and arc-length step in reciprocal lattice coordinates. Use a Boolean source mask or a predicate of corrected image HKL for exclusions. Apply background and conversion on the measured grid. Return the ratio `Σ(D-N*b)/Σ(N*g)`, its standard error, coverage, effective denominator, numerator variance, and normalized off-diagonal covariance.
"""
function cut_path(h::NormalizedHistogram,path;width::Real,step::Real,fold=nothing,background=nothing,conv=nothing,
                  exclude=nothing,indices=nothing,translations::Bool=true,origin=nothing,reciprocal_basis=nothing)
    return _cut_path(h,_path_geometry(path,width,step),fold,background,conv,exclude,indices,translations,origin,reciprocal_basis)
end

function _cut_path(h,path,fold,background,conv,exclude,indices,translations,origin,basis)
    g = _momentum_geometry(h)
    _check_bin_array(conv,h,"conv")
    _check_bin_array(background,h,"background")
    exclude isa AbstractArray && !(exclude isa AbstractArray{Bool}) && throw(ArgumentError("exclude mask must be Boolean"))
    mask = exclude isa AbstractArray{Bool} ? exclude : nothing
    _check_bin_array(mask,h,"exclude")
    predicate = mask === nothing ? exclude : nothing
    norms = _background_norms(h,g,background,basis)
    dim = Val(length(first(path.points)))
    ops = fold === nothing ? nothing : _point_operations(fold,dim)
    target = q -> predicate !== nothing && predicate(q) ? 0 : _path_bin(q,path)
    targets = _momentum_targets(target,g,indices,ops,_zone_origin(origin,dim),translations)
    shape = (length(path.edges)-1,size(h.D,g.energy))
    numerator,variance,weight,cov = _scatter_reciprocal(h,g,targets,shape,conv,background,norms,mask)
    valid = weight .> 0
    values = map((d,n) -> n>0 ? d/n : NaN,numerator,weight)
    σ = variance === nothing ? nothing : map((v,n) -> n>0 ? sqrt(v)/n : NaN,variance,weight)
    covariance = _CovarianceEntry[(i,j,v/(weight[i]*weight[j])) for (i,j,v) in cov if valid[i] && valid[j]]
    s = [(path.edges[i]+path.edges[i+1])/2 for i in 1:shape[1]]
    return (;values,σ,weight,s,ω=h.centres[g.energy],numerator,variance,valid,covariance)
end
