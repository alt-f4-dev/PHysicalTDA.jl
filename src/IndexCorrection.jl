function _momentum_geometry(h::NormalizedHistogram{R}) where {R}
    energy = energy_axis(h)
    dims = ntuple(i -> i < energy ? i : i+1,Val(R-1))
    directions = map(d -> SVector{3,Float64}(axis_vector(h.names[d])),dims)
    centres = map(d -> h.centres[d],dims)
    return (;energy,dims,directions,centres,shape=map(length,centres))
end

function _stored_hkl(g,i::CartesianIndex{R}) where {R}
    q = zero(SVector{3,Float64})
    for d in 1:R
        q += g.centres[d][i[d]] * g.directions[d]
    end
    return q
end

_momentum_index(g,i) = CartesianIndex(map(d -> i[d],g.dims))
_index_components(q,::Val{P}) where {P} = SVector{P,Float64}(ntuple(i -> q[i],Val(P)))

function _hkl_point(p)
    length(p) in (2,3) && all(isfinite,p) || throw(ArgumentError("expected finite 2D or 3D HKL coordinates"))
    return SVector{length(p),Float64}(p)
end

"""
    find_bragg_peaks(h; positions, window=0.15, energy_window, min_contrast=5)

Search elastic intensity near expected HKL positions. Return nominal positions, observed centroids weighted by excess above the local median, and peak contrasts relative to local robust noise and statistical uncertainty.
"""
function find_bragg_peaks(h::NormalizedHistogram;positions,window=0.15,energy_window,min_contrast::Real=5)
    targets = map(_hkl_point,collect(positions))
    isempty(targets) && throw(ArgumentError("supply Bragg positions"))
    p = length(first(targets))
    all(q -> length(q)==p,targets) || throw(DimensionMismatch("Bragg positions must share a dimension"))
    widths = SVector{p,Float64}(window isa Real ? ntuple(_ -> window,p) : window)
    all(x -> isfinite(x) && x>0,widths) || throw(ArgumentError("invalid search width"))
    isfinite(min_contrast) && min_contrast>=0 || throw(ArgumentError("invalid minimum contrast"))
    g = _momentum_geometry(h)
    return _bragg_observations(integrate(h;dims=g.energy,range=energy_window),g,targets,widths,min_contrast,Val(p))
end

function _bragg_observations(h,g,targets,widths,min_contrast,::Val{P}) where {P}
    points = map(i -> _index_components(_stored_hkl(g,i),Val(P)),CartesianIndices(g.shape))
    out = NamedTuple{(:P,:observed,:contrast),Tuple{NTuple{P,Float64},NTuple{P,Float64},Float64}}[]
    for target in targets
        ids,values = Int[],Float64[]
        for i in eachindex(points)
            h.N[i]>0 && all(abs.(points[i]-target) .<= widths) || continue
            x = h.D[i]/h.N[i]
            isfinite(x) || continue
            push!(ids,i)
            push!(values,x)
        end
        length(values)>=3 || continue
        baseline = median(values)
        height,j = findmax(values)
        height>baseline || continue
        noise = 1.4826*median(abs.(values .- baseline))
        if h.D2 !== nothing
            i = ids[j]
            noise = max(noise,sqrt(h.D2[i])/h.N[i])
        end
        contrast = (height-baseline)/max(noise,eps(max(1.0,abs(height))))
        contrast>=min_contrast || continue
        centre = zero(target)
        total = 0.0
        for (i,x) in zip(ids,values)
            w = max(0.0,x-baseline)
            centre += w*points[i]
            total += w
        end
        push!(out,(;P=Tuple(target),observed=Tuple(centre/total),contrast))
    end
    return out
end

"""
    fit_index_correction(peaks; weights=:contrast)

Fit `observed = M * P + t` using contrast, uniform, or explicit nonnegative weights. Return `M`, `t`, and the weighted RMS residuals before and after correction.
"""
function fit_index_correction(peaks;weights=:contrast)
    records = collect(peaks)
    isempty(records) && throw(ArgumentError("no Bragg observations"))
    p = length(_hkl_point(first(records).P))
    w = if weights === :contrast
        Float64[r.contrast for r in records]
    elseif weights === :uniform
        ones(length(records))
    elseif weights isa AbstractVector
        collect(Float64,weights)
    else
        throw(ArgumentError("invalid fit weights"))
    end
    length(w)==length(records) || throw(DimensionMismatch("one weight per peak required"))
    all(x -> isfinite(x) && x>=0,w) && count(>(0),w)>=p+1 ||
        throw(ArgumentError("affine fit needs at least $(p+1) finite positive weights"))
    X,Y = ones(length(records),p+1),zeros(length(records),p)
    for (i,r) in enumerate(records)
        x,y = _hkl_point(r.P),_hkl_point(r.observed)
        length(x)==length(y)==p || throw(DimensionMismatch("inconsistent Bragg dimensions"))
        X[i,1:p] .= x
        Y[i,:] .= y
    end
    root_w = sqrt.(w ./ maximum(w))
    WX = X .* root_w
    rank(WX)==p+1 || throw(ArgumentError("Bragg positions do not span an affine basis"))
    coeff = qr(WX) \ (Y .* root_w)
    M = permutedims(coeff[1:p,:])
    rank(M)==p || throw(ArgumentError("fitted affine map is singular"))
    rms(A) = norm(A .* root_w)/norm(root_w)
    return (;M,t=coeff[end,:],rms=rms(X*coeff-Y),rms_uncorrected=rms(X[:,1:p]-Y))
end

"""
    corrected_indices(h, M, t)

Apply the inverse affine correction to stored HKL. Return a momentum-shaped array of HKL tuples. A 2D fit corrects H and K and preserves L.
"""
function corrected_indices(h::NormalizedHistogram,M::AbstractMatrix,t)
    p = length(t)
    p in (2,3) && size(M)==(p,p) || throw(DimensionMismatch("expected a 2D or 3D affine map"))
    all(isfinite,M) && all(isfinite,t) && rank(M)==p || throw(ArgumentError("invalid affine map"))
    return _corrected_indices(_momentum_geometry(h),SMatrix{p,p,Float64}(inv(Matrix{Float64}(M))),SVector{p,Float64}(t))
end

function _corrected_indices(g,A::SMatrix{P,P,Float64},t) where {P}
    return map(CartesianIndices(g.shape)) do i
        q = _stored_hkl(g,i)
        x = A * (_index_components(q,Val(P))-t)
        P == 2 ? (x[1],x[2],q[3]) : Tuple(x)
    end
end
