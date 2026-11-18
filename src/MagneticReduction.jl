function _moment_direction(axis)
    m = SVector{3,Float64}(axis)
    n = norm(m)
    isfinite(n) && n>0 || throw(ArgumentError("moment_axis must be finite and nonzero"))
    return m/n
end

function _polarization_channel(channel)
    channel in (:transverse,:longitudinal,:isotropic) || throw(ArgumentError("invalid polarization channel"))
    return Val(channel)
end

_polarization(x,::Val{:transverse}) = (1+x^2)/2
_polarization(x,::Val{:longitudinal}) = 1-x^2
_polarization_vector(Q,n,m,::Val{:isotropic}) = 2/3
_polarization_vector(Q,n,m,mode) = iszero(n) ? NaN : _polarization(clamp(dot(Q,m)/n,-1,1),mode)

"""
    polarization_factor(Qhat, moment_axis; channel=:transverse)

Return the projection relative to a Cartesian moment axis for the transverse, longitudinal, or isotropic channel. The transverse factor averages two components. Zero momentum has undefined anisotropic projection.
"""
function polarization_factor(Qhat,moment_axis;channel::Symbol=:transverse)
    Q = SVector{3,Float64}(Qhat)
    all(isfinite,Q) || throw(ArgumentError("Q must be finite"))
    return _polarization_vector(Q,norm(Q),_moment_direction(moment_axis),_polarization_channel(channel))
end

function _reciprocal_basis(h,basis)
    if basis !== nothing
        size(basis)==(3,3) && all(isfinite,basis) && rank(basis)==3 || throw(ArgumentError("invalid reciprocal basis"))
        return SMatrix{3,3,Float64}(basis)
    end
    cell = lattice(h)
    cell === nothing && throw(ArgumentError("histogram has no lattice; supply reciprocal_basis"))
    all(x -> x isa Real && isfinite(x) && x>0,values(cell)) || throw(ArgumentError("invalid lattice parameters"))
    a,b,c,α,β,γ = values(cell)
    all(<(180),(α,β,γ)) || throw(ArgumentError("invalid lattice angles"))
    x = cosd(β)
    y = (cosd(α)-x*cosd(γ))/sind(γ)
    z2 = 1-x^2-y^2
    z2>0 || throw(ArgumentError("cell volume must be positive"))
    direct = @SMatrix [a b*cosd(γ) c*x; 0.0 b*sind(γ) c*y; 0.0 0.0 c*sqrt(z2)]
    return 2π*transpose(inv(direct))
end

"""
    magnetic_conversion(h; ion=:Mn2, moment_axis=(0,0,1), channel=:transverse, average_over=nothing, reciprocal_basis=nothing)

Evaluate squared form factor times polarization at measured HKL. The reciprocal basis includes `2π` in Å⁻¹. The default basis places direct `a` along Cartesian x and `a×b` along z. Provide Mn²⁺, form factor coefficients, or a function of `|Q|`. Use `average_over` for exposure-weighted averaging.
"""
function magnetic_conversion(h::NormalizedHistogram;ion=:Mn2,moment_axis=(0,0,1),channel::Symbol=:transverse,
                             average_over=nothing,reciprocal_basis=nothing)
    g = _momentum_geometry(h)
    basis = _reciprocal_basis(h,reciprocal_basis)
    F = _form_factor_function(ion === :Mn2 ? MN2_J0 : ion)
    factors = _momentum_factors(g,basis,F,_moment_direction(moment_axis),_polarization_channel(channel))
    field = i -> factors[_momentum_index(g,i)]
    return average_over === nothing ? map(field,CartesianIndices(h.D)) : weighted_average(h,field;dims=average_over)
end

function _momentum_factors(g,basis,F::T,m,mode) where {T}
    return map(CartesianIndices(g.shape)) do i
        Q = basis*_stored_hkl(g,i)
        n = norm(Q)
        abs2(F(n))*_polarization_vector(Q,n,m,mode)
    end
end

_background_value(::Nothing,i,q,ω) = 0.0
_background_value(b::Real,i,q,ω) = b
_background_value(b::AbstractArray,i,q,ω) = b[i]
_background_value(b,i,q,ω) = b(q,ω)
_background_parts(b::Real) = (Float64(b),0.0)
_background_parts(b::Tuple{Real,Real}) = (Float64(b[1]),Float64(b[2]))
_background_parts(b::NamedTuple{(:value,:σ)}) = _background_parts(Tuple(b))

_uses_background_geometry(b) = !(b isa Union{Nothing,Real,AbstractArray})

function _background_norms(h,g,b,basis)
    _uses_background_geometry(b) || return nothing
    B = _reciprocal_basis(h,basis)
    return map(i -> norm(B*_stored_hkl(g,i)),CartesianIndices(g.shape))
end

function _conversion_bin(D,N,conv,background,i,q,ω)
    n = N[i]
    n>0 && isfinite(n) && isfinite(D[i]) || return nothing
    g = _bin_value(conv,i)
    isfinite(g) || return nothing
    g>=0 || throw(ArgumentError("conversion must be nonnegative"))
    b,σ = _background_parts(_background_value(background,i,q,ω))
    isfinite(b) && isfinite(σ) || return nothing
    σ>=0 || throw(ArgumentError("background uncertainty must be nonnegative"))
    return D[i]-n*b,n*g,abs2(n*σ)
end

"""
    to_smag(h::NormalizedHistogram, g; background=nothing, reciprocal_basis=nothing)

Convert to `(D-N*b)/(N*g)`. Background may be a scalar, array, or callable of measured `(|Q|, ω)` returning a value or `(value, σ)`. Background errors are treated as independent across source bins. Return `nothing` for total uncertainty when source variance is unavailable.
"""
function to_smag(h::NormalizedHistogram,g;background=nothing,reciprocal_basis=nothing)
    _check_bin_array(g,h,"g")
    _check_bin_array(background,h,"background")
    geom = _uses_background_geometry(background) ? _momentum_geometry(h) : nothing
    norms = _background_norms(h,geom,background,reciprocal_basis)
    coords = geom === nothing ? i -> (0.0,0.0) :
        i -> (norms[_momentum_index(geom,i)],h.centres[geom.energy][i[geom.energy]])
    return _to_smag(h.D,h.D2,h.N,coords,g,background)
end

function _to_smag(D,D2::V,N,coords::F,g,background) where {V,F}
    S = fill(NaN,size(D))
    σ = D2 === nothing ? nothing : fill(NaN,size(D))
    for i in CartesianIndices(D)
        q,ω = coords(i)
        bin = _conversion_bin(D,N,g,background,i,q,ω)
        bin === nothing && continue
        numerator,weight,var = bin
        weight>0 || continue
        S[i] = numerator/weight
        σ === nothing || (σ[i] = sqrt(D2[i]+var)/weight)
    end
    return S,σ
end
