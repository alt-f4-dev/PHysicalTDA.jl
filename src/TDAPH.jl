"""
Median Absolute Deviation (MAD)
Applies elementwise relative to the median of `t`. 
Suitable for thresholding outliers in intensity arrays.
Returns an array of MAD values with the same shape as `t`.
"""
@inline MAD(t) = 1.48 .* median.(abs.(t .- median(vec(t))))
"""
    persistence_entropy(PDs; dims = 0:1, tol = 0.0)

Compute per-dimension persistence entropy Sₚ and total persistence Eₚ.

- `PDs` is a vector of Ripserer persistence diagrams, grouped by homology degree as returned by `pd_array_intensities()` or `pd_sunny_intensities()`.
- `dims` selects homology degrees (default 0:1).
- `tol` discards `lifetimes ≤ tol` (guards numerical noise).

Returns `(S::Dict{Int,Float64}, E::Dict{Int,Float64})` keyed by p.
"""
function persistence_entropy(PDs::Vector{PersistenceDiagrams.PersistenceDiagram}; dims = 0:1, tol::Real = 0.0)
	S = Dict{Int,Float64}()
	E = Dict{Int,Float64}()
	for p in dims
		bd = p+1 <= length(PDs) ? PDs[p+1] : ()
		isempty(bd) && (S[p] = 0.0; E[p] = 0.0; continue)
		# lifetimes
		τ = Float64[death(x) - birth(x) for x in bd]
		τ = filter(t -> isfinite(t) && t > tol, τ)
		isempty(τ) && (S[p] = 0.0; E[p] = 0.0; continue)
		Ep = sum(τ)
		q  = τ ./ Ep
		# entropy; add tiny epsilon to avoid log(0) if needed
		H  = -sum(@view(q[:]) .* log.(@view(q[:]) .+ eps()))  # natural log
		S[p] = H
		E[p] = Ep
	end
	return S, E
end
"""
    betti_curve(PDs, τ; dims = 0:1)

Compute per-dimension Betti curves β_p(τ_j) on a user-provided grid `τ`.

- `PDs` is a vector of Ripserer PD (grouped by degree).
- `τ` is a sorted vector of thresholds (monotone increasing).

Returns `Dict{Int,Vector{Int}}` mapping p ↦ β_p(τ).
"""
function betti_curve(PDs::Vector{PersistenceDiagrams.PersistenceDiagram}, τ::AbstractVector{<:Real}; dims = 0:1)
    @assert issorted(τ) "τ must be sorted ascending"
    β = Dict{Int, Vector{Int}}()
    for p in dims
        bd = p+1 <= length(PDs) ? PDs[p+1] : ()
        if isempty(bd)
            β[p] = zeros(Int, length(τ))
            continue
        end
        b = sort!(Float64[birth(x) for x in bd])
        d = sort!(Float64[death(x) for x in bd])
        # β(τ) = count(b ≤ τ) - count(d < τ)  (half-open [b,d))
        βp = Vector{Int}(undef, length(τ))
        @inbounds for j in eachindex(τ)
            t  = τ[j]
            nb = searchsortedlast(b, t)             # # of births ≤ t
            nd = searchsortedfirst(d, t) - 1        # # of deaths  < t
            βp[j] = nb - nd
        end
        β[p] = βp
    end
    return β
end

"""
    betti_curvature(pd, τ; dims = 0:1, scheme = :forward)

Compute β_p(τ) and κ_p(τ) = dβ_p/dτ on a grid `τ`.

- `scheme = :forward` (default) uses forward differences and sets κ[end]=κ[end-1].
  Use `:central` for central differences on interior points.

Returns `(β::Dict{Int,Vector{Int}}, κ::Dict{Int,Vector{Float64}})`.
"""
function betti_curvature(pd, τ::AbstractVector{<:Real}; dims = 0:1, scheme::Symbol = :forward)
    @assert length(τ) ≥ 2 "τ grid must have at least 2 points"
    β = betti_curve(pd, τ; dims)
    κ = Dict{Int, Vector{Float64}}()
    Δτ = diff(collect(τ))
    @assert all(Δτ .> 0) "τ grid must be strictly increasing"
    for p in dims
        βp = β[p]
        κp = zeros(Float64, length(βp))
        if scheme === :forward
            @inbounds for j in 1:length(Δτ)
                κp[j] = (βp[j+1] - βp[j]) / Δτ[j]
            end
            κp[end] = κp[end-1]
        elseif scheme === :central
            κp[1] = (βp[2] - βp[1]) / Δτ[1]
            @inbounds for j in 2:length(βp)-1
                dt = τ[j+1] - τ[j-1]
                κp[j] = (βp[j+1] - βp[j-1]) / dt
            end
            κp[end] = (βp[end] - βp[end-1]) / Δτ[end-1]
        else
            error("Unknown scheme=$scheme (use :forward or :central)")
        end
        κ[p] = κp
    end
    return β, κ
end
"""
Computes persistance diagrams of an intensity array via cubical complexes.
Returns PDs with birth–death scatter per dimension, grouped by homology
degree up to `maxdim`.

Optional superlevel filtration (invert intensities) and normalization to [0,1].
Used to analyze the topology of I(Q,ω) = |S(Q,ω)|².
"""
function pd_array_intensity(A::AbstractArray{<:Real,N};
			    maxdim::Int=1, superlevel::Bool=true, 
			    normalize::Bool=false) where {N}
    if any(x -> !isfinite(x), A)
        A = map(x -> isfinite(x) ? x : 0.0, A)
    end

    if normalize
            amin, amax = extrema(A)
            if amax == amin
                Z = zeros(size(A))
            else
                A = (A .- amin) ./ (amax - amin)
                Z = superlevel ? 1 .- A : A
            end
    else
	    Z = superlevel ? 1 .- A : A
    end
	
    PD = ripserer(Cubical(Z); dim_max=maxdim)
    for d in 0:maxdim
	    bd = PD[d+1]; isempty(bd) && continue
	    bth = birth.(bd); dth = death.(bd)
    end
    return PD
end
"""
Wrapper: Sunny.Intensities -> dense Array -> pd_array_intensity

Converts `I.data` to a dense array and calls `pd_array_intensity`.
Respects `maxdim` and `superlevel`. Returns (PD, Figures).

Convenient for topology on LSWT or phonon intensity objects.
"""
function pd_sunny_intensities(I::Sunny.Intensities{T,G,N};
			    maxdim::Int=1, superlevel::Bool=true, normalize::Bool=true) where {T<:Real,G,N}
	A = Array(I.data) # whatever rank Sunny gives here
	return pd_array_intensity(A; maxdim=maxdim, superlevel=superlevel, normalize=normalize)
end
