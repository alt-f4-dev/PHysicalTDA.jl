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
	    Z = superlevel ? -A : A
    end
	
    PD = ripserer(Cubical(Z); dim_max=maxdim)
    return PD
end
"""
Wrapper: Sunny.Intensities -> dense Array -> pd_array_intensity

Converts `I.data` to a dense array and calls `pd_array_intensity`.
Respects `maxdim` and `superlevel`. Returns PD.

Convenient for topology on LSWT or phonon intensity objects.
"""
function pd_sunny_intensities(I::Sunny.Intensities{T,G,N};
			    maxdim::Int=1, superlevel::Bool=true, normalize::Bool=true) where {T<:Real,G,N}
    A = Array(I.data) # whatever rank Sunny gives here
    return pd_array_intensity(A; maxdim=maxdim, superlevel=superlevel, normalize=normalize)
end
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
    S = Dict{Int,Float64}(); E = Dict{Int,Float64}()
    for p in dims
	bd = p+1 ≤ length(PDs) ? PDs[p+1] : ()
	isempty(bd) && (S[p] = 0.0; E[p] = 0.0; continue)
	# lifetimes
	ts = Float64[death(x) - birth(x) for x in bd]
	ts = filter(t -> isfinite(t) && t > tol, ts)
	isempty(ts) && (S[p] = 0.0; E[p] = 0.0; continue)
	# total persistence & fractional lifetimes 
        Eₚ = sum(ts); Tₚ  = ts ./ Eₚ
	# persistence entropy; add tiny epsilon to avoid log(0) if needed
	Sₚ  = -sum(@view(Tₚ[:]) .* log.(@view(Tₚ[:]) .+ eps()))  # natural log
	S[p] = Sₚ; E[p] = Eₚ
    end
    return S, E
end
"""
    persistence_entropy_curve(PDs, τs; dims = 0:1, tol = 0.0)

Filtration-dependent persistence entropy Sₚ(τⱼ) and total truncated persistence Eₚ(τⱼ).

- `PDs`  :: Vector of Ripserer persistence diagrams grouped by homology degree.
- `τs`   :: sorted vector of filtration parameters (monotone increasing).
- `dims` :: homology degrees to use (default 0:1).
- `tol`  :: discard contributions with truncated lifetime ≤ tol.

For each p and τⱼ, lifetimes are truncated as
    ℓᵢ(τⱼ) = max(0, min(dᵢ, τⱼ) - bᵢ)
and only ℓᵢ(τⱼ) > tol and finite dᵢ are used.

Returns `(S::Dict{Int,Vector{Float64}}, E::Dict{Int,Vector{Float64}})` keyed by p,
where S[p][j] = Sₚ(τⱼ), E[p][j] = Eₚ(τⱼ).
"""
function persistence_entropy_curve(PDs::Vector{PersistenceDiagrams.PersistenceDiagram}, τs::AbstractVector{<:Real}; dims = 0:1, tol::Real = 0.0)
    @assert issorted(τs) "τs must be sorted ascending"

    S = Dict{Int,Vector{Float64}}()
    E = Dict{Int,Vector{Float64}}()

    for p in dims
        bd = p + 1 <= length(PDs) ? PDs[p + 1] : ()
        if isempty(bd)
            S[p] = zeros(Float64, length(τs))
            E[p] = zeros(Float64, length(τs))
            continue
        end

        births = Float64[birth(x) for x in bd]
        deaths = Float64[death(x) for x in bd]

        Sₚ = zeros(Float64, length(τs))
        Eₚ = zeros(Float64, length(τs))

        @inbounds for (j, τ) in pairs(τs)
            ts = Float64[]
            for (b, d) in zip(births, deaths)
                isfinite(d) || continue
                ℓ = min(d, τ) - b
                ℓ > tol && push!(ts, ℓ)
            end

            if isempty(ts)
                Sₚ[j] = 0.0
                Eₚ[j] = 0.0
            else
                Eτ = sum(ts)
                Tτ  = ts ./ Eτ
                Sτ = -sum(@view(Tτ[:]) .* log.(@view(Tτ[:])))
                Sₚ[j] = Sτ
                Eₚ[j] = Eτ
            end
        end

        S[p] = Sₚ
        E[p] = Eₚ
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
        βₚ = Vector{Int}(undef, length(τ))
        @inbounds for j in eachindex(τ)
            t  = τ[j]
            nb = searchsortedlast(b, t)             # # of births ≤ t
            nd = searchsortedfirst(d, t) - 1        # # of deaths  < t
            βₚ[j] = nb - nd
        end
        β[p] = βₚ
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
    β = betti_curve(pd, τ; dims); κ = Dict{Int, Vector{Float64}}(); Δτ = diff(collect(τ))
    @assert all(Δτ .> 0) "τ grid must be strictly increasing"
    for p in dims
        βₚ = β[p]; κₚ = zeros(Float64, length(βₚ))
        if scheme === :forward
            @inbounds for j in 1:length(Δτ)
                κₚ[j] = (βₚ[j+1] - βₚ[j]) / Δτ[j]
            end
            κₚ[end] = κₚ[end-1]
        elseif scheme === :central
            κₚ[1] = (βₚ[2] - βₚ[1]) / Δτ[1]
            @inbounds for j in 2:length(βₚ)-1
                dτ = τ[j+1] - τ[j-1]
                κₚ[j] = (βₚ[j+1] - βₚ[j-1]) / dτ
            end
            κₚ[end] = (βₚ[end] - βₚ[end-1]) / Δτ[end-1]
        else
            error("Unknown scheme=$scheme (use :forward or :central)")
        end
        κ[p] = κₚ
    end
    return β, κ
end




