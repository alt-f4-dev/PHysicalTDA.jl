"""
    pd_array_intensities(A; maxdim=1, threshold=nothing, superlevel=true, normalize=false,
                         make_death_finite=true, essential_cap=:max_death)

Compute persistence diagrams of an intensity array using a cubical complex.

Input:

- `A`             :: N-dimensional intensity array (non-finite entries are replaced by 0)
- `maxdim`        :: maximum homology degree `p`
- `threshold`     :: fractional filtration cutoff (see below); `nothing` uses the full range
- `superlevel`    :: superlevel (`true`) or sublevel (`false`) filtration
- `normalize`     :: min–max normalize `A` to [0, 1] before filtering
- `make_death_finite` :: replace infinite deaths by a finite cap
- `essential_cap` :: cap used for infinite deaths when `make_death_finite=true`

Returns:
- `PDs` :: persistence diagrams `{Dₚ}` for `p ∈ [0, maxdim]`

# Filtration values

The filtration is computed on `Z`, not on `A`:

| `superlevel` | `normalize` | `Z`              |
|:-------------|:------------|:-----------------|
| `true`       | `false`     | `-A`             |
| `true`       | `true`      | `1 - Â`          |
| `false`      | `false`     | `A`              |
| `false`      | `true`      | `Â`              |

where `Â = (A - min A)/(max A - min A)`. Births and deaths are reported in `Z`
units, so for superlevel filtrations they are negated (or complemented)
intensities and lifetimes are `death - birth ≥ 0`.

# Threshold

`threshold` is a fraction (not a percent):

- `normalize=false`: superlevel keeps cells with `A ≥ threshold * maximum(A)`;
  sublevel keeps cells with `A ≤ threshold * maximum(A)` (for `A ≥ 0`).
- `normalize=true`: the fraction refers to the normalized range, i.e. superlevel
  keeps `Â ≥ threshold` and sublevel keeps `Â ≤ threshold`.

# Infinite deaths

With `make_death_finite=true`, each infinite death is replaced by a cap:

- `essential_cap=:max_death` (default, historical behaviour): the largest finite
  death in that degree, or `maximum(births) + max(1, 0.05max(|bmax|, 1))` when
  no finite death exists. The essential class then depends on the other
  features of the diagram.
- `essential_cap=:threshold`: the filtration value at which the complex is
  complete (the full-range maximum of `Z`, or the threshold). For a
  superlevel filtration this is the background (minimum) intensity level, so
  the essential `H₀` lifetime equals the peak-to-floor contrast.

The returned diagrams are rebuilt from (birth, death) pairs and do not carry
Ripserer's representative/simplex metadata.
"""
function pd_array_intensities(A::AbstractArray{<:Real,N};
                              maxdim::Int=1,
                              threshold::Union{Real,Nothing}=nothing,
                              superlevel::Bool=true,
                              normalize::Bool=false,
                              make_death_finite::Bool=true,
                              essential_cap::Symbol=:max_death) where {N}
    essential_cap in (:max_death, :threshold) ||
        throw(ArgumentError("essential_cap must be :max_death or :threshold; got $essential_cap."))
    if any(x -> !isfinite(x), A) 
        A = map(x -> isfinite(x) ? x : 0.0, A)
    end
    θ = threshold; amin, amax = extrema(A)
    if normalize
        if amax == amin
            Z = zeros(size(A))
            τ = 0.0
        else
            A = (A .- amin) ./ (amax - amin)
            Z = superlevel ? 1 .- A : A
            if !isnothing(θ)
                τ = superlevel ? 1 - θ : θ
            else
                #τ = superlevel ? minimum(Z) : maximum(Z)
                τ = maximum(Z)
            end
        end
    else
	Z = superlevel ? -A : A
        if !isnothing(θ)
            Z⁻, Z⁺ = extrema(Z)
            τ = superlevel ? θ * Z⁻ : θ * Z⁺
        else
            #τ = superlevel ? minimum(Z) : maximum(Z)
            τ = maximum(Z)
        end
    end
    PD = ripserer(Cubical(Z, threshold=τ); dim_max=maxdim)

    if make_death_finite
        PD = map(PD) do Δ
            isempty(Δ) && return Δ

            births = Float64[birth(x) for x in Δ if isfinite(birth(x))]
            deaths = Float64[death(x) for x in Δ if isfinite(death(x))]

            bmax = isempty(births) ? 0.0 : maximum(births)

            cap = if essential_cap === :threshold
                # Filtration value at which the complex is complete.
                Float64(τ)
            elseif isempty(deaths)
                # All deaths are Inf: choose a finite cap above births.
                bmax + max(1.0, 0.05 * max(abs(bmax), 1.0))
            else
                dmax = maximum(deaths)
                dmax
            end

            pts = map(Δ) do x
                b = birth(x)
                d = death(x)
                d2 = isfinite(d) ? d : cap
                (b, max(d2, b))  # keep valid ordering
            end

            PersistenceDiagrams.PersistenceDiagram(pts)
        end
    end

    return PD
end
"""
    pd_sunny_intensities(I; maxdim=1, threshold=nothing, superlevel=true, normalize=true,
                         make_death_finite=true, essential_cap=:max_death)

Wrapper: `Sunny.Intensities` → dense `Array` → [`pd_array_intensities`](@ref).

Note that `normalize` defaults to `true` here but to `false` in
`pd_array_intensities`. All keywords are forwarded. Returns the vector of PDs.
"""
function pd_sunny_intensities(I::Sunny.Intensities{T,G,N};
			      maxdim::Int=1,
                              threshold::Union{Real, Nothing}=nothing,
                              superlevel::Bool=true,
                              normalize::Bool=true,
                              make_death_finite::Bool=true,
                              essential_cap::Symbol=:max_death) where {T<:Real,G,N}
    A = Array(I.data) # whatever rank Sunny gives here
    return pd_array_intensities(A; maxdim=maxdim,
                                 superlevel=superlevel,
                                 threshold=threshold,
                                 normalize=normalize,
                                 make_death_finite=make_death_finite,
                                 essential_cap=essential_cap)
end
"""
    MAD(t; elementwise=false)

Scaled median absolute deviation, `1.4826 * median(abs.(t .- median(t)))`,
a robust estimate of the standard deviation of Gaussian noise. Returns a
scalar.

With `elementwise=true`, returns the scaled absolute deviations
`1.4826 .* abs.(t .- median(t))` with the same shape as `t`, without taking the
outer median.

!!! warning "Changed behaviour"
    Earlier versions returned the elementwise deviations by default. Use
    `MAD(t; elementwise=true)` to reproduce them.
"""
function MAD(t; elementwise::Bool=false)
    dev = 1.4826 .* abs.(t .- median(vec(t)))
    return elementwise ? dev : median(vec(dev))
end
"""
    persistence_entropy(PDs; dims = 0:1, tol = 0.0)

Compute per-dimension persistence entropy Sₚ and total persistence Eₚ.

- `PDs` :: A vector of persistence diagrams, grouped by homology degree `p`.
- `dims`:: Homology degrees (default 0:1).
- `tol` :: Discards `lifetimes ≤ tol` (guards numerical noise).

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
        Sₚ  = -sum(@view(Tₚ[:]) .* log2.(@view(Tₚ[:]) .+ eps()))  # S(T) = -∑ₚTₚ⋅log₂(Tₚ)
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
                Sτ = -sum(@view(Tτ[:]) .* log2.(@view(Tτ[:])))
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
    betti_curve(PDs, τ; dims = 0:1, interval = :closed)

Compute per-dimension Betti curves β_p(τ_j) on a user-provided grid `τ`.

- `PDs`      :: Vector of Ripserer PD (grouped by degree).
- `τ`        :: Sorted vector of filtration values (monotone increasing).
- `interval` :: Convention at the endpoints of each bar.
  - `:closed` (default, historical behaviour): a feature is alive for `b ≤ τ ≤ d`.
  - `:halfopen`: a feature is alive for `b ≤ τ < d`, the standard
    persistence-module convention `rank Hₚ(𝒦^τ)`.

The two conventions differ only on grid points that coincide exactly with a
death value (common for integer counts or capped essential deaths).

Returns `Dict{Int,Vector{Int}}` mapping p ↦ β_p(τ).
"""
function betti_curve(PDs::Vector{PersistenceDiagrams.PersistenceDiagram}, τ::AbstractVector{<:Real}; dims = 0:1,
                     interval::Symbol = :closed)
    @assert issorted(τ) "τ must be sorted ascending"
    interval in (:closed, :halfopen) ||
        throw(ArgumentError("interval must be :closed or :halfopen; got $interval."))
    β = Dict{Int, Vector{Int}}()
    for p in dims
        bd = p+1 <= length(PDs) ? PDs[p+1] : ()
        if isempty(bd)
            β[p] = zeros(Int, length(τ))
            continue
        end
        b = sort!(Float64[birth(x) for x in bd])
        d = sort!(Float64[death(x) for x in bd])
        # :closed   β(τ) = #(b ≤ τ) - #(d < τ)   alive on [b,d]
        # :halfopen β(τ) = #(b ≤ τ) - #(d ≤ τ)   alive on [b,d)
        halfopen = interval === :halfopen
        βₚ = Vector{Int}(undef, length(τ))
        @inbounds for j in eachindex(τ)
            t  = τ[j]
            nb = searchsortedlast(b, t)                                  # births ≤ t
            nd = halfopen ? searchsortedlast(d, t) : searchsortedfirst(d, t) - 1
            βₚ[j] = nb - nd
        end
        β[p] = βₚ
    end
    return β
end

"""
    betti_curvature(pd, τ; dims = 0:1, scheme = :forward, interval = :closed)

Compute β_p(τ) and κ_p(τ) = dβ_p/dτ on a grid `τ`.

- `scheme = :forward` (default) uses forward differences and sets κ[end]=κ[end-1].
  Use `:central` for central differences on interior points.
- `interval` is forwarded to [`betti_curve`](@ref).

Returns `(β::Dict{Int,Vector{Int}}, κ::Dict{Int,Vector{Float64}})`.
"""
function betti_curvature(pd, τ::AbstractVector{<:Real}; dims = 0:1, scheme::Symbol = :forward,
                         interval::Symbol = :closed)
    @assert length(τ) ≥ 2 "τ grid must have at least 2 points"
    β = betti_curve(pd, τ; dims, interval); κ = Dict{Int, Vector{Float64}}(); Δτ = diff(collect(τ))
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

#=========================================#
# PH-TDA ~ Automatic Multi-Band Tracking  # Note, Sunny.jl outputs (nω × nq) 
#                                         # Also, persistence homology analyzes Sqω 
# Sqω :: Matrix{Float64} (nq × nω)        # differently than Sωq.
# ωs  :: Vector{Float64}                  #  
#                                         # Sqω seems to be good for phonons, 
#                                         # while Sωq seems to be good for spins.  
#=========================================# 
"""
    merge_energy_bands!(peaks_sorted::Vector{Int}, ωs; merge_tol = 1.0)

Given a sorted list of peak indices for one q-slice, merge peaks whose
ω-values differ by less than merge_tol (meV).

Returns a new vector of representative peak indices.
"""
function merge_energy_bands!(peaks_sorted::Vector{Int}, ωs; merge_tol = 1.0)
    isempty(peaks_sorted) && return peaks_sorted
    merged = Int[]
    current = peaks_sorted[1]
    push!(merged, current)
    for j in 2:length(peaks_sorted)
        Δω = abs(ωs[peaks_sorted[j]] - ωs[current])
        if Δω ≤ merge_tol
            # same physical band → skip
            continue
        else
            # new physical band
            current = peaks_sorted[j]
            push!(merged, current)
        end
    end

    return merged
end
# Birth locations (indices) of the `nb` most persistent superlevel H₀ classes of
# a 1D slice, sorted by index. The essential class (global maximum) ranks first.
function _persistent_peak_indices(spec::AbstractVector{<:Real}, nb::Integer)
    y = Float64[isfinite(v) ? Float64(v) : 0.0 for v in spec]
    D0 = ripserer(Cubical(-y); dim_max=0)[1]
    isempty(D0) && return Int[]
    pers = Float64[death(x) - birth(x) for x in D0]          # Inf for the essential class
    order = sortperm(pers; rev=true)
    k = min(nb, length(order))
    idx = Int[]
    for i in order[1:k]
        v = only(Ripserer.vertices(Ripserer.birth_simplex(D0[i])))
        push!(idx, v isa CartesianIndex ? v[1] : Int(v))
    end
    return sort!(idx)
end
"""
    slice_partition_marginals(Sqω, ωs; pers_ratio_min=2.0, row_spacing=1.0, peak_selection=:height)

Partition each row of `Sqω` along the column-coordinate grid `ωs` using 0D persistent homology. Formerly `naivesplitbands`, which remains an alias. Return `(bandIq, bandIω)` with dimensions `(Nbands, size(Sqω,1))` and `(Nbands, size(Sqω,2))`, respectively.

`bandIq` integrates over columns with `Δω = ωs[2] - ωs[1]` (unchanged). `bandIω` integrates over rows with the optional `row_spacing` factor. The default `row_spacing=1.0` reproduces the historical unweighted row sum.

For physical integration, use uniformly spaced coordinates with a positive row spacing; `partition_marginals(...; grid_weighting=:physical)` validates both grids.

# What the partition counts

For each row, the superlevel `H₀` lifetimes are sorted, `ℓ₁ ≥ ℓ₂ ≥ …`, and the
local count is the position `k` of the largest consecutive ratio `ℓₖ/ℓₖ₊₁`
(or 1 if that ratio is below `pers_ratio_min`). Consequences:

- A row with exactly `n` noise-free peaks and no further bars returns at most
  `n - 1` (a single clean two-peak row returns 1), because there is no `ℓₙ₊₁`.
- In rows dominated by noise the largest ratio can fall between two noise
  bars, so the local count can be large.
- `Nbands` is the maximum local count over all rows, so it can vary from sample
  to sample and rows `b` beyond a row's local count are zero.

The output is therefore a slice-wise, persistence-gated partition of the full
2D intensity distribution rather than a set of tracked physical branches; band
index `b` is a local ordinal in each row. This is the behaviour used in the
manuscript benchmark and is preserved by default.

# Peak placement

- `peak_selection=:height` (default, historical): the `k` partition centres are
  the `k` highest local maxima of the row.
- `peak_selection=:persistence`: the centres are the birth locations of the `k`
  most persistent `H₀` classes, so low-persistence maxima on the shoulder of a
  strong band are not chosen over a weaker but persistent band. The local
  count `k` is unchanged.
"""
function slice_partition_marginals(Sqω::AbstractMatrix{<:Real}, ωs; pers_ratio_min::Real=2.0, row_spacing::Real=1.0,
                         peak_selection::Symbol=:height)
    isfinite(row_spacing) && row_spacing > 0 || throw(ArgumentError("row_spacing must be finite and positive; got $row_spacing."))
    peak_selection in (:height, :persistence) ||
        throw(ArgumentError("peak_selection must be :height or :persistence; got $peak_selection."))
    # Rows and columns exchange their physical interpretation after transposition.
    nq, nω = size(Sqω)
    nbands_local = zeros(Int, nq)

    # ================================================= #
    #          FIRST PASS: detect # bands at qᵢ         #
    # ================================================= #
    #   ∀qᵢ ∃ { bands } ⟺   [ D₀(I) ≠ 0 && {ℓ} ≠ ∅ ]    #
    #                                                   #
    # such that I ≡ I(qᵢ, ω) and D₀(I) denotes the      #
    # persistence diagram for connected components      #
    # with homological dimension p = 0.                 #
    #                                                   #
    #                                                   #
    # Note, {ℓ} = {ℓ₀, ℓ₁, ℓ₂,...} denotes the set of   #
    # lifetimes (barcodes) in D₀(I) at qᵢ.              #
    #                                                   #
    # The condition (ℓₖ/ℓₖ₊₁) ≥ ρ₀ determines the       #
    # number of bands, where ρ₀ is a constant value.    #
    #---------------------------------------------------#
    for iq in 1:nq
        spec = @view Sqω[iq, :] #spectral density
        PD = pd_array_intensities(spec; maxdim=0, superlevel=true,
                                  threshold=nothing, normalize=false)
        D0 = PD[1] #persistence diagram
        if isempty(D0)
            nbands_local[iq] = 0
            continue
        end

        ℓs = Float64[] #lifetimes
        for x in D0
            b = birth(x); d = death(x)
            isfinite(d) && d > b && push!(ℓs, d - b)
        end

        isempty(ℓs) && (nbands_local[iq] = 0; continue)
        sort!(ℓs, rev=true)
        
        #Count bands
        if length(ℓs) == 1
            nbands_local[iq] = 1
        else
            ratios = Float64[]
            for k in 1:length(ℓs)-1
                ℓs[k+1] > 0 && push!(ratios, ℓs[k] / ℓs[k+1])
            end
            if isempty(ratios)
                nbands_local[iq] = 1
            else
                k′ = argmax(ratios)
                nbands_local[iq] = (ratios[k′] ≥ pers_ratio_min ? k′ : 1)
            end
        end
    end

    # Maximum number of global bands
    Nbands = maximum(nbands_local)
    Nbands == 0 && return (bandIq = zeros(Float64, 0, nq), 
                           bandIω = zeros(Float64, 0, nω))

    # Store the ω-ranges (bands) per q so that we can integrate over q later
    bands_local = [Vector{UnitRange{Int}}() for _ in 1:nq]

    # ================================================= #
    #          SECOND PASS: determine ω-ranges          #
    # ================================================= #
    #   ∀qᵢ ∃{𝒾} : {𝒾} ≠ ∅ ⟺   ∀𝒾 ∃𝓌  : 𝓌  ∈ [ω⁻, ω⁺]   #
    #                                                   #
    # At each qᵢ check if intensity peaks {𝒾} exists.   #
    # If {𝒾} ≠ ∅, then an energy range 𝓌  must exist    #
    # such that 𝒾 ∈ 𝓌 . Tesselate the set of intensity  #
    # peaks {𝒾} using a Voronoi partition.              #
    #                                                   #
    # Returns energy ranges {𝓌 } between each peak 𝒾.   #
    #---------------------------------------------------#
    for iq in 1:nq
        spec = Sqω[iq, :]
        nb = nbands_local[iq]
        if nb == 0
            continue
        end

        # --- peak finding
        if peak_selection === :persistence
            peaks = _persistent_peak_indices(spec, nb)
        else
            peaks = Int[]
            for j in 2:nω-2
                spec[j] ≥ spec[j-1] && spec[j] ≥ spec[j+1] && spec[j] > spec[j-1] && spec[j] > spec[j+2] && push!(peaks, j)
            end
        end

        if isempty(peaks)
            # only one band covering all ω
            bands_local[iq] = [1:nω]
            continue
        end

        # choose nb strongest peaks (already the nb most persistent for :persistence)
        if peak_selection === :height
            nb = min(nb, length(peaks))
            sort!(peaks, by = j -> -spec[j])
            peaks = sort(peaks[1:nb])
        end

        # optional merging step (as in your code)
        peaks = merge_energy_bands!(peaks, ωs; merge_tol=mean(diff(ωs)))
        nb = length(peaks)

        # --- Voronoi partition between peaks
        bounds = Int[]
        for k in 1:nb-1
            mid = (peaks[k] + peaks[k+1]) ÷ 2
            push!(bounds, mid)
        end

        ranges = UnitRange{Int}[]
        start = 1
        for b in bounds
            push!(ranges, start:b)
            start = b + 1
        end
        push!(ranges, start:nω)

        bands_local[iq] = ranges
    end

    # ================================================================= #
    #          THIRD PASS: calculate integrated intensities             #
    # ================================================================= #
    #                                                                   #
    # For each band b and each energy index j:                          #
    #     bandIω[b,j] = sum over q such that j ∈ band-range(q,b)        #
    #                                                                   #
    # For each band b with band-range(q,b) = r:                         #
    #     bandIq[b,r] = sum(Sqω[qᵢ, r])*Δω                              #
    #                                                                   #
    # Generally,                                                        #
    #           Iᵇ(qₙ) = ∫I(qₙ, ω)dω : ω ∈ band-range                   #
    #           Iᵇ(ωₙ) = ∫I(q, ωₙ)dq : q ∈ band-range                   #
    #                                                                   #
    # Discretely,                                                       #
    #           Iᵇ(qₙ) = ∑ₘ I(qₙ,ωₘ)⋅Δω                                 #
    #           Iᵇ(ωₙ) = ∑ₘ I(qₘ,ωₙ)⋅row_spacing                        #
    #-------------------------------------------------------------------#

    # momentum-integrated spectrum
    bandIω = zeros(Float64, Nbands, nω)
    # energy-integrated spectrum
    bandIq = zeros(Float64, Nbands, nq)
    Δω = ωs[2] - ωs[1]
    for iq in 1:nq
        ranges = bands_local[iq]
        if isempty(ranges)
            continue
        end
        @inbounds for b in 1:length(ranges)
            r = ranges[b]
            # energy-integrated sum
            bandIq[b, iq] = sum(Sqω[iq,r])*Δω
            @inbounds for j in r
                # momentum-integrated sum
                bandIω[b, j] += Sqω[iq, j] * row_spacing
            end
        end
    end

    return (bandIq = bandIq, bandIω = bandIω)
end
"""
    partition_marginals(Iqω, ωs; pers_ratio_min=2.0, grid_weighting=:legacy, peak_selection=:height)

Apply [`slice_partition_marginals`](@ref) in both orientations using an index-space momentum axis. Formerly `autosplitbands`, which remains an alias. The default `:legacy` preserves historical numerical outputs. `peak_selection` is forwarded to [`slice_partition_marginals`](@ref).

Physical weighting requires explicit momentum coordinates via the three-argument method and is rejected here rather than assuming a momentum spacing.
"""
function partition_marginals(Iqω::AbstractMatrix{<:Real}, ωs; pers_ratio_min::Real=2.0, grid_weighting::Symbol=:legacy,
                        peak_selection::Symbol=:height)
    grid_weighting === :legacy || throw(ArgumentError("The two-argument partition_marginals method supports only grid_weighting=:legacy; supply qs and ωs for :physical."))
    nq, nω = size(Iqω)

    nq ≥ 2 || throw(ArgumentError("partition_marginals requires at least two momentum samples; got nq=$nq."))

    nω ≥ 2 || throw(ArgumentError("partition_marginals requires at least two energy samples; got nω=$nω."))

    length(ωs) == nω || throw(DimensionMismatch("Energy-axis length $(length(ωs)) does not match " *
                                                "the I(q,ω) column count nω=$nω."))

    # The legacy naive route is index based. A uniform index axis preserves
    # the momentum-grid topology without requiring a physical q coordinate.
    qs_index = collect(1.0:nq)

    return partition_marginals(Iqω, qs_index, ωs; pers_ratio_min=pers_ratio_min, peak_selection=peak_selection)
end


function _uniform_grid_spacing(coords, label::AbstractString)
    length(coords) >= 2 || throw(ArgumentError("$label grid requires at least two coordinates."))
    all(x -> x isa Real && isfinite(x), coords) || throw(ArgumentError("$label grid must contain finite real coordinates."))
    steps = diff(coords)
    all(x -> x > 0, steps) || throw(ArgumentError("$label grid must be strictly increasing."))
    spacing = Float64(steps[1])
    all(x -> isapprox(x, spacing; rtol=1e-8, atol=32 * eps(spacing)), steps) ||
        throw(ArgumentError("$label grid must be uniformly spaced for grid_weighting=:physical."))
    return spacing
end


"""
    partition_marginals(Iqω, qs, ωs; pers_ratio_min=2.0, grid_weighting=:legacy, peak_selection=:height)

Apply [`slice_partition_marginals`](@ref) in both orientations with explicit momentum `qs` and energy `ωs` coordinates, returning `(splitE, splitQ)`. Formerly `autosplitbands`, which remains an alias. `grid_weighting=:legacy` preserves historical marginal amplitudes. `peak_selection` is forwarded to [`slice_partition_marginals`](@ref).

Note on `:legacy` weighting: each orientation always multiplies by the spacing of its own column coordinate (`ωs[2]-ωs[1]` for `splitE.bandIq`, `qs[2]-qs[1]` for `splitQ.bandIω`) but uses a row spacing of 1. The two-argument method passes `qs = 1:nq`, so its `splitQ.bandIω` differs from the three-argument `:legacy` result by the factor `Δq`. Complete-array L¹ normalization removes these constant factors.

With `grid_weighting=:physical`, each marginal is integrated using the spacing of its summed coordinate, as in Eqs. (7)–(8) of the manuscript. Physical mode requires strictly increasing, finite, uniformly spaced coordinate grids.
"""
function partition_marginals(Iqω::AbstractMatrix{<:Real}, qs, ωs; pers_ratio_min::Real=2.0, grid_weighting::Symbol=:legacy,
                        peak_selection::Symbol=:height)
    grid_weighting in (:legacy, :physical) || throw(ArgumentError("grid_weighting must be :legacy or :physical; got $grid_weighting."))
    nq, nω = size(Iqω)

    nq ≥ 2 || throw(ArgumentError("partition_marginals requires at least two momentum samples; got nq=$nq."))

    nω ≥ 2 || throw(ArgumentError("partition_marginals requires at least two energy samples; got nω=$nω."))

    length(qs) == nq || throw(DimensionMismatch("Momentum-axis length $(length(qs)) does not match " *
                                                "the I(q,ω) row count nq=$nq."))

    length(ωs) == nω || throw(DimensionMismatch("Energy-axis length $(length(ωs)) does not match " *
                                                "the I(q,ω) column count nω=$nω."))

    Δq = grid_weighting === :physical ? _uniform_grid_spacing(qs, "Momentum") : 1.0
    Δω = grid_weighting === :physical ? _uniform_grid_spacing(ωs, "Energy") : 1.0

    # Each call weights its rows by the other coordinate's physical spacing.
    splitω = slice_partition_marginals(Iqω, ωs; pers_ratio_min=pers_ratio_min, row_spacing=Δq, peak_selection=peak_selection)
    splitq = slice_partition_marginals(permutedims(Iqω), qs; pers_ratio_min=pers_ratio_min, row_spacing=Δω, peak_selection=peak_selection)

    splitE = (bandIq = splitω.bandIq, bandIω = splitω.bandIω)

    # Restore the original physical naming after transposition.
    splitQ = (bandIq = splitq.bandIω, bandIω = splitq.bandIq)

    return (splitE = splitE, splitQ = splitQ)
end

#---------------------------------------------------------#
# Historical names (aliases of the same function objects) #
#---------------------------------------------------------#
"""
    naivesplitbands(Sqω, ωs; kwargs...)

Historical name of [`slice_partition_marginals`](@ref); the same function.
"""
const naivesplitbands = slice_partition_marginals

"""
    autosplitbands(Iqω, [qs,] ωs; kwargs...)

Historical name of [`partition_marginals`](@ref); the same function.
"""
const autosplitbands = partition_marginals
