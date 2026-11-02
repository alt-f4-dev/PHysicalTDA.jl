#=
Quantum Fisher information (QFI) from the dynamical structure factor.

Conventions (Supplemental Material, Sec. S2):

    F_Q^αα(q, T) = ∫₀^∞ K_T(ω) S^αα(q, ω) dω,
    K_T(ω)       = 4 tanh(βω/2) (1 − e^{−βω}),     β = 1/(k_B T),

with S^αα normalized so that ∫ S^αα(q, ω) dω = ⟨Ô_q† Ô_q⟩ and inversion symmetry
assumed. Summing over components, F_Q = Σ_α F_Q^αα and S_mag = Σ_α S^αα, and

    nQFI(q, T) = F_Q(q, T) / (4 n_c S² N),

with n_c = 3 components for the component-summed response (12 S² N). nQFI > m
certifies at least (m + 1)-partite entanglement.

Band-resolved contributions restrict the ω > 0 integral to disjoint supports
C^(b) (Eq. S29); F_Q = Σ_b F_Q^(b) + F_Q^(res) holds exactly for the measured bins.
=#

"""
    KB_MEV_PER_K

Boltzmann constant in meV/K (CODATA, exact): `0.08617333262`.
"""
const KB_MEV_PER_K = 0.08617333262

#----------------#
# Kernel / grids #
#----------------#

"""
    qfi_kernel(ω, T; kB=KB_MEV_PER_K)

QFI kernel `K_T(ω) = 4 tanh(βω/2) (1 − e^{−βω})` with `β = 1/(kB*T)`.

`ω` is an energy transfer and `T` a temperature; with the default `kB` they are in
meV and K. The kernel is the positive-frequency form, so it returns `0` for
`ω ≤ 0`. It satisfies `0 ≤ K_T ≤ 4`, tends to `4` as `T → 0` (`T == 0` returns
exactly `4` for `ω > 0`), and behaves as `2(βω)²` for `βω ≪ 1`.
"""
function qfi_kernel(ω::Real, T::Real; kB::Real=KB_MEV_PER_K)
    T >= 0 || throw(ArgumentError("temperature must be non-negative; got T = $T."))
    kB > 0 || throw(ArgumentError("kB must be positive; got kB = $kB."))
    ω > 0 || return 0.0
    T == 0 && return 4.0
    x = Float64(ω) / (Float64(kB) * Float64(T))
    return 4.0 * tanh(x / 2) * (-expm1(-x))
end

"""
    energy_bin_widths(ωs)

Widths `Δω_j` of energy bins whose centres are `ωs` (strictly increasing).

Interior widths are `(ω_{j+1} − ω_{j−1})/2`; the first and last bins take the
spacing to their single neighbour. For a uniform grid every width equals the
spacing.
"""
function energy_bin_widths(ωs)
    n = length(ωs)
    n >= 2 || throw(ArgumentError("need at least two energy bins to infer widths."))
    w = Vector{Float64}(undef, n)
    @inbounds for j in 1:n-1
        ωs[j+1] > ωs[j] || throw(ArgumentError("energy centres must be strictly increasing."))
    end
    w[1] = ωs[2] - ωs[1]
    w[n] = ωs[n] - ωs[n-1]
    @inbounds for j in 2:n-1
        w[j] = (ωs[j+1] - ωs[j-1]) / 2
    end
    return w
end

function _energy_weights(ωs, Δω)
    if Δω === nothing
        return energy_bin_widths(ωs)
    elseif Δω isa Real
        Δω > 0 || throw(ArgumentError("Δω must be positive; got $Δω."))
        return fill(Float64(Δω), length(ωs))
    else
        length(Δω) == length(ωs) ||
            throw(DimensionMismatch("Δω has length $(length(Δω)); expected $(length(ωs))."))
        all(>(0), Δω) || throw(ArgumentError("all bin widths Δω must be positive."))
        return Float64.(collect(Δω))
    end
end

# Columns that enter the positive-frequency integral.
function _included_columns(ωs, ωrange)
    lo, hi = ωrange
    lo <= hi || throw(ArgumentError("ωrange must satisfy lo ≤ hi; got $ωrange."))
    return [ω > max(lo, 0) && ω <= hi for ω in ωs]
end

function _check_grid(S, ωs)
    size(S, 2) == length(ωs) ||
        throw(DimensionMismatch("S has $(size(S, 2)) energy columns but ωs has length $(length(ωs))."))
end

function _check_same_size(A, S, name)
    A === nothing && return nothing
    size(A) == size(S) ||
        throw(DimensionMismatch("$name has size $(size(A)); expected $(size(S))."))
    return nothing
end

"""
    qfi_density(S, ωs, T; kB=KB_MEV_PER_K)

QFI spectral density `ℵ_T(q, ω) = K_T(ω) S(q, ω)` for an `nq × nω` matrix `S` whose
columns are the energies `ωs`. Columns with `ω ≤ 0` are zero; `NaN` (unmeasured)
bins stay `NaN`.
"""
function qfi_density(S::AbstractMatrix{<:Real}, ωs, T::Real; kB::Real=KB_MEV_PER_K)
    _check_grid(S, ωs)
    K = [qfi_kernel(ω, T; kB=kB) for ω in ωs]
    out = Matrix{Float64}(undef, size(S))
    @inbounds for j in axes(S, 2), i in axes(S, 1)
        s = Float64(S[i, j])
        out[i, j] = isnan(s) ? NaN : K[j] * s
    end
    return out
end

#-----------#
# Total QFI #
#-----------#

"""
    QFIResult

Momentum-resolved QFI returned by [`qfi`](@ref).

- `FQ`    :: `F_Q(q_i, T)` for each row of `S`.
- `σFQ`   :: propagated one-standard-deviation uncertainty (zeros when no `σ` was given).
- `nbins` :: number of measured (non-`NaN`) bins entering each `F_Q(q_i, T)`.
- `T`     :: temperature used in the kernel.
"""
struct QFIResult
    FQ::Vector{Float64}
    σFQ::Vector{Float64}
    nbins::Vector{Int}
    T::Float64
end

"""
    qfi(S, ωs, T; σ=nothing, mask=nothing, ωrange=(0, Inf), Δω=nothing, kB=KB_MEV_PER_K) -> QFIResult

Momentum-resolved QFI `F_Q(q_i, T) = Σ_j K_T(ω_j) S(q_i, ω_j) Δω_j` over the bins with
`ω_j > max(ωrange[1], 0)` and `ω_j ≤ ωrange[2]`.

- `S`      :: `nq × nω` dynamical structure factor; rows are momenta, columns the energies `ωs`.
              Pass the component-summed `S_mag` for the component-summed QFI.
- `σ`      :: optional `nq × nω` standard deviations of `S`; bins are treated as
              independent, so `σFQ² = Σ_j (K_T Δω_j σ_ij)²`.
- `mask`   :: optional `nq × nω` Bool matrix; only `true` bins are summed.
- `ωrange` :: energy window; set `ωrange[1]` above the elastic line to exclude it.
- `Δω`     :: bin widths (scalar or vector); inferred from `ωs` by default
              ([`energy_bin_widths`](@ref)).

Bins where `S` (or `σ`) is `NaN` are treated as unmeasured and skipped. Because
`K_T ≥ 0`, restricting the sum to measured bins of a non-negative magnetic response
gives a lower bound on the full integral; subtracted backgrounds that leave negative
bins are summed as they are.
"""
function qfi(S::AbstractMatrix{<:Real}, ωs, T::Real;
             σ::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
             mask::Union{Nothing,AbstractMatrix{Bool}}=nothing,
             ωrange::Tuple{<:Real,<:Real}=(0, Inf),
             Δω=nothing,
             kB::Real=KB_MEV_PER_K)
    _check_grid(S, ωs)
    _check_same_size(σ, S, "σ")
    _check_same_size(mask, S, "mask")
    K = [qfi_kernel(ω, T; kB=kB) for ω in ωs]
    w = _energy_weights(ωs, Δω)
    cols = _included_columns(ωs, ωrange)
    nq = size(S, 1)
    FQ = zeros(Float64, nq); var = zeros(Float64, nq); nb = zeros(Int, nq)
    @inbounds for (jj, j) in enumerate(axes(S, 2))
        cols[jj] || continue
        c = K[jj] * w[jj]
        for (ii, i) in enumerate(axes(S, 1))
            mask === nothing || mask[i, j] || continue
            s = Float64(S[i, j])
            isnan(s) && continue
            e = σ === nothing ? 0.0 : Float64(σ[i, j])
            isnan(e) && continue
            FQ[ii] += c * s
            var[ii] += (c * e)^2
            nb[ii] += 1
        end
    end
    return QFIResult(FQ, sqrt.(var), nb, Float64(T))
end

#-------------------------#
# Normalization & witness #
#-------------------------#

function _qfi_bound_scale(S_spin::Real, N::Real, ncomponents::Integer)
    S_spin > 0 || throw(ArgumentError("spin S must be positive; got $S_spin."))
    N > 0 || throw(ArgumentError("N must be positive; got $N."))
    1 <= ncomponents <= 3 || throw(ArgumentError("ncomponents must be 1, 2 or 3; got $ncomponents."))
    return 4 * ncomponents * Float64(S_spin)^2 * Float64(N)
end

"""
    nqfi(FQ; S, N=1, ncomponents=3)
    nqfi(r::QFIResult; S, N=1, ncomponents=3)

Normalized QFI `F_Q / (4 n_c S² N)`. With `ncomponents = 3` (component-summed
response) this is `F_Q / (12 S² N)` (Eq. S28).

Use `N = 1` when `S(q, ω)` is normalized per spin (per magnetic ion); use the number
of spins when it is normalized to the whole system.
"""
nqfi(FQ::AbstractVector{<:Real}; S::Real, N::Real=1, ncomponents::Integer=3) =
    Float64.(FQ) ./ _qfi_bound_scale(S, N, ncomponents)
nqfi(r::QFIResult; kwargs...) = nqfi(r.FQ; kwargs...)

"""
    nqfi_lower_bound(r::QFIResult; S, N=1, ncomponents=3, nσ=2)

Conservative normalized QFI `(F_Q − nσ σ_FQ) / (4 n_c S² N)`.
"""
function nqfi_lower_bound(r::QFIResult; S::Real, N::Real=1, ncomponents::Integer=3, nσ::Real=2)
    nσ >= 0 || throw(ArgumentError("nσ must be non-negative; got $nσ."))
    return (r.FQ .- nσ .* r.σFQ) ./ _qfi_bound_scale(S, N, ncomponents)
end

"""
    qfi_witness(r::QFIResult; S, m=1, N=1, ncomponents=3, nσ=2)

Entanglement-witness criterion `G_m(q, T)` (Eq. S30). Returns a named tuple

- `nqfi_lb`   :: conservative normalized QFI ([`nqfi_lower_bound`](@ref)).
- `certified` :: `G_m = 1` (`true`) where `nqfi_lb > m`, i.e. at least
                 `(m + 1)`-partite entanglement is certified at that momentum.
- `m`         :: the depth threshold used.

`false` means only that entanglement was not certified.
"""
function qfi_witness(r::QFIResult; S::Real, m::Integer=1, N::Real=1,
                     ncomponents::Integer=3, nσ::Real=2)
    m >= 1 || throw(ArgumentError("m must be a positive integer; got $m."))
    lb = nqfi_lower_bound(r; S=S, N=N, ncomponents=ncomponents, nσ=nσ)
    return (nqfi_lb=lb, certified=BitVector(lb .> m), m=Int(m))
end

#---------------#
# Band supports #
#---------------#

"""
    band_supports(result, nq, nω; support=:components)
    band_supports(domains, nq, nω)

One `nq × nω` `BitMatrix` per parent band `b`, marking the support `C^(b)` used for
band-resolved QFI.

For an [`adaptive_autosplitbands2d`](@ref) result:

- `support = :components` (default): the union over `ℓ` of the retained connected
  components `C^(b,ℓ)` of parent `b` (exact pixel membership), as in Eq. S29.
- `support = :domains`: the q-first energy windows `result.qfirst.domains`.

A domain list is indexed as `domains[iq][b]`, an energy-index range for parent `b`
at momentum row `iq`.
"""
function band_supports(result, nq::Integer, nω::Integer; support::Symbol=:components)
    if support === :components
        nb = result.qfirst.nparents
        masks = [falses(nq, nω) for _ in 1:nb]
        for c in result.components
            1 <= c.parent_band <= nb || continue
            M = masks[c.parent_band]
            for li in c.linear_indices
                1 <= li <= nq * nω ||
                    throw(DimensionMismatch("component pixel $li lies outside an $nq × $nω grid; " *
                                            "S must be on the grid the segmentation was run on."))
                M[li] = true
            end
        end
        return masks
    elseif support === :domains
        return band_supports(result.qfirst.domains, nq, nω)
    else
        throw(ArgumentError("support must be :components or :domains; got $support."))
    end
end

function band_supports(domains::AbstractVector, nq::Integer, nω::Integer)
    length(domains) == nq ||
        throw(DimensionMismatch("domains has $(length(domains)) rows; expected $nq."))
    nb = maximum((length(d) for d in domains); init=0)
    masks = [falses(nq, nω) for _ in 1:nb]
    for iq in 1:nq, (b, r) in enumerate(domains[iq])
        isempty(r) && continue
        (first(r) >= 1 && last(r) <= nω) ||
            throw(DimensionMismatch("domain range $r at row $iq lies outside 1:$nω."))
        masks[b][iq, r] .= true
    end
    return masks
end

#--------------------#
# Band-resolved QFI  #
#--------------------#

"""
    BandQFIResult

Band-resolved QFI returned by [`band_qfi`](@ref). With `nb` bands and `nq` momenta:

- `total` :: [`QFIResult`](@ref) over all included measured bins.
- `FQ`, `σFQ` :: `nb × nq` band contributions `F_Q^(b)(q, T)` and uncertainties.
- `FQres`, `σres` :: residual `F_Q^(res)` from included bins outside every support.
- `ξ`, `σξ` :: `nb × nq` relative contributions `F_Q^(b) / F_Q` (Eq. S31) and their
  uncertainties; `NaN` where `F_Q ≤ 0`. The denominator is the QFI over the
  measured domain, i.e. `ξ_b^obs` when only part of the spectrum was measured.

`total.FQ == vec(sum(FQ; dims=1)) .+ FQres` up to rounding.
"""
struct BandQFIResult
    total::QFIResult
    FQ::Matrix{Float64}
    σFQ::Matrix{Float64}
    FQres::Vector{Float64}
    σres::Vector{Float64}
    ξ::Matrix{Float64}
    σξ::Matrix{Float64}
end

"""
    band_qfi(S, supports, ωs, T; σ=nothing, ωrange=(0, Inf), Δω=nothing, kB=KB_MEV_PER_K)
    band_qfi(S, result, ωs, T; support=:components, kwargs...)

Band-resolved QFI (Eq. S29) for disjoint supports `C^(b)`, given either as a vector
of `nq × nω` Bool masks or as an [`adaptive_autosplitbands2d`](@ref) result (see
[`band_supports`](@ref)). Keywords are as for [`qfi`](@ref).

The decomposition is linear in `S` for fixed supports and exact over the measured
bins: `F_Q = Σ_b F_Q^(b) + F_Q^(res)`. Overlapping supports raise an error.

The supports must be on the grid of `S`. The segmentation itself is
data-dependent, so it is usually run on the intensity map and the supports are
then applied to the converted `S_mag` on the same grid.
"""
function band_qfi(S::AbstractMatrix{<:Real}, supports::AbstractVector{<:AbstractMatrix{Bool}},
                  ωs, T::Real;
                  σ::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
                  ωrange::Tuple{<:Real,<:Real}=(0, Inf),
                  Δω=nothing,
                  kB::Real=KB_MEV_PER_K)
    _check_grid(S, ωs)
    for (b, M) in enumerate(supports)
        size(M) == size(S) ||
            throw(DimensionMismatch("support $b has size $(size(M)); expected $(size(S))."))
    end
    coverage = zeros(Int, size(S))
    for M in supports
        coverage .+= M
    end
    any(>(1), coverage) && throw(ArgumentError("band supports overlap; band-resolved QFI requires disjoint supports."))

    kw = (σ=σ, ωrange=ωrange, Δω=Δω, kB=kB)
    total = qfi(S, ωs, T; kw...)
    nb, nq = length(supports), size(S, 1)
    FQ = zeros(Float64, nb, nq); σFQ = zeros(Float64, nb, nq)
    for (b, M) in enumerate(supports)
        rb = qfi(S, ωs, T; mask=M, kw...)
        FQ[b, :] .= rb.FQ
        σFQ[b, :] .= rb.σFQ
    end
    rres = qfi(S, ωs, T; mask=coverage .== 0, kw...)

    ξ = fill(NaN, nb, nq); σξ = fill(NaN, nb, nq)
    @inbounds for i in 1:nq
        F = total.FQ[i]
        F > 0 || continue
        varF = total.σFQ[i]^2
        for b in 1:nb
            x = FQ[b, i] / F
            ξ[b, i] = x
            # ξ = a / (a + rest) with a and rest from disjoint (independent) bins.
            vara = σFQ[b, i]^2
            varrest = max(varF - vara, 0.0)
            σξ[b, i] = sqrt((1 - x)^2 * vara + x^2 * varrest) / F
        end
    end
    return BandQFIResult(total, FQ, σFQ, rres.FQ, rres.σFQ, ξ, σξ)
end

function band_qfi(S::AbstractMatrix{<:Real}, result::NamedTuple, ωs, T::Real;
                  support::Symbol=:components, kwargs...)
    masks = band_supports(result, size(S, 1), size(S, 2); support=support)
    return band_qfi(S, masks, ωs, T; kwargs...)
end

"""
    band_dependence(br::BandQFIResult; S, m=1, N=1, ncomponents=3, nσ=2)

Leave-one-band-out test of a certified witness. For each band `b`,

    nQFI_lb^(−b) = nQFI_lb − F_Q^(b) / (4 n_c S² N),

using the point estimate `F_Q^(b)`. Returns a named tuple

- `nqfi_lb`           :: `nQFI_lb(q)` of the full measured response.
- `certified`         :: `nqfi_lb .> m`.
- `nqfi_lb_without`   :: `nb × nq` values of `nQFI_lb^(−b)`.
- `certified_without` :: `nb × nq`, whether the certificate survives removing band `b`.

A certificate that is lost when `b` is removed depends on that band. This
attributes contribution, not an entanglement depth per band.
"""
function band_dependence(br::BandQFIResult; S::Real, m::Integer=1, N::Real=1,
                         ncomponents::Integer=3, nσ::Real=2)
    w = qfi_witness(br.total; S=S, m=m, N=N, ncomponents=ncomponents, nσ=nσ)
    scale = _qfi_bound_scale(S, N, ncomponents)
    without = w.nqfi_lb' .- br.FQ ./ scale
    return (nqfi_lb=w.nqfi_lb, certified=w.certified,
            nqfi_lb_without=Matrix(without), certified_without=BitMatrix(without .> m))
end

#--------------------------------------#
# Converting intensity to S_mag(q, ω)  #
#--------------------------------------#

"""
    MN2_J0

Dipole-approximation `⟨j₀⟩` coefficients `(A, a, B, b, C, c, D)` for Mn²⁺
(P. J. Brown, International Tables for Crystallography Vol. C, Sec. 4.4.5):
`(0.4220, 17.6840, 0.5948, 6.0050, 0.0043, -0.6090, -0.0219)`.
"""
const MN2_J0 = (0.4220, 17.6840, 0.5948, 6.0050, 0.0043, -0.6090, -0.0219)

"""
    magnetic_form_factor(Q, coeffs=MN2_J0)

Dipole-approximation magnetic form factor
`⟨j₀⟩(s) = A e^{−a s²} + B e^{−b s²} + C e^{−c s²} + D` with `s = Q/(4π)` and `Q` in Å⁻¹.
`coeffs = (A, a, B, b, C, c, D)`; the default is Mn²⁺ ([`MN2_J0`](@ref)).
"""
function magnetic_form_factor(Q::Real, coeffs=MN2_J0)
    A, a, B, b, C, c, D = coeffs
    s2 = (Float64(Q) / (4π))^2
    return A * exp(-a * s2) + B * exp(-b * s2) + C * exp(-c * s2) + D
end

"""
    to_smag(I, Qmag; σ=nothing, form_factor=Q -> magnetic_form_factor(Q), g=2,
            polarization=3/2, scale=1) -> (S, σ)

Convert a background-subtracted magnetic intensity `I(q, ω)` (`nq × nω`) into the
component-summed `S_mag(q, ω)`:

    S_mag = scale × polarization × I / ((g/2)² F(|Q|)²).

- `Qmag`         :: `|Q|` in Å⁻¹, either one value per row (`nq`) or one per bin (`nq × nω`).
- `form_factor`  :: function of `|Q|` returning the magnetic form factor.
- `g`            :: Landé g-factor.
- `polarization` :: ratio of the component-summed response to the unpolarized cross
                    section. Unpolarized neutrons measure `Σ_α (1 − Q̂_α²) S^αα`; for an
                    isotropic response this is `2/3` of `S_mag`, hence the default `3/2`.
                    Supply the appropriate factor for an anisotropic response.
- `scale`        :: absolute-normalization factor. With the unpolarized magnetic cross
                    section `d²σ/dΩdE = (γr₀)² (k_f/k_i) [(g/2)F(Q)]² e^{−2W} Σ_αβ (δ_αβ − Q̂_αQ̂_β) S^αβ`,
                    `scale` converts `I` (in absolute units) to `S` per magnetic ion: divide by
                    `(γr₀)² ≈ 0.2905 b`, by `k_f/k_i` and `e^{−2W}` if not already removed, and by
                    the number of magnetic ions in the normalization volume.

Returns the converted `S` and, when `σ` is given, its uncertainties (scaled the same
way; the uncertainty of `scale` itself is not propagated). `NaN` bins stay `NaN`.
Applying the conversion is the caller's choice: none of the QFI functions apply it
implicitly.
"""
function to_smag(I::AbstractMatrix{<:Real}, Qmag;
                 σ::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
                 form_factor=Q -> magnetic_form_factor(Q),
                 g::Real=2, polarization::Real=3/2, scale::Real=1)
    _check_same_size(σ, I, "σ")
    nq, nω = size(I)
    Qm = if Qmag isa AbstractVector
        length(Qmag) == nq || throw(DimensionMismatch("Qmag has length $(length(Qmag)); expected $nq."))
        repeat(Float64.(collect(Qmag)), 1, nω)
    else
        size(Qmag) == size(I) || throw(DimensionMismatch("Qmag has size $(size(Qmag)); expected $(size(I))."))
        Float64.(Qmag)
    end
    factor = similar(Qm)
    @inbounds for k in eachindex(Qm)
        F = form_factor(Qm[k])
        factor[k] = scale * polarization / ((g / 2)^2 * F^2)
    end
    S = Float64.(I) .* factor
    σS = σ === nothing ? nothing : abs.(Float64.(σ) .* factor)
    return (S, σS)
end
