#-------------------------------------#
# Resolution-Matched Rebinning        #
#-------------------------------------#
# Coarsens an I(q, ω) map to blocks of about one resolution width before segmentation, and maps
# masks on the rebinned grid back to the original bins (e.g. parent-band supports for the QFI).

"""
    RebinResult

Result of [`rebin_to_resolution`](@ref).

Fields:

- `I`, `σ` :: rebinned intensity and uncertainty (`σ === nothing` when none was given);
  blocks without measured bins hold `fillvalue` and `σ = NaN`.
- `valid` :: `BitMatrix`, true where the block contains at least one measured bin.
- `q`, `ω` :: block centres (means of the original bin centres in each block).
- `qfactor`, `ωfactor` :: block size in original bins.
- `resolution_q`, `resolution_ω` :: resolution widths used to choose the block size.
- `source_size` :: `(nq, nω)` of the original grid. Trailing bins that do not fill a whole
  block (`source_size .- (qfactor, ωfactor) .* size(I)`) are dropped.
- `fillvalue`, `nfilled` :: value given to unmeasured bins and how many original bins had it.
- `average` :: `:valid` or `:all` (see [`rebin_to_resolution`](@ref)).
"""
struct RebinResult
    I::Matrix{Float64}
    σ::Union{Nothing,Matrix{Float64}}
    valid::BitMatrix
    q::Vector{Float64}
    ω::Vector{Float64}
    qfactor::Int
    ωfactor::Int
    resolution_q::Float64
    resolution_ω::Float64
    source_size::Tuple{Int,Int}
    fillvalue::Float64
    nfilled::Int
    average::Symbol
end

function Base.show(io::IO, r::RebinResult)
    print(io, "RebinResult($(r.source_size[1])×$(r.source_size[2]) → $(size(r.I, 1))×$(size(r.I, 2)), ",
          "factors=($(r.qfactor), $(r.ωfactor)), Γq=$(round(r.resolution_q; sigdigits=4)), ",
          "Γω=$(round(r.resolution_ω; sigdigits=4)), average=:$(r.average))")
end

"""
    rebin_to_resolution(I, qs, ωs; σ=nothing, valid=nothing, resolution_q=nothing,
                        resolution_ω=nothing, factors=nothing, fill=:noise_floor,
                        average=:valid) -> RebinResult

Average `I` (`nq × nω`, rows = momentum, columns = energy) over blocks of `qfactor × ωfactor`
bins, where each factor is the resolution width in bins, `clamp(round(Γ/Δ), 1, n)`.

The widths `Γq`, `Γω` come from `resolution_q`, `resolution_ω` when given, and are otherwise
estimated from the data as in [`adaptive_autosplitbands2d`](@ref). `factors = (qfactor, ωfactor)`
sets the block size directly. Trailing bins that do not fill a block are dropped.

Unmeasured bins are those with `valid[i, j] == false`; by default every non-finite `I` (or `σ`).
They are given the value set by `fill`: `:noise_floor` (the smallest positive measured value) or
`:zero`. The resolution estimate is made on this filled map.

`average` selects how a block is averaged:

- `:valid` (default) averages only the measured bins, so unmeasured bins do not bias a block.
  A block with no measured bins holds the fill value.
- `:all` averages every bin of the filled map, which reproduces the manuscript workflow
  (fill first, then block means).

With `σ`, the block uncertainty is `sqrt(Σ σ²) / n` over the measured bins, where `n` is the number
of bins averaged; blocks without measured bins have `σ = NaN`.

Use [`upsample_mask`](@ref) or `band_supports(result, rebin)` to map masks on the rebinned grid
back to the original bins.
"""
function rebin_to_resolution(I::AbstractMatrix{<:Real}, qs::AbstractVector{<:Real}, ωs::AbstractVector{<:Real};
                             σ::Union{Nothing,AbstractMatrix{<:Real}}=nothing,
                             valid::Union{Nothing,AbstractMatrix{Bool}}=nothing,
                             resolution_q=nothing, resolution_ω=nothing,
                             factors::Union{Nothing,Tuple{Integer,Integer}}=nothing,
                             fill::Symbol=:noise_floor, average::Symbol=:valid)
    nq, nω = size(I)
    length(qs) == nq || throw(DimensionMismatch("qs has $(length(qs)) entries; I has $nq rows."))
    length(ωs) == nω || throw(DimensionMismatch("ωs has $(length(ωs)) entries; I has $nω columns."))
    nq > 1 && nω > 1 || throw(ArgumentError("at least two bins are needed along each axis."))
    σ === nothing || size(σ) == size(I) || throw(DimensionMismatch("σ must have the size of I."))
    fill in (:noise_floor, :zero) || throw(ArgumentError("fill must be :noise_floor or :zero; got $fill."))
    average in (:valid, :all) || throw(ArgumentError("average must be :valid or :all; got $average."))

    if valid === nothing
        mask = BitMatrix(isfinite.(I))
        σ === nothing || (mask .&= isfinite.(σ))
    else
        size(valid) == size(I) || throw(DimensionMismatch("valid must have the size of I."))
        mask = BitMatrix(valid)
        any(i -> mask[i] && !isfinite(I[i]), eachindex(I)) &&
            throw(ArgumentError("valid marks non-finite intensities as measured."))
    end

    Ifill = Matrix{Float64}(I)
    fillvalue = if fill === :noise_floor
        pos = [Float64(I[i]) for i in eachindex(I) if mask[i] && I[i] > 0]
        isempty(pos) ? 0.0 : minimum(pos)
    else
        0.0
    end
    nfilled = count(!, mask)
    Ifill[.!mask] .= fillvalue

    qv, ωv = Float64.(qs), Float64.(ωs)
    Δq, Δω = mean(abs.(diff(qv))), mean(abs.(diff(ωv)))
    _, _, Γq, Γω = adaptive_resolution_radii(qv, ωv, Ifill; resolution_q=resolution_q, resolution_ω=resolution_ω)
    qf, ωf = if factors === nothing
        (clamp(floor(Int, Γq / Δq + 0.5), 1, nq), clamp(floor(Int, Γω / Δω + 0.5), 1, nω))
    else
        all(>(0), factors) || throw(ArgumentError("factors must be positive; got $factors."))
        (min(Int(factors[1]), nq), min(Int(factors[2]), nω))
    end
    mq, mω = fld(nq, qf), fld(nω, ωf)

    Ib = Matrix{Float64}(undef, mq, mω)
    σb = σ === nothing ? nothing : Matrix{Float64}(undef, mq, mω)
    vb = falses(mq, mω)
    qb = [mean(@view qv[(i-1)*qf+1:i*qf]) for i in 1:mq]
    ωb = [mean(@view ωv[(j-1)*ωf+1:j*ωf]) for j in 1:mω]
    @inbounds for i in 1:mq, j in 1:mω
        rq, rω = (i-1)*qf+1:i*qf, (j-1)*ωf+1:j*ωf
        nvalid = count(@view mask[rq, rω])
        vb[i, j] = nvalid > 0
        if average === :all
            Ib[i, j] = mean(@view Ifill[rq, rω])
            n = qf * ωf
        elseif nvalid > 0
            s = 0.0
            for jj in rω, ii in rq
                mask[ii, jj] && (s += Ifill[ii, jj])
            end
            Ib[i, j] = s / nvalid
            n = nvalid
        else
            Ib[i, j] = fillvalue
            n = 0
        end
        if σb !== nothing
            if nvalid == 0
                σb[i, j] = NaN
            else
                v = 0.0
                for jj in rω, ii in rq
                    mask[ii, jj] && (v += Float64(σ[ii, jj])^2)
                end
                σb[i, j] = sqrt(v) / n
            end
        end
    end
    return RebinResult(Ib, σb, vb, qb, ωb, qf, ωf, Γq, Γω, (nq, nω), fillvalue, nfilled, average)
end

"""
    upsample_mask(M, rebin::RebinResult) -> BitMatrix

Map a mask `M` on the rebinned grid back to the original grid: every original bin in a block
takes the block's value. Bins dropped by the rebinning are `false`. Disjoint masks stay disjoint.
"""
function upsample_mask(M::AbstractMatrix{Bool}, rebin::RebinResult)
    size(M) == size(rebin.I) ||
        throw(DimensionMismatch("mask is $(size(M)); the rebinned grid is $(size(rebin.I))."))
    qf, ωf = rebin.qfactor, rebin.ωfactor
    F = falses(rebin.source_size)
    @inbounds for j in axes(M, 2), i in axes(M, 1)
        M[i, j] && (F[(i-1)*qf+1:i*qf, (j-1)*ωf+1:j*ωf] .= true)
    end
    return F
end
