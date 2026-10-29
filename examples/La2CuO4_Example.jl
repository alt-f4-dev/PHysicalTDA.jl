# LSWT spectra of La₂CuO₄, a 4D intensity cube, and its persistence summaries.
#
# Run from an environment containing PHysicalTDA and Sunny:
#     julia --project examples/La2CuO4_Example.jl
using PHysicalTDA, Sunny, Statistics

#--------------------------------#
# Load Crystal and LSWT Spectrum #
#--------------------------------#
cryst, sys = La2CuO4()

# q-paths through the zone centre along h, k and ℓ
H = q_space_path(cryst, [[h, 0, 0] for h in -1:1], 500)
K = q_space_path(cryst, [[0, k, 0] for k in -1:1], 500)
L = q_space_path(cryst, [[0, 0, ℓ] for ℓ in -1:1], 500)
qpaths = @name [H, K, L]
LCO = LSWT(cryst, sys, qpaths.paths)
ωaxis = LCO.energies                      # shared energy axis of every path

#------------------------------#
# Construct 4D Intensity Field #
#------------------------------#
Qs = [[h, k, ℓ] for h in -1:0.2:1, k in -1:0.2:1, ℓ in -1:0.2:1]
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
result2D = intensities(swt, vec(Qs); energies=ωaxis, kernel=gaussian(fwhm=35))
result4D = convert4D(result2D, Qs, ωaxis)  # axes (h, k, ℓ, ω)

hkl = collapse(result4D; over=:ω, op=sum)          # energy-integrated (h,k,ℓ)
hω  = collapse(result4D; over=(:k, :ℓ), op=mean)   # (h, ω) projection

#------------------------------#
# Compute Persistence Diagrams #
#------------------------------#
PD = pd_array_intensities(hω; maxdim=1, superlevel=true, normalize=true)

#----------------------------#
# Topological Fingerprinting #
#----------------------------#
Sₚ, Eₚ = persistence_entropy(PD; dims=0:1, tol=0.0)

# τ grid spanning the finite births and deaths in PD
function taugrid(PD; N::Int=256)
    births = Float64[]; deaths = Float64[]
    for D in PD
        isempty(D) && continue
        append!(births, filter(isfinite, Float64.(birth.(D))))
        append!(deaths, filter(isfinite, Float64.(death.(D))))
    end
    @assert !isempty(births) && !isempty(deaths) "Empty PD: no features!"
    tmin, tmax = minimum(births), maximum(deaths)
    δ = 1e-9 * max(1.0, abs(tmax - tmin))
    return collect(range(tmin + δ, tmax - δ; length=N))
end
τ = taugrid(PD)

β, κ = betti_curvature(PD, τ; dims=0:1, scheme=:forward)

#-------------------------------------#
# Validation Check: Eₚ ≈ ∫βₚ(τ) dτ    #
#-------------------------------------#
function trapz(x, y)
    s = 0.0
    for i in 1:length(x)-1
        s += (x[i+1] - x[i]) * (y[i+1] + y[i]) / 2
    end
    return s
end

for p in 0:1
    est = trapz(τ, float.(β[p]))
    ref = Eₚ[p]
    rel = iszero(ref) ? 0.0 : abs(est - ref) / ref
    @info "p=$p : ∫βₚ dτ ≈ $(round(est; digits=4)) vs Eₚ = $(round(ref; digits=4)) | rel.err = $(round(rel; digits=3))"
end
