using PHysicalTDA, Sunny, Ripserer, GLMakie, Statistics
import PersistenceDiagrams: birth, death


#--------------------------------#
# Load Crystal and LSWT Spectrum #
#--------------------------------#
#Pre-built Crystal Structure
cryst, sys = La2CuO4() 

#q-paths (wavevector paths)
H = [[h,0,0] for h in -1:1]; H = Sunny.q_space_path(cryst, H, 500) #convert to
K = [[0,k,0] for k in -1:1]; K = Sunny.q_space_path(cryst, K, 500) #Sunny type
L = [[0,0,ℓ] for ℓ in -1:1]; L = Sunny.q_space_path(cryst, L, 500) #for LSWT

#Linear response theory calculation (per path)
qpaths = [H,K,L]; LCO = LSWT(cryst, sys, qpaths)

#Every q-path has the same energy axis
ωaxis = LCO.energies

#Obtaining Results for L-axis
path₃ = LCO.qpaths[3]
IntensityMatrix₃ = LCO.intensities[3]

#------------------------------#
# Construct 4D Intensity Field #
#------------------------------#

# 3D Grid of Q-points
Qs = [[h,k,ℓ] for h in -1:0.1:1, k in -1:0.1:1, ℓ in -1:0.1:1]

# Perform LSWT
swt = SpinWaveTheory(sys; measure=ssf_perp(sys))
result2D = intensities(swt, Qs[:]; energies, kernel=Gaussian(fwhm=35))

# Convert to 4D
result4D = convert4D(result2D, Qs, Es)

#Slice projections
hkl= collapse(result4D; over=:w, op=sum)
hw = collapse(result4D; over=(:k,:l), op=mean)

#------------------------------#
# Compute Persistence Diagrams #
#------------------------------#
PD = pd_sunny_intensities(result2D; maxdim=2, superlevel=true, normalize=true)

#------------------------------#
# Topological Fingerprinting   #
#------------------------------#

# Persistence Entropy & Total Persistence
Sₚ, Eₚ = persistence_entropy(PD; dims=0:2, tol=0.0, N::Int=256)

#Construct tau-grid spanning births..deaths in PD
function taugrid(PD)
	births = Float64[]; deaths = Float64[]; lifetimes = Float64[]
	for d in 1:length(PD)
		bd = PD[d]; isempty(bd) && continue
		b = filter(isfinite, Float64.(birth.(bd)))
		d = filter(isfinite, Float64.(death.(bd)))
		append!(births, b); append!(deaths, d)
		@inbounds for i in eachindex(b)
			bi = Float64(b[i])
			di = Float64(d[i])
			isfinite(bi) && isfinite(di) && push!(lifetimes, di-bi)
		end
	end
	@assert !isempty(births) || !isempty(deaths) "Empty PD: no features!"

	tmin = minimum(births); tmax = maximum(deaths)
	if !(isfinite(tmin) && isfinite(tmax)) || !(tmax > tmin)
		if !isempty(lifetimes)
			L = quantile(lifetimes, 0.9)
			tmin = isfinite(tmin) ? tmin : 0.0
			tmax = tmin + max(L, eps())
		else
			tmin, tmax = 0.0, 1.0
		end
	end
	delta = 1e-9 * max(1.0, abs(tmax-tmin))
	return range(tmin + delta, tmax - delta; length=N)
end
tau = taugrid(PD)

# Betti curves and Betti curvature on tau-grid
β, κ = betti_curvature(PD, tau; dims=0:2, scheme=:forward)

#-----------------------------------#
# Validation Check: εₚ = ∫βₚ(τ)dτ   #
#-----------------------------------#
function trapz(x,y)
	s = 0.0
	@inbounds for i in 1:length(x)-1
		s += (x[i+1] - x[i]) * (y[i+1] + y[i]) * 0.5
	end
	return s
end

for p in 0:2
	βₚ = get(β, p, zeros(length(tau)))
	Eₚ_estimate = trapz(tau, float.(betti_p))
	Eₚ_numerical = get(Eₚ, p, 0.0)
	rel_err = isapprox(Eₚ_numerical, 0) ? 0.0 : abs(Eₚ_estimate - Eₚ_numerical)/max(Eₚ_numerical, 1e-12)
	@info "p=$p : ∫βₚ dτ ≈ $(round(Eₚ_estimate,digits=4)) vs Ep = $(round(Eₚ_numerical,digits=4)) | rel.err = $(round(rel_err,digits=3))"
end
