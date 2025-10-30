using .PHysicalTDA, Sunny, Ripserer, GLMakie, Statistics
import PersistenceDiagrams: birth, death


#--------------------------------#
# Load Crystal and LSWT Spectrum #
#--------------------------------#
cryst, sys = La2CuO4()
(result2D, results_paths, Es) = LSWT(cryst, sys)


#------------------------------#
# Construct 4D Intensity Field #
#------------------------------#
A4 = convert4D(result2D[1], result2D[2], Es)
hkl= collapse(A4; over=:w, op=sum)
hw = collapse(A4; over=(:k,:l), op=mean)

#------------------------------#
# Compute Persistence Diagrams #
#------------------------------#
PD, fig = pd_sunny_intensities(result2D[1]; maxdim=2, superlevel=true, normalize=true)
s1 = display(fig)
wait(s1)

#------------------------------#
# Topological Fingerprinting   #
#------------------------------#

# Persistence Entropy & Total Persistence
Sp, Ep = persistence_entropy(PD; dims=0:2, tol=0.0, N::Int=256)

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
betti, kappa = betti_curvature(PD, tau; dims=0:2, scheme=:forward)

#------------------------------#
# Validation Checks	       #
#------------------------------#
function trapz(x,y)
	s = 0.0
	@inbounds for i in 1:length(x)-1
		s += (x[i+1] - x[i]) * (y[i+1] + y[i]) * 0.5
	end
	return s
end

for p in 0:2
	betti_p = get(betti, p, zeros(length(tau)))
	Ep_estimate = trapz(tau, float.(betti_p))
	Ep_numerical = get(Ep, p, 0.0)
	rel_err = isapprox(Ep_numerical, 0) ? 0.0 : abs(Ep_estimate - Ep_numerical)/max(Ep_numerical, 1e-12)
	@info "p=$p ∫β_p dτ ≈ $(round(Ep_estimate,digits=4)) vs Ep = $(round(Ep_numerical,digits=4)) | rel.err = $(round(rel_err,digits=3))"
end
