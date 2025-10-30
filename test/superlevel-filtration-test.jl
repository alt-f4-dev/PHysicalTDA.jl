using Ripserer, PersistenceDiagrams, GLMakie
import Base: searchsortedlast, searchsortedfirst
#------------------#
# Helper Functions #
#------------------#
finite_PD(PD) = [isempty(PD[d]) ? PD[d] :
                 PD[d][map(x -> isfinite(birth(x)) && isfinite(death(x)) && death(x) > birth(x), PD[d])]
                 for d in 1:length(PD)]

#Exact (Finite) Integration
function exact_integral_beta(PDfin, p::Int)
    d = p + 1
    (d > length(PDfin) || isempty(PDfin[d])) && return 0.0
    b   = Float64[birth(x) for x in PDfin[d]]
    dth = Float64[death(x) for x in PDfin[d]]
    # events: +1 at births, -1 at deaths (half-open [b,d))
    t = vcat(b, dth)
    s = vcat(ones(Int, length(b)), -ones(Int, length(dth)))
    ord = sortperm(t); t = t[ord]; s = s[ord]
    area = 0.0; β = 0
    @inbounds for i in 1:length(t)-1
        β += s[i]
        dt = t[i+1] - t[i]
        area += β * dt
    end
    return area
end

#Trapezoid Integration
function trapz(x::AbstractVector, y::AbstractVector)
    s = 0.0
    @inbounds for i in 1:length(x)-1
        s += (x[i+1]-x[i]) * (y[i+1] + y[i]) * 0.5
    end
    return s
end

# show trapezoid bias shrinks with denser τ grids
function grid_bias(PDfin, p::Int; n::Int=256)
    b = birth.(PDfin[p+1]); d = death.(PDfin[p+1])
    τmin, τmax = minimum(b), maximum(d)
    τ = range(τmin, τmax; length=n)
    # β(τ) from events (right-continuous)
    bsort = sort(Float64.(b)); dsort = sort(Float64.(d))
    β = similar(τ, Int)
    @inbounds for j in eachindex(τ)
        t = τ[j]
        nb = searchsortedlast(bsort, t)
        nd = searchsortedfirst(dsort, t) - 1
        β[j] = nb - nd
    end
    trapz(collect(τ), float.(β))
end

#---------------------------------#
# Test: Superlevel Set Filtration #
#---------------------------------#
A = zeros(64,64); A[32,32] = 10.0; A[20,45] = 6.0
PD_sub, _ = pd_array_intensity(A; maxdim=1, superlevel=false, normalize=true)
PD_sup, _ = pd_array_intensity(A; maxdim=1, superlevel=true,  normalize=true)

# Expect: the dominant 0D pair in PD_sup has earlier birth (and typically larger lifetime)
b_sub = birth.(PD_sub[1]); d_sub = death.(PD_sub[1])
b_sup = birth.(PD_sup[1]); d_sup = death.(PD_sup[1])
@info "sublevel:  max Δ0 lifetime = $(maximum(d_sub .- b_sub))"
@info "superlevel: max Δ0 lifetime = $(maximum(d_sup .- b_sup))"

# finite-only views (aligns with persistence_entropy’s finite-lifetime use)
PDs = finite_PD(PD_sub)
PDu = finite_PD(PD_sup)

# exact integral vs lifetime sum for p=0
E0_sub_exact = exact_integral_beta(PDs, 0)
E0_sup_exact = exact_integral_beta(PDu, 0)
E0_sub_sum   = sum(death.(PDs[1]) .- birth.(PDs[1]))
E0_sup_sum   = sum(death.(PDu[1]) .- birth.(PDu[1]))
@info "SUB: exact ∫β0 dτ = $E0_sub_exact   vs   Σ(d-b) = $E0_sub_sum   rel.err = $(abs(E0_sub_exact-E0_sub_sum)/max(E0_sub_sum,1e-12))"
@info "SUP: exact ∫β0 dτ = $E0_sup_exact   vs   Σ(d-b) = $E0_sup_sum   rel.err = $(abs(E0_sup_exact-E0_sup_sum)/max(E0_sup_sum,1e-12))"

#Trapezoid Integration Bias Should Shrink
for N in (64, 256, 1024, 4096)
    approx = grid_bias(PDu, 0; n=N)
    exact  = E0_sup_exact
    rel    = abs(approx - exact)/max(exact,1e-12)
    @info "τ-grid N=$N: trapezoid ≈ $approx vs exact = $exact  rel.err = $rel"
end

# PD_super(A) ≈ PD_sub(-A)
PD_sub, _ = pd_array_intensity(A; superlevel=false, normalize=true)
PD_sup, _ = pd_array_intensity(A; superlevel=true,  normalize=true)
PD_sub_neg, _ = pd_array_intensity(-A; superlevel=false, normalize=true)

# finite-pair identity (event sweep)
PDf = finite_PD(PD_sup)
E_sum = sum([death(x) - birth(x) for x in PDf[1]])
E_exact = exact_integral_beta(PDf, 0)
@assert isapprox(E_exact, E_sum; rtol=1e-10)
