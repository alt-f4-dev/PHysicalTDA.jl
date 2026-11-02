using Test
using PHysicalTDA

#------------------#
# Helper Functions #
#------------------#
finite_PD(PD) = [isempty(PD[d]) ? PD[d] :
                 PD[d][map(x -> isfinite(birth(x)) && isfinite(death(x)) && death(x) > birth(x), PD[d])]
                 for d in 1:length(PD)]

# Exact ∫β_p dτ by an event sweep over (half-open) bars.
function exact_integral_beta(PDfin, p::Int)
    d = p + 1
    (d > length(PDfin) || isempty(PDfin[d])) && return 0.0
    b   = Float64[birth(x) for x in PDfin[d]]
    dth = Float64[death(x) for x in PDfin[d]]
    t = vcat(b, dth)
    s = vcat(ones(Int, length(b)), -ones(Int, length(dth)))
    ord = sortperm(t); t = t[ord]; s = s[ord]
    area = 0.0; β = 0
    @inbounds for i in 1:length(t)-1
        β += s[i]
        area += β * (t[i+1] - t[i])
    end
    return area
end

function trapz(x::AbstractVector, y::AbstractVector)
    s = 0.0
    @inbounds for i in 1:length(x)-1
        s += (x[i+1] - x[i]) * (y[i+1] + y[i]) * 0.5
    end
    return s
end

# Trapezoid ∫β_p dτ on a uniform grid using the package Betti curve.
function grid_integral(PDfin, p::Int; n::Int=256)
    D = PDfin[p+1]
    τ = collect(range(minimum(birth.(D)), maximum(death.(D)); length=n))
    β = betti_curve(PDfin, τ; dims=p:p)[p]
    return trapz(τ, float.(β))
end

@testset "Superlevel filtration" begin
    A = zeros(64, 64); A[32, 32] = 10.0; A[20, 45] = 6.0
    PD_sub = pd_array_intensities(A; maxdim=1, superlevel=false, normalize=true)
    PD_sup = pd_array_intensities(A; maxdim=1, superlevel=true,  normalize=true)

    # Superlevel resolves the two isolated peaks as two H0 classes.
    @test length(PD_sup[1]) == 2
    @test sort(birth.(PD_sup[1])) ≈ [0.0, 0.4]                 # 1 - Â at the peaks
    # Sublevel sees only the connected zero background as one H0 class.
    @test length(PD_sub[1]) == 1

    # Total persistence equals the exact integral of β₀.
    PDu = finite_PD(PD_sup)
    E_sum = sum(death.(PDu[1]) .- birth.(PDu[1]))
    @test exact_integral_beta(PDu, 0) ≈ E_sum rtol=1e-10

    # Trapezoid estimate converges to the exact integral as the τ grid refines.
    exact = exact_integral_beta(PDu, 0)
    errs = [abs(grid_integral(PDu, 0; n=N) - exact) / exact for N in (64, 256, 1024, 4096)]
    @test errs[end] < errs[1]
    @test errs[end] < 1e-2

    # Superlevel PD of A matches sublevel PD of -A (normalized → identical Z).
    PD_neg = pd_array_intensities(-A; maxdim=1, superlevel=false, normalize=true)
    bars(D) = sort([(birth(x), death(x)) for x in D])
    @test first.(bars(PD_sup[1])) ≈ first.(bars(PD_neg[1]))
    @test last.(bars(PD_sup[1]))  ≈ last.(bars(PD_neg[1]))
end
