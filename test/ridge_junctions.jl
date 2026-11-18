using Random

const RJ_Γq, RJ_Γω = 0.02, 0.25
rj_gauss(x, Γ) = (s = Γ / (2sqrt(2log(2))); exp(-x^2 / (2s^2)) / (sqrt(2π) * s))
rj_x(q) = q - 0.5
rj_avoided(v1, v2, Δ) = (q -> 4 + (v1 + v2) / 2 * rj_x(q) - sqrt(((v1 - v2) * rj_x(q) / 2)^2 + (Δ / 2)^2),
                         q -> 4 + (v1 + v2) / 2 * rj_x(q) + sqrt(((v1 - v2) * rj_x(q) / 2)^2 + (Δ / 2)^2))

# Two branches with momentum resolution (Gaussian of FWHM Γq along q), 3 % noise and a 1 % floor.
function rj_map(E, W; seed=1, qs=collect(range(0, 1; length=121)), ωs=collect(range(0, 8; length=321)))
    δs = range(-RJ_Γq, RJ_Γq; length=7)
    wδ = exp.(-δs .^ 2 ./ (2 * (RJ_Γq / (2sqrt(2log(2))))^2)); wδ ./= sum(wδ)
    I = [sum(wδ[k] * sum(W[b](q + δs[k]) * rj_gauss(ω - E[b](q + δs[k]), RJ_Γω) for b in 1:2) for k in eachindex(δs))
         for q in qs, ω in ωs]
    σ = sqrt.((0.03 .* I) .^ 2 .+ (0.01maximum(I))^2)
    return I .+ σ .* randn(MersenneTwister(seed), size(I)), σ, qs, ωs
end

# Seeds from the true branches on resolved columns either side of [qa, qb].
function rj_seeds(E, qs, ωs, I, qa, qb)
    left = findall(q -> qa - 0.07 <= q < qa, qs); right = findall(q -> qb < q <= qb + 0.07, qs)
    seg(cols, f) = (q=qs[cols], ω=f.(qs[cols]), I=[I[i, argmin(abs.(ωs .- f(qs[i])))] for i in cols])
    lo(q) = min(E[1](q), E[2](q)); hi(q) = max(E[1](q), E[2](q))
    return (seg(left, lo), seg(left, hi), seg(right, lo), seg(right, hi)), first(left):last(right)
end

const RJ_X = ((q -> 4 + 5rj_x(q), q -> 4 - 5rj_x(q)), (q -> 1.0, q -> 1.0))
const RJ_SWAP = ((q -> 4 + 5rj_x(q), q -> 4 - 5rj_x(q)), (q -> 1 + 0.8tanh(rj_x(q) / 0.05), q -> 1 - 0.8tanh(rj_x(q) / 0.05)))
const RJ_GAP6 = (rj_avoided(5, -5, 0.6), (q -> 1.0, q -> 1.0))
const RJ_GAP1 = (rj_avoided(5, -5, 0.1), (q -> 1.0, q -> 1.0))

@testset "Two-level junction fit" begin
    for (case, E, W, expect) in (("crossing", RJ_X..., :cross), ("weights swap", RJ_SWAP..., :cross),
                                 ("resolved gap", RJ_GAP6..., :avoid), ("unresolved gap", RJ_GAP1..., :avoid))
        I, σ, qs, ωs = rj_map(E, W)
        seeds, cols = rj_seeds(E, qs, ωs, I, 0.47, 0.53)
        f = fit_junction(I, σ, qs, ωs, cols, (2.0, 6.0), seeds; resolution=(RJ_Γq, RJ_Γω))
        @test f.verdict === expect
        @test f.χ2ν < 1.5
        @test isapprox(f.Γq, RJ_Γq; rtol=0.25)
        e1, e2, g = diabatic_energies(f, 0.5)
        @test abs((e1 + e2) / 2 - 4) < 0.05
        case == "resolved gap" && @test isapprox(f.gap, 0.6; rtol=0.1)
        case == "unresolved gap" && @test 0.05 < f.gap < 0.2
        expect === :cross && @test f.gap < 0.1RJ_Γω
    end
    I, σ, qs, ωs = rj_map(RJ_X...)
    seeds, cols = rj_seeds(RJ_X[1], qs, ωs, I, 0.47, 0.53)
    f = fit_junction(I, σ, qs, ωs, cols, (2.0, 6.0), seeds; resolution=(RJ_Γq, RJ_Γω))
    M = junction_model(f, qs[cols], ωs)
    @test size(M) == (length(cols), length(ωs))
    @test sum(abs2, (M .- I[cols, :]) ./ σ[cols, :]) / length(M) < 1.5
    @test occursin("cross", sprint(show, f))
    @test_throws ArgumentError fit_junction(I, σ, qs, ωs, 1:2, (2.0, 6.0), seeds; resolution=(RJ_Γq, RJ_Γω))
    @test_throws ArgumentError fit_junction(I, σ, qs, ωs, cols, (2.0, 6.0), seeds[1:3]; resolution=(RJ_Γq, RJ_Γω))
    @test_throws DimensionMismatch fit_junction(I, σ[:, 2:end], qs, ωs, cols, (2.0, 6.0), seeds; resolution=(RJ_Γq, RJ_Γω))
end

@testset "Branch identity through junctions" begin
    params(qs, ωs) = RidgeRecoveryParams(ωs, qs; max_slope=40.0, grid=:unrebinned, touch_crossings=true, min_significance=5.0)
    follows(r, f) = maximum(count(i -> isfinite(r.ω[b, i]) && abs(r.ω[b, i] - f(r.q[i])) < 0.15, eachindex(r.q))
                            for b in axes(r.ω, 1)) / length(r.q)
    for (E, W, verdict) in ((RJ_X..., :cross), (RJ_SWAP..., :cross), (RJ_GAP1..., :avoid))
        I, σ, qs, ωs = rj_map(E, W)
        ch = recover_child_ridges(I, qs, ωs; params=params(qs, ωs), valid=trues(size(I)), σ)
        r = resolve_branches(ch, σ; resolution=(RJ_Γq, RJ_Γω))
        @test size(r.ω, 1) == 2 && length(r.junctions) == 1
        @test r.junctions[1].verdict === verdict
        @test all(abs.(r.unresolved[1] .- 0.5) .< 0.05)
        @test follows(r, E[1]) > 0.9 && follows(r, E[2]) > 0.9   # the 0.2-weight side of the swap drops a few columns
        @test any(r.from_model)
    end
    # labelling: energy order (adiabatic) never exchanges; diabatic follows the uncoupled modes through a gap
    I, σ, qs, ωs = rj_map(RJ_X...)
    ch = recover_child_ridges(I, qs, ωs; params=params(qs, ωs), valid=trues(size(I)), σ)
    ad = resolve_branches(ch, σ; resolution=(RJ_Γq, RJ_Γω), labelling=:adiabatic)
    @test all(i -> !(isfinite(ad.ω[1, i]) && isfinite(ad.ω[2, i])) ||
                   sign(ad.ω[1, i] - ad.ω[2, i]) == sign(ad.ω[1, 1] - ad.ω[2, 1]), eachindex(qs))
    I, σ, qs, ωs = rj_map(RJ_GAP1...)
    ch = recover_child_ridges(I, qs, ωs; params=params(qs, ωs), valid=trues(size(I)), σ)
    di = resolve_branches(ch, σ; resolution=(RJ_Γq, RJ_Γω), labelling=:diabatic)
    @test follows(di, q -> 4 + 5rj_x(q)) > 0.9 && follows(di, q -> 4 - 5rj_x(q)) > 0.9
    # trajectories and supports
    t = branch_trajectories(I, di, ωs; halfwidth=RJ_Γω / 2)
    @test size(t.ω) == size(di.ω) && length(t.supports) == 2
    @test t.imputed[di.from_model] == trues(count(di.from_model))
    @test any(t.supports[1] .& t.supports[2])          # the supports overlap at the junction
    @test occursin("junction", sprint(show, di))
    @test_throws ArgumentError resolve_branches(ch, σ; resolution=(RJ_Γq, RJ_Γω), labelling=:energy)
    @test_throws DimensionMismatch resolve_branches(ch, σ[2:end, :]; resolution=(RJ_Γq, RJ_Γω))
    @test_throws ArgumentError branch_supports(di, ωs; halfwidth=0)
end
