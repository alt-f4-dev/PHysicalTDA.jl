using Random

et_gauss(x, Γ) = (s = Γ / (2sqrt(2log(2))); exp(-x^2 / (2s^2)) / (sqrt(2π) * s))

function et_map(f, qs, ωs; noise=0.03, seed=1)
    I0 = [f(q, ω) for q in qs, ω in ωs]
    σ = noise .* I0 .+ 0.01maximum(I0)
    return I0 .+ σ .* randn(MersenneTwister(seed), size(I0)), σ
end

@testset "Edge-pinned child tracks" begin
    qs = collect(range(0, 1; length=101)); Γ = 0.25
    params(ωs) = RidgeRecoveryParams(ωs, qs; max_slope=40.0, grid=:unrebinned, min_significance=5.0)
    # an elastic tail entering the bottom of the window next to a dispersive mode
    ωs = collect(0.3:0.02:6.0)
    I, σ = et_map((q, ω) -> 300et_gauss(ω, Γ) * exp(-(q - 0.5)^2 / (2 * 0.03^2)) + et_gauss(ω - (2 + 3q), Γ), qs, ωs)
    r = recover_child_ridges(I, qs, ωs; params=params(ωs), valid=trues(size(I)), σ)
    p = prune_edge_tracks(r; σ, span=Γ / 2)
    @test sum(r.global_children) == sum(p.result.global_children) + 1
    @test length(p.dropped) == 1 && p.dropped[1].convex > p.dropped[1].concave
    @test all(abs.(p.dropped[1].q .- 0.5) .< 0.1) && isempty(p.at_edge)
    kept = [T[b, :] for T in p.result.tracks for b in axes(T, 1)]
    @test any(t -> count(i -> isfinite(t[i]) && abs(t[i] - (2 + 3qs[i])) < 0.1, eachindex(qs)) > 90, kept)
    viakw = recover_child_ridges(I, qs, ωs; params=params(ωs), valid=trues(size(I)), σ, edge_tracks=:drop, edge_span=Γ / 2)
    @test viakw.global_children == p.result.global_children
    @test prune_edge_tracks(r; σ, span=Γ / 2, action=:flag).result.global_children == r.global_children
    @test_throws ArgumentError recover_child_ridges(I, qs, ωs; params=params(ωs), valid=trues(size(I)), σ, edge_tracks=:drop)
    @test_throws ArgumentError prune_edge_tracks(r; span=Γ / 2, action=:delete)
    # a flat mode whose maximum lies on the window edge is located there, not dropped
    ωs = collect(1.0:0.02:6.0)
    I, σ = et_map((q, ω) -> et_gauss(ω - 1.0, Γ) + et_gauss(ω - (2 + 3q), Γ), qs, ωs)
    p = prune_edge_tracks(recover_child_ridges(I, qs, ωs; params=params(ωs), valid=trues(size(I)), σ); σ, span=Γ / 2)
    @test isempty(p.dropped) && !isempty(p.at_edge)
    # a mode leaving through the top edge is kept and its edge columns are reported
    ωs = collect(0.5:0.02:6.0)
    I, σ = et_map((q, ω) -> et_gauss(ω - (2 + 7q), Γ) + et_gauss(ω - 1.5, Γ), qs, ωs)
    p = prune_edge_tracks(recover_child_ridges(I, qs, ωs; params=params(ωs), valid=trues(size(I)), σ); σ, span=Γ / 2)
    @test isempty(p.dropped) && !isempty(p.truncated)
end

@testset "Significance candidates and touching crossings" begin
    qs = collect(range(0, 1; length=81)); ωs = collect(0.0:0.02:8.0); Γ = 0.25
    # a branch 60 times weaker than its neighbour falls below the relative floor but not below 5σ
    I, _ = et_map((q, ω) -> 60et_gauss(ω - (1 + 2q), Γ) + et_gauss(ω - (6 - q), Γ), qs, ωs; noise=0.0)
    σ = fill(1e-3maximum(I), size(I))
    base = RidgeRecoveryParams(ωs, qs; max_slope=20.0, grid=:unrebinned, relative_score_floor=0.05)
    sig = RidgeRecoveryParams(ωs, qs; max_slope=20.0, grid=:unrebinned, relative_score_floor=0.05, min_significance=5.0)
    nweak(c) = count(n -> abs(n.ω - (6 - n.q)) < 0.1, c.nodes)
    @test nweak(find_ridge_candidates(I, qs, ωs; params=base, valid=trues(size(I)), σ)) == 0
    @test nweak(find_ridge_candidates(I, qs, ωs; params=sig, valid=trues(size(I)), σ)) > 70
    @test nweak(find_ridge_candidates(I, qs, ωs; params=sig, valid=trues(size(I)))) == 0   # needs σ
    # two tracks exchange order through a stretch where they touch (|Δ| ≤ tol): one crossing, at its middle
    q = Float32.(range(0, 1; length=11))
    T = Float32[1 2 3 4 4 4 4 5 6 7 8; 8 7 6 4.05 4 4 4.02 3 2 1 0]
    found = PHysicalTDA._touching_crossings!(NTuple{2,Float32}[], [T], q, 0.1)
    @test length(found) == 1 && isapprox(found[1][1], (q[3] + q[8]) / 2; atol=1e-6)
    @test isempty(PHysicalTDA._touching_crossings!(NTuple{2,Float32}[], [T[[1, 1], :]], q, 0.1))
end

@testset "Energy-limited exclusion sites" begin
    hs = collect(0.0:0.05:1.0); ω = collect(0.0:0.25:2.0); n = length(hs)
    h = rg_histogram((hs, hs, ω), ones(n, n, length(ω)), ones(n, n, length(ω)))
    full = exclusion_mask(h; sites=[((0.5, 0.5), 0.06)])
    low = exclusion_mask(h; sites=[((0.5, 0.5), 0.06, (-Inf, 0.6))])
    i = findfirst(==(0.5), hs)
    @test all(full[i, i, :]) && low[i, i, :] == (ω .<= 0.6)
    @test !any(low[1, 1, :]) && count(low) == count(ω .<= 0.6) * count(full[:, :, 1])
    @test exclusion_mask(h; sites=[((0.5, 0.5), 0.06, (-Inf, Inf))]) == full
    @test_throws ArgumentError exclusion_mask(h; sites=[((0.5, 0.5), 0.06, (1.0, 0.0))])
end

@testset "Elastic-line subtraction" begin
    ω = collect(-1.0:0.02:4.0); Γ = 0.25; q = range(0, 1; length=30)
    amp = [200exp(-(x - 0.5)^2 / (2 * 0.1^2)) + 1 for x in q]
    inel = [et_gauss(w - 2.0, Γ) for x in q, w in ω]
    I0 = [a * et_gauss(w - 0.01, Γ) for a in amp, w in ω] .+ inel .+ 0.02
    σ = 0.01 .* I0 .+ 1e-3
    I = I0 .+ σ .* randn(MersenneTwister(3), size(I0))
    r = subtract_elastic(I, ω; σ)
    @test isapprox(r.fwhm, Γ; rtol=0.02) && isapprox(r.centre, 0.01; atol=0.01)
    @test isapprox(r.amplitude ./ (amp .* et_gauss(0.0, Γ)), ones(length(q)); rtol=0.03)
    edge = findfirst(>=(0.4), ω)
    tail = I0[:, edge] .- inel[:, edge] .- 0.02
    resid = abs.(r.I[:, edge] .- inel[:, edge] .- 0.02)
    @test maximum(resid) < 0.1maximum(tail) && all(resid .< 5 .* r.σ[:, edge])
    peak = findfirst(>=(2.0), ω)
    @test isapprox(r.I[:, peak], I[:, peak] .- r.model[:, peak]; rtol=1e-12) && all(r.σ .>= σ)
    @test isapprox(r.I[:, peak], inel[:, peak] .+ 0.02; rtol=0.05)
    shape = x -> exp(-4log(2) * (x / Γ)^2)
    s = subtract_elastic(I, ω; σ, lineshape=shape, centre=0.01, window=Γ)
    @test isapprox(s.I, subtract_elastic(I, ω; σ, fwhm=Γ, centre=0.01).I; rtol=1e-10) && isnan(s.fwhm)
    valid = trues(size(I)); valid[1, :] .= false
    @test subtract_elastic(I, ω; σ, valid).I[1, :] == I[1, :]
    @test_throws ArgumentError subtract_elastic(I, ω; lineshape=shape)
    @test_throws DimensionMismatch subtract_elastic(I, ω[2:end])
end
