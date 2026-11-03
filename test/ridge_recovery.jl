using Test
using PHysicalTDA

@testset "Child-ridge recovery" begin
    Lor(ω, ω0, Γ) = 1 / (1 + ((ω - ω0) / Γ)^2)
    qs = collect(range(0.0, 1.0; length=101))
    ωs = collect(range(0.0, 10.0; length=201))
    Δω = ωs[2] - ωs[1]

    @testset "parameters" begin
        p = RidgeRecoveryParams()
        @test (p.support_fraction, p.support_max_half_width, p.support_touch) == (0.95f0, 0.22f0, 0.020f0)
        @test (p.max_slope, p.prediction_scale, p.min_child_q_fraction) == (16.0f0, 0.05f0, 0.025f0)
        c = RidgeRecoveryParams(0.0:0.4:8.0; max_slope=48)
        @test c.support_fraction == 0.5f0 && c.max_slope == 48.0f0
        @test c.support_max_half_width ≈ 0.8f0 && c.support_touch ≈ 0.4f0 && c.prediction_scale ≈ 0.4f0
        @test RidgeRecoveryParams(0.0:0.4:8.0; max_slope=48, support_fraction=0.7).support_fraction == 0.7f0
        @test RidgeRecoveryParams(; max_children=3).max_children == 3
    end

    @testset "single dispersive band" begin
        truth = 5 .+ 2 .* cospi.(2 .* qs)                  # max |dω/dq| ≈ 12.6 < max_slope
        I = [Lor(ω, ω0, 0.4) for ω0 in truth, ω in ωs]
        c = find_ridge_candidates(I, qs, ωs)
        @test c isa RidgeCandidates && length(c.byq) == length(qs)
        @test all(0 .<= c.A .<= 1)
        t = recover_child_ridges(c)
        @test length(t.parent_groups) == 1 && t.global_children == [1]
        @test isempty(t.split_nodes) && isempty(t.merge_nodes) && isempty(t.crossing_points)
        @test maximum(abs.(vec(t.tracks[1]) .- truth)) <= Δω
        d = child_ridge_dispersions(t)
        @test length(d.segments) == 1 && d.segments[1].q ≈ Float32.(qs)
        @test sprint(show, t) isa String
        @test isequal(recover_child_ridges(I, qs, ωs).tracks, t.tracks)
    end

    @testset "narrow steep band: defaults fragment, coarse-grid preset connects" begin
        # With the manuscript defaults a support is about one bin wide for a narrow peak, so a ridge
        # that moves by more than a bin between momentum rows breaks into several parents.
        truth = 5 .+ 2 .* cospi.(2 .* qs)
        I = [Lor(ω, ω0, 0.15) for ω0 in truth, ω in ωs]
        @test length(recover_child_ridges(I, qs, ωs).parent_groups) > 1
        t = recover_child_ridges(I, qs, ωs; params=RidgeRecoveryParams(ωs; max_slope=20))
        @test length(t.parent_groups) == 1 && t.global_children == [1]
        @test maximum(abs.(vec(t.tracks[1]) .- truth)) <= 2Δω
    end

    @testset "separate flat bands" begin
        I = [Lor(ω, 3.0, 0.15) + 0.6Lor(ω, 7.0, 0.15) for _ in qs, ω in ωs]
        t = recover_child_ridges(I, qs, ωs)
        @test length(t.parent_groups) == 2 && t.global_children == [1, 1]
        @test maximum(abs.(vec(t.tracks[1]) .- 3)) <= Δω && maximum(abs.(vec(t.tracks[2]) .- 7)) <= Δω
        @test !any(t.support_mask[:, abs.(ωs .- 5) .< 1])
    end

    @testset "crossing branches" begin
        up(q) = 5 + 3(q - 0.5); dn(q) = 5 - 3(q - 0.5)
        I = [Lor(ω, up(q), 0.4) + Lor(ω, dn(q), 0.4) for q in qs, ω in ωs]
        t = recover_child_ridges(I, qs, ωs)
        @test length(t.parent_groups) == 1 && t.global_children == [2]
        # Inside the unresolved region an exchange of order cannot be told apart from touching
        # branches, so only its location is checked.
        @test all(abs(q - 0.5) < 0.15 && abs(ω - 5) < 0.3 for (q, ω) in t.crossing_points)
        far = findall(q -> abs(q - 0.5) > 0.15, qs)
        lo, hi = vec(minimum(t.tracks[1]; dims=1)), vec(maximum(t.tracks[1]; dims=1))
        @test maximum(abs.(lo[far] .- min.(up.(qs[far]), dn.(qs[far])))) <= Δω
        @test maximum(abs.(hi[far] .- max.(up.(qs[far]), dn.(qs[far])))) <= Δω
    end

    @testset "argument checks" begin
        @test_throws DimensionMismatch find_ridge_candidates(zeros(5, 5), qs, ωs)
    end
end
