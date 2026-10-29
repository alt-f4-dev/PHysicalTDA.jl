using Test
using PHysicalTDA
using Random, Statistics
const PersistenceDiagram = PHysicalTDA.PersistenceDiagrams.PersistenceDiagram

const Lc(ω, ω0, Γ) = 1 / (1 + ((ω - ω0) / Γ)^2)

@testset "Opt-in conventions (defaults unchanged)" begin
    ωs = collect(range(0.0, 100.0; length=256))
    A2 = [Lc(ω, 30, 4) + 0.5Lc(ω, 70, 4) for ω in ωs]

    @testset "essential_cap" begin
        dflt = pd_array_intensities(A2; maxdim=0, superlevel=true, normalize=true)
        mx   = pd_array_intensities(A2; maxdim=0, superlevel=true, normalize=true, essential_cap=:max_death)
        thr  = pd_array_intensities(A2; maxdim=0, superlevel=true, normalize=true, essential_cap=:threshold)
        @test [(birth(x), death(x)) for x in dflt[1]] == [(birth(x), death(x)) for x in mx[1]]
        ess(pd) = only(x for x in pd[1] if birth(x) == minimum(birth.(pd[1])))
        @test death(ess(thr)) ≈ 1.0                 # complement of the normalized minimum
        @test death(ess(dflt)) < 1.0                # historical cap: largest finite death
        raw = pd_array_intensities(A2; maxdim=0, superlevel=true, essential_cap=:threshold)
        @test death(ess(raw)) ≈ -minimum(A2)        # background level in -A units
        @test_throws ArgumentError pd_array_intensities(A2; essential_cap=:bogus)
        @test length(pd_array_intensities(Float32.(A2); threshold=0.1f0)) == 2
    end

    @testset "betti_curve interval" begin
        D = PersistenceDiagram([(0.0, 1.0), (0.5, 2.0)])
        τ = [0.25, 0.75, 1.0, 1.5]
        @test betti_curve([D], τ; dims=0:0)[0] == [1, 2, 2, 1]
        @test betti_curve([D], τ; dims=0:0, interval=:halfopen)[0] == [1, 2, 1, 1]
        β, κ = betti_curvature([D], τ; dims=0:0, interval=:halfopen)
        @test β[0] == [1, 2, 1, 1]
        @test_throws ArgumentError betti_curve([D], τ; interval=:open)
    end

    @testset "naive peak_selection" begin
        rng = MersenneTwister(1)
        wrong = Dict(:height => 0, :persistence => 0); twoband = 0
        for t in 1:100
            spec = max.([Lc(ω, 30, 4) + 0.4Lc(ω, 70, 4) + 0.03randn(rng) for ω in ωs], 0)
            S = reshape(spec, 1, :)
            h = naivesplitbands(S, ωs)
            @test h == naivesplitbands(S, ωs; peak_selection=:height)
            p = naivesplitbands(S, ωs; peak_selection=:persistence)
            @test size(p.bandIq) == size(h.bandIq)              # count rule unchanged
            count(>(0), h.bandIq[:, 1]) == 2 || continue
            twoband += 1
            for (sel, out) in ((:height, h), (:persistence, p))
                hi = ωs[findlast(>(0), out.bandIω[1, :])]
                (40 <= hi <= 65) || (wrong[sel] += 1)
            end
        end
        @test twoband > 0
        @test wrong[:persistence] == 0
        @test wrong[:persistence] <= wrong[:height]

        qs = collect(range(0.0, 0.4; length=32))
        I = [Lc(ω, 20 + 60h, 3) + 0.5Lc(ω, 80 - 20h, 3) for h in qs, ω in ωs]
        @test autosplitbands(I, qs, ωs) == autosplitbands(I, qs, ωs; peak_selection=:height)
        @test autosplitbands(I, ωs) == autosplitbands(I, ωs; peak_selection=:height)
        @test autosplitbands(I, qs, ωs; peak_selection=:persistence).splitE.bandIq isa Matrix{Float64}
        @test_throws ArgumentError naivesplitbands(I, ωs; peak_selection=:bogus)
    end

    @testset "MAD is elementwise" begin
        x = [1.0, 2.0, 3.0, 4.0, 100.0]
        @test MAD(x) ≈ 1.4826 .* abs.(x .- median(x))
        @test median(MAD(x)) ≈ PHysicalTDA.robust_mad(x)
    end
end
