using Test
using PHysicalTDA

@testset "Naive band marginal grid weighting" begin
    nq, nω = 24, 31
    qs = collect(range(0.0, 0.4; length=nq))
    ωs = collect(range(0.0, 100.0; length=nω))
    Δq, Δω = qs[2] - qs[1], ωs[2] - ωs[1]
    Iqω = [0.2 + (1.0 + 0.4sin(4π * (i - 1) / (nq - 1))) *
                  (exp(-((j - 9) / 2)^2) + 0.7exp(-((j - 24) / 3)^2))
           for i in 1:nq, j in 1:nω]

    legacy = autosplitbands(Iqω, qs, ωs)
    legacy_2arg = autosplitbands(Iqω, ωs)
    physical = autosplitbands(Iqω, qs, ωs; grid_weighting=:physical)
    @test sum(abs, legacy.splitE.bandIq) > 0
    @test sum(abs, legacy.splitQ.bandIω) > 0
    @test legacy_2arg.splitE == legacy.splitE
    @test size(physical.splitE.bandIq, 2) == nq
    @test size(physical.splitQ.bandIω, 2) == nω

    # Orientation-dependent masks remain unchanged; only the missing
    # integration factors are introduced in the corresponding arrays.
    @test physical.splitE.bandIq == legacy.splitE.bandIq
    @test physical.splitE.bandIω ≈ Δq .* legacy.splitE.bandIω
    @test physical.splitQ.bandIq ≈ Δω .* legacy.splitQ.bandIq
    @test physical.splitQ.bandIω == legacy.splitQ.bandIω

    for split in (physical.splitE, physical.splitQ)
        @test sum(split.bandIq) * Δq ≈ sum(split.bandIω) * Δω
    end

    direct = naivesplitbands(Iqω, ωs)
    weighted = naivesplitbands(Iqω, ωs; row_spacing=Δq)
    @test direct.bandIq == weighted.bandIq
    @test weighted.bandIω ≈ direct.bandIω .* Δq

    # Four complete-array normalizations are invariant to constant spacing.
    for orientation in (:splitE, :splitQ), marginal in (:bandIq, :bandIω)
        x = getproperty(getproperty(legacy, orientation), marginal)
        y = getproperty(getproperty(physical, orientation), marginal)
        sum(abs, x) > 0 || continue
        @test x ./ sum(abs, x) ≈ y ./ sum(abs, y)
    end

    @test_throws ArgumentError autosplitbands(Iqω, ωs; grid_weighting=:physical)
    @test_throws ArgumentError autosplitbands(Iqω, qs, ωs; grid_weighting=:unknown)
    @test_throws ArgumentError naivesplitbands(Iqω, ωs; row_spacing=-1.0)
    @test_throws ArgumentError autosplitbands(Iqω, reverse(qs), ωs; grid_weighting=:physical)
    @test_throws ArgumentError autosplitbands(Iqω, qs, reverse(ωs); grid_weighting=:physical)
    @test_throws ArgumentError autosplitbands(Iqω, qs .^ 2, ωs; grid_weighting=:physical)
    @test_throws ArgumentError autosplitbands(Iqω, qs, ωs .^ 2; grid_weighting=:physical)
    @test_throws DimensionMismatch autosplitbands(Iqω, qs[1:end-1], ωs; grid_weighting=:physical)
end
