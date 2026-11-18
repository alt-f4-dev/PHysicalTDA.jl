using Random, Statistics

# A steep band ω = 1 + 6|sin πq| and a flat band at 9.5, with Gaussian resolution (Γq, Γω) and k pixels per
# resolution element. The noise per bin grows as k, so each resolution element holds the same counts.
function phys_bands(k; Γq=0.04, Γω=0.5, seed=1, noise=0.04)
    rng = MersenneTwister(seed)
    qs = collect(range(0, 1; step=Γq / k))
    ωs = collect(range(0, 12; step=Γω / k))
    σq, σω = Γq / 2.355, Γω / 2.355
    I = [exp(-(ω - (1 + 6abs(sinpi(q))))^2 / (2σω^2)) + 0.6exp(-(ω - 9.5)^2 / (2σω^2)) for q in qs, ω in ωs]
    I .+= noise * k .* randn(rng, size(I))
    return I, qs, ωs
end

@testset "Segmentation in resolution elements" begin
    @testset "same parents on 1×, 2× and 4× grids" begin
        for k in (1, 2, 4)
            I, qs, ωs = phys_bands(k)
            r = adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_q=0.04, resolution_ω=0.5,
                                          physical_units=true, max_slope=20.0)
            @test r.qfirst.nparents == 2
            @test r.component_counts == [1, 1]
            @test r.radius_ω == k && r.radius_q == k
            @test 1 <= r.ridge_long_beta0 <= r.ridge_beta0
        end
    end

    @testset "slope-aware parent count" begin
        I, qs, ωs = phys_bands(2)
        r = adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_q=0.04, resolution_ω=0.5, max_slope=20.0)
        @test r.qfirst.nparents == 2
        @test r.qfirst.track_lengths[1:2] == [length(qs), length(qs)]
        @test r.max_slope == 20.0 && !r.physical_units
        # A short feature crossing a long track is linked one to one, so it stays a separate short track.
        peaks = [[10], [12], [14, 30], [16, 31], [18, 32], [20], [22], [24]]
        keep, lengths = PHysicalTDA.long_track_peaks(peaks, 3, 0.5)
        @test lengths == [8, 3]
        @test keep == [[10], [12], [14], [16], [18], [20], [22], [24]]
        keep_all, _ = PHysicalTDA.long_track_peaks(peaks, 3, 0.0)
        @test keep_all == peaks
        @test PHysicalTDA.long_track_peaks([[5], [6], [5]], 1, 1.0; periodic_q=true)[2] == [3]
        @test_throws ArgumentError PHysicalTDA.qfirst_parent_domains(I, qs, ωs; max_slope=-1.0, verbose=false)
    end

    @testset "slope-aware component labelling" begin
        mask = falses(4, 20)
        for (iq, iω) in enumerate((2, 7, 12, 17))
            mask[iq, iω:iω+1] .= true
        end
        @test PHysicalTDA.label_support_components(mask)[2] == 4
        @test PHysicalTDA.label_support_components(mask; max_jump_ω=5)[2] == 1
        @test PHysicalTDA.slope_neighbors(:full, 1) === PHysicalTDA.FULL_PIXEL_NEIGHBORS
        @test length(PHysicalTDA.slope_neighbors(:edge, 3)) == 4 + 2 * 7 - 2
    end

    @testset "prominence noise lag" begin
        s = Float64.(sin.(range(0, 2π; length=200)) .+ 0.01 .* randn(MersenneTwister(3), 200))
        PHysicalTDA.normalize01!(s)
        @test PHysicalTDA.adaptive_prominence(s) == PHysicalTDA.adaptive_prominence(s; lag=1)
        @test PHysicalTDA.adaptive_prominence(s; lag=500) >= 0.025
    end

    @testset "resolution report" begin
        I, qs, ωs = phys_bands(2)
        est = adaptive_autosplitbands2d(I, qs, ωs; verbose=false)
        rep = resolution_report(est)
        @test rep.source_ω === :estimated && rep.source_q === :default
        @test rep.Γω == rep.Γω_apparent && isnan(rep.ω_evaluated)
        sup = resolution_report(adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_q=0.04, resolution_ω=0.5))
        @test sup.Γω == 0.5 && sup.Γq == 0.04 && sup.source_ω === :supplied && sup.source_q === :supplied
        @test sup.Γω_apparent == rep.Γω_apparent
        fn = resolution_report(adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_ω=ω -> 0.1 + 0.05ω))
        @test fn.source_ω === :function
        @test fn.Γω ≈ 0.1 + 0.05fn.ω_evaluated
        @test_throws ArgumentError adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_ω=ω -> -1.0)
        @test_throws ArgumentError adaptive_autosplitbands2d(I, qs, ωs; verbose=false, resolution_estimate=:mean)
        @test_throws ArgumentError resolution_report((; radius_q=1))
        # With momentum resolution, a steep band broadens in energy; its flat columns are narrower.
        σq, σω = 0.04 / 2.355, 0.5 / 2.355
        δ = range(-2σq, 2σq; length=9)
        w = exp.(-δ .^ 2 ./ (2σq^2))
        qf, ωf = collect(0:0.01:1), collect(0:0.125:8)
        clean = [sum(w[i] * exp(-(ω - (1 + 6abs(sinpi(q + δ[i]))))^2 / (2σω^2)) for i in eachindex(δ)) for q in qf, ω in ωf]
        @test PHysicalTDA.estimate_flat_energy_resolution(clean, ωf) < PHysicalTDA.estimate_energy_resolution(clean, ωf)
    end

    @testset "energy-edge maxima in the ridge graph" begin
        qs = collect(range(0, 1; length=41))
        ωs = collect(range(0.45, 8.0; length=40))
        truth = ωs[1] .+ 12 .* abs.(qs .- 0.5)
        I = [exp(-(ω - ω0)^2 / (2 * 0.15^2)) for ω0 in truth, ω in ωs]
        t = build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=8)
        @test any(n -> n.iω == 1, t.nodes)
        @test t.beta0 == 1
        @test build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=8, include_energy_edges=false).beta0 == 2
        # Isolated maxima in the last energy bin are not ridge nodes.
        J = copy(I)
        J[10, end] = J[30, end] = 0.8
        tj = build_ridge_topology(J, qs, ωs; smooth_radius_ω=1, max_jump_bins=8)
        @test tj.diagnostics[:isolated_edge_nodes_dropped] == 2
        @test tj.beta0 == 1 && !any(n -> n.iω == length(ωs), tj.nodes)
        @test tj.diagnostics[:component_rows] == [length(qs)]
        r = adaptive_autosplitbands2d(I, qs, ωs; verbose=false)
        @test any(n -> n.iω == 1, r.ridge_topology.nodes)
    end
end
