using Test
using PHysicalTDA

const Lor(ω, ω0, Γ) = 1 / (1 + ((ω - ω0) / Γ)^2)

@testset "Physical band segmentation" begin
    qs = collect(range(0.0, 1.0; length=96))
    ωs = collect(range(0.0, 100.0; length=192))
    Δq, Δω = qs[2] - qs[1], ωs[2] - ωs[1]

    two = [Lor(ω, 30, 2.5) + 0.6Lor(ω, 70, 3) for h in qs, ω in ωs]
    one = [Lor(ω, 20 + 50h, 2.5) for h in qs, ω in ωs]
    split(h) = 0.3 < h < 0.7 ? 12sin(π * (h - 0.3) / 0.4) : 0.0
    lens = [Lor(ω, 50 + split(h), 2.0) + Lor(ω, 50 - split(h), 2.0) for h in qs, ω in ωs]

    @testset "range and vector axes" begin
        rv = adaptive_autosplitbands2d(two, qs, ωs; verbose=false, periodic_q=false)
        rr = adaptive_autosplitbands2d(two, range(0.0, 1.0; length=96), range(0.0, 100.0; length=192);
                                       verbose=false, periodic_q=false)
        @test rr.qfirst.nparents == rv.qfirst.nparents
        @test [c.linear_indices for c in rr.components] == [c.linear_indices for c in rv.components]
        @test [c.weight for c in rr.components] ≈ [c.weight for c in rv.components]
        @test splitbands(two, range(0.0, 1.0; length=96), ωs; verbose=false).bandIq isa Matrix{Float64}
    end

    @testset "two flat bands" begin
        r = adaptive_autosplitbands2d(two, qs, ωs; verbose=false, periodic_q=false)
        @test r.qfirst.nparents == 2
        @test r.component_counts == [1, 1]
        @test r.ridge_beta1 == 0
        d = ridge_dispersions(two, qs, ωs; verbose=false, periodic_q=false)
        @test all(abs.(d.ω[1, :] .- 30) .<= Δω)
        @test all(abs.(d.ω[2, :] .- 70) .<= Δω)
        @test !any(d.imputed)
    end

    @testset "dispersion tracking" begin
        d = ridge_dispersions(one, qs, ωs; verbose=false, periodic_q=false)
        @test size(d.ω, 1) == 1
        @test maximum(abs.(d.ω[1, :] .- (20 .+ 50 .* qs))) <= Δω
    end

    @testset "split-merge lens" begin
        r = adaptive_autosplitbands2d(lens, qs, ωs; verbose=false, periodic_q=false)
        @test r.ridge_beta1 == 1
        @test r.ridge_beta1 == r.ridge_cycle_count
        @test r.ridge_status === :splitmerge_detected
        @test r.ridge_diagnostics[:raw_beta1] >= r.ridge_beta1
    end

    @testset "component and domain marginals" begin
        r = adaptive_autosplitbands2d(two, qs, ωs; verbose=false, periodic_q=false)
        for c in r.components
            @test sum(c.bandIq) * Δq ≈ c.weight * Δq * Δω
            @test sum(c.bandIω) * Δω ≈ c.weight * Δq * Δω
            @test c.npixels == length(c.linear_indices)
        end
        m = components_to_band_marginals(r)
        @test size(m.bandIq) == (length(r.components), length(qs))
        @test size(m.bandIω) == (length(r.components), length(ωs))
        # Supports are band cores: they keep a substantial but partial share.
        frac = sum(c.weight for c in r.components) / sum(two)
        @test 0.4 < frac < 1.0

        # Energy domains partition every q-row, so domain marginals conserve weight.
        dm = domain_marginals(two, r, qs, ωs)
        @test sum(dm.domainIq) * Δq ≈ sum(two) * Δq * Δω
        @test sum(dm.domainIω) * Δω ≈ sum(two) * Δq * Δω

        A = restrict_to_domain(two, r.qfirst.domains, 1)
        for iq in eachindex(qs)
            rng = r.qfirst.domains[iq][1]
            @test all(iszero, A[iq, setdiff(1:length(ωs), rng)])
            @test A[iq, rng] == two[iq, rng]
        end
    end

    @testset "periodic smoothing flag" begin
        I = zeros(16, 8); I[end, 4] = 1.0
        B = zeros(16, 8)
        PHysicalTDA.separable_smooth2d!(B, I; radius_q=1, radius_ω=0, periodic_q=true)
        @test B[1, 4] > 0
        PHysicalTDA.separable_smooth2d!(B, I; radius_q=1, radius_ω=0, periodic_q=false)
        @test B[1, 4] == 0
    end

    @testset "domain topology summary" begin
        r = adaptive_autosplitbands2d(two, qs, ωs; verbose=false, periodic_q=false)
        t = analyze_domain_topologies(two, qs, ωs, r)
        @test length(t.topologies) == r.qfirst.nparents
        @test t.summary.branch_counts == [1, 1]
        @test all(==(:single), t.summary.status)
    end
end
