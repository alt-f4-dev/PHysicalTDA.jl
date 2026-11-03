using Test, Random, Statistics
using PHysicalTDA

@testset "Resolution-matched rebinning" begin
    rng = MersenneTwister(7)
    qs = collect(range(0.0, 1.0; length=40))
    ωs = collect(range(-2.0, 10.0; length=121))
    I = 0.1 .+ rand(rng, 40, 121)
    σ = 0.05 .+ 0.01 .* rand(rng, 40, 121)
    I[3:5, 10:20] .= NaN                      # unmeasured patch
    σ[3:5, 10:20] .= NaN

    @testset "block means and factors" begin
        rb = rebin_to_resolution(I, qs, ωs; σ=σ, factors=(2, 4), average=:all)
        @test (rb.qfactor, rb.ωfactor) == (2, 4)
        @test size(rb.I) == (20, 30)                     # 121 = 30·4 + 1: last energy bin dropped
        @test rb.source_size == (40, 121)
        @test rb.fillvalue == minimum(filter(x -> isfinite(x) && x > 0, I))
        @test rb.nfilled == 33
        Ifill = replace(I, NaN => rb.fillvalue)
        @test rb.I[7, 9] ≈ mean(Ifill[13:14, 33:36])
        @test rb.q[7] ≈ mean(qs[13:14]) && rb.ω[9] ≈ mean(ωs[33:36])
        @test rb.σ[7, 9] ≈ sqrt(sum(abs2, σ[13:14, 33:36])) / 8
        @test sprint(show, rb) isa String
    end

    @testset ":valid excludes unmeasured bins" begin
        rb = rebin_to_resolution(I, qs, ωs; σ=σ, factors=(2, 4))
        @test rb.average === :valid
        blk = I[5:6, 9:12]                               # row 5 unmeasured for columns 10:12
        good = filter(isfinite, blk)
        @test rb.I[3, 3] ≈ mean(good)
        @test rb.σ[3, 3] ≈ sqrt(sum(abs2, filter(isfinite, σ[5:6, 9:12]))) / length(good)
        @test rb.valid[3, 3]
        ra = rebin_to_resolution(I, qs, ωs; σ=σ, factors=(2, 4), average=:all)
        @test ra.I[3, 3] < rb.I[3, 3]                    # the noise floor pulls the :all mean down
        # fully measured blocks agree
        @test rb.I[10, 20] ≈ ra.I[10, 20]
    end

    @testset "empty blocks" begin
        J = copy(I); J[1:2, 1:4] .= NaN
        rb = rebin_to_resolution(J, qs, ωs; σ=σ, factors=(2, 4), fill=:zero)
        @test rb.fillvalue == 0.0
        @test !rb.valid[1, 1] && rb.I[1, 1] == 0.0 && isnan(rb.σ[1, 1])
        @test rebin_to_resolution(J, qs, ωs; factors=(2, 4)).σ === nothing
    end

    @testset "factors from resolution" begin
        Δq, Δω = qs[2] - qs[1], ωs[2] - ωs[1]
        rb = rebin_to_resolution(I, qs, ωs; resolution_q=3Δq, resolution_ω=5Δω)
        @test (rb.qfactor, rb.ωfactor) == (3, 5)
        @test rb.resolution_q ≈ 3Δq && rb.resolution_ω ≈ 5Δω
        @test size(rb.I) == (13, 24)
        rr = rebin_to_resolution(I, range(0.0, 1.0; length=40), range(-2.0, 10.0; length=121);
                                 resolution_q=3Δq, resolution_ω=5Δω)
        @test isequal(rr.I, rb.I)
    end

    @testset "upsample_mask" begin
        rb = rebin_to_resolution(I, qs, ωs; factors=(3, 5))     # 40 = 13·3 + 1, 121 = 24·5 + 1
        M = falses(size(rb.I)); M[2, 3] = true; M[13, 24] = true
        F = upsample_mask(M, rb)
        @test size(F) == (40, 121)
        @test count(F) == 2 * 15
        @test all(F[4:6, 11:15]) && all(F[37:39, 116:120])
        @test !any(F[40, :]) && !any(F[:, 121])                 # dropped bins
        A = trues(size(rb.I)); B = .!A; B[1, 1] = true; A[1, 1] = false
        @test !any(upsample_mask(A, rb) .& upsample_mask(B, rb))
        @test_throws DimensionMismatch upsample_mask(falses(2, 2), rb)
    end

    @testset "remainder = :merge" begin
        rd = rebin_to_resolution(I, qs, ωs; σ=σ, factors=(3, 5))
        rm = rebin_to_resolution(I, qs, ωs; σ=σ, factors=(3, 5), remainder=:merge)
        @test rd.remainder === :drop && rm.remainder === :merge
        @test size(rm.I) == size(rd.I) == (13, 24)
        @test last(rd.qblocks) == 37:39 && last(rm.qblocks) == 37:40
        @test last(rd.ωblocks) == 116:120 && last(rm.ωblocks) == 116:121
        @test rm.I[1:end-1, 1:end-1] == rd.I[1:end-1, 1:end-1]       # only the last blocks change
        @test rm.I[end, end] ≈ mean(I[37:40, 116:121])
        @test rm.σ[end, end] ≈ sqrt(sum(abs2, σ[37:40, 116:121])) / 24
        @test rm.q[end] ≈ mean(qs[37:40]) && rm.ω[end] ≈ mean(ωs[116:121])
        @test all(upsample_mask(trues(size(rm.I)), rm))                # every bin belongs to a block
        @test count(upsample_mask(trues(size(rd.I)), rd)) == 39 * 120
        ra = rebin_to_resolution(I, qs, ωs; factors=(3, 5), average=:all, remainder=:merge)
        @test ra.I[end, end] ≈ mean(replace(I[37:40, 116:121], NaN => ra.fillvalue))
        # an exact multiple is unaffected
        @test isequal(rebin_to_resolution(I, qs, ωs; factors=(2, 11), remainder=:merge).I,
                      rebin_to_resolution(I, qs, ωs; factors=(2, 11)).I)
        @test occursin("remainder=:merge", sprint(show, rm))
    end

    @testset "argument checks" begin
        @test_throws ArgumentError rebin_to_resolution(I, qs, ωs; fill=:mean)
        @test_throws ArgumentError rebin_to_resolution(I, qs, ωs; average=:median)
        @test_throws ArgumentError rebin_to_resolution(I, qs, ωs; remainder=:pad)
        @test_throws DimensionMismatch rebin_to_resolution(I, qs[1:10], ωs)
        @test_throws DimensionMismatch rebin_to_resolution(I, qs, ωs; σ=σ[1:5, :])
        @test_throws ArgumentError rebin_to_resolution(I, qs, ωs; valid=trues(size(I)))
        @test_throws ArgumentError rebin_to_resolution(I, qs, ωs; factors=(0, 2))
    end

    @testset "segmentation on the rebinned grid" begin
        qf = collect(range(0.0, 1.0; length=96)); ωf = collect(range(0.0, 100.0; length=384))
        G(ω, ω0, s) = exp(-(ω - ω0)^2 / 2s^2)
        Iqω = [G(ω, 30, 2.0) + 0.6G(ω, 70, 2.5) for h in qf, ω in ωf]
        rb = rebin_to_resolution(Iqω, qf, ωf; factors=(2, 4))
        seg = adaptive_autosplitbands2d(rb.I, rb.q, rb.ω; verbose=false)
        coarse = band_supports(seg, size(rb.I)...)
        fine = band_supports(seg, rb)
        @test length(fine) == seg.qfirst.nparents == 2
        @test fine == [upsample_mask(M, rb) for M in coarse]
        @test all(size(M) == size(Iqω) for M in fine)
        @test !any(fine[1] .& fine[2])
        @test count(fine[1]) == 8 * count(coarse[1])
        # most of the weight lies inside the fine-grid supports (the segmentation leaves out the
        # weakest tails, about 6% here)
        br = band_qfi(Iqω, fine, ωf, 10.0; ωrange=(1.0, Inf))
        @test all(br.FQres .< 0.1 .* br.total.FQ)
        @test all(abs.(sum(br.ξ; dims=1) .+ br.FQres' ./ br.total.FQ' .- 1) .< 1e-10)
    end
end

@testset "Parent-band dispersions from a result" begin
    Lor(ω, ω0, Γ) = 1 / (1 + ((ω - ω0) / Γ)^2)
    qs = collect(range(0.0, 1.0; length=96))
    ωs = collect(range(0.0, 100.0; length=192))
    Δω = ωs[2] - ωs[1]
    truth = 20 .+ 50 .* qs
    one = [Lor(ω, 20 + 50h, 2.5) for h in qs, ω in ωs]
    r = adaptive_autosplitbands2d(one, qs, ωs; verbose=false)
    d0 = ridge_dispersions(one, qs, ωs; verbose=false)
    db = ridge_dispersions(r, qs, ωs)
    dw = ridge_dispersions(r, qs, ωs; centre=:weighted)
    @test db.ω == d0.ω && db.center_indices == d0.center_indices && db.imputed == d0.imputed
    @test maximum(abs.(dw.ω[1, :] .- truth)) <= Δω
    @test sum(abs, dw.ω[1, :] .- truth) < sum(abs, db.ω[1, :] .- truth)
    @test_throws ArgumentError ridge_dispersions(r, qs, ωs; centre=:mode)
    @test_throws DimensionMismatch ridge_dispersions(r, qs[1:10], ωs)
end
