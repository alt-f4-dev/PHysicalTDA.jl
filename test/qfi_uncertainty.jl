@testset "QFI uncertainties" begin
    ωs = collect(0.25:0.25:6.0)
    T = 5.0
    S = [exp(-(ω - (1 + q))^2 / 0.3) for q in 0:0.5:2, ω in ωs]
    σ = fill(0.01, size(S))
    σb = fill(0.004, size(S))
    K = [qfi_kernel(ω, T) for ω in ωs]
    Δ = energy_bin_widths(ωs)

    r = qfi(S, ωs, T; σ, σ_sys=σb)
    @test r.FQ ≈ qfi(S, ωs, T; σ).FQ
    @test r.σsys ≈ fill(sum(K .* Δ .* 0.004), size(S, 1))
    @test r.σFQ ≈ fill(sqrt(sum((K .* Δ .* 0.01) .^ 2)), size(S, 1))
    @test qfi(S, ωs, T).σsys == zeros(size(S, 1))
    @test QFIResult([1.0], [0.1], [3], 5.0).σsys == [0.0]

    @test qfi_uncertainty(r) ≈ sqrt.(r.σFQ .^ 2 .+ r.σsys .^ 2)
    @test qfi_uncertainty(r; systematic=false) ≈ r.σFQ
    @test qfi_uncertainty(r; scale_σ=0.05) ≈ sqrt.(r.σFQ .^ 2 .+ r.σsys .^ 2 .+ (0.05 .* r.FQ) .^ 2)
    @test_throws ArgumentError qfi_uncertainty(r; scale_σ=-0.1)
    B = qfi_bound(; S=0.5)
    @test nqfi_lower_bound(qfi(S, ωs, T; σ); S=0.5) ≈ (r.FQ .- 2 .* r.σFQ) ./ B
    @test nqfi_lower_bound(r; S=0.5, scale_σ=0.05, nσ=1) ≈ (r.FQ .- qfi_uncertainty(r; scale_σ=0.05)) ./ B
    w0 = qfi_witness(r; S=0.5, systematic=false)
    w1 = qfi_witness(r; S=0.5, scale_σ=0.2)
    @test all(w1.nqfi_lb .< w0.nqfi_lb)

    sens = qfi_sensitivity(S, ωs, T; σ, σ_sys=σb, lower_limits=[0.0, 1.0, 2.0], S=0.5, scale_σ=0.05)
    @test size(sens.FQ) == (3, size(S, 1))
    @test sens.FQ[1, :] ≈ r.FQ
    @test sens.FQ[2, :] ≈ qfi(S, ωs, T; ωrange=(1.0, Inf)).FQ
    @test all(diff(sens.FQ; dims=1) .<= 0)
    @test sens.nqfi_lb[2, :] ≈ qfi_witness(qfi(S, ωs, T; σ, σ_sys=σb, ωrange=(1.0, Inf)); S=0.5, scale_σ=0.05).nqfi_lb
    @test sens.certified == (sens.nqfi_lb .> 1)
    @test_throws ArgumentError qfi_sensitivity(S, ωs, T; lower_limits=Float64[], S=0.5)

    low = falses(size(S)); low[:, ωs .< 2.5] .= true
    a = band_qfi(S, [low, .!low], ωs, T; σ, σ_sys=σb)
    b = band_qfi(S, [low .& (ωs' .> 0.5), .!low], ωs, T; σ, σ_sys=σb)
    @test a.total.σsys ≈ r.σsys
    @test vec(sum(a.σsys; dims=1)) ≈ r.σsys
    c = compare_band_qfi(a, a)
    @test all(iszero, c.Δξ) && all(c.consistent)
    d = compare_band_qfi(a, b)
    @test d.Δξ ≈ a.ξ .- b.ξ && d.σΔξ ≈ sqrt.(a.σξ .^ 2 .+ b.σξ .^ 2)
    @test any(!, compare_band_qfi(a, b; nσ=0).consistent)
    @test_throws DimensionMismatch compare_band_qfi(a, band_qfi(S, [low], ωs, T))
    dep = band_dependence(a; S=0.5, scale_σ=0.05)
    @test dep.nqfi_lb ≈ qfi_witness(a.total; S=0.5, scale_σ=0.05).nqfi_lb
end
