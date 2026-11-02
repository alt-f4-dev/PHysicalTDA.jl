using Test, PHysicalTDA, LinearAlgebra, Random

# Two spin-1/2 sites: exact QFI from the definition (Eq. S19) for a Hermitian probe.
function _dimer(J, h)
    sx = [0 1; 1 0] / 2; sy = [0 -im; im 0] / 2; sz = [1 0; 0 -1] / 2
    id = Matrix{ComplexF64}(I, 2, 2)
    S1 = [kron(s, id) for s in (sx, sy, sz)]
    S2 = [kron(id, s) for s in (sx, sy, sz)]
    H = J * sum(S1[a] * S2[a] for a in 1:3) + h * (S1[3] + S2[3]) + 0.3h * S1[1]
    return Hermitian(H), S1, S2
end

function _qfi_definition(H, O, β)
    E, V = eigen(H)
    p = exp.(-β .* (E .- minimum(E))); p ./= sum(p)
    Om = V' * O * V
    F = 0.0
    for λ in eachindex(E), μ in eachindex(E)
        s = p[λ] + p[μ]
        s > 0 || continue
        F += 2 * (p[λ] - p[μ])^2 / s * abs2(Om[λ, μ])
    end
    return F
end

# S(ω) = Σ p_λ |⟨μ|O|λ⟩|² δ(ω − (E_μ − E_λ)), as (energies, weights) of positive-ω lines.
function _positive_lines(H, O, β)
    E, V = eigen(H)
    p = exp.(-β .* (E .- minimum(E))); p ./= sum(p)
    Om = V' * O * V
    ωs = Float64[]; w = Float64[]
    for λ in eachindex(E), μ in eachindex(E)
        ω = E[μ] - E[λ]
        ω > 1e-12 || continue
        push!(ωs, ω); push!(w, p[λ] * abs2(Om[μ, λ]))
    end
    return ωs, w
end

@testset "QFI" begin
    kB = KB_MEV_PER_K

    @testset "kernel" begin
        @test qfi_kernel(-1.0, 10.0) == 0.0
        @test qfi_kernel(0.0, 10.0) == 0.0
        @test qfi_kernel(1.0, 0.0) == 4.0
        @test qfi_kernel(50.0, 1.0) ≈ 4.0
        for (ω, T) in ((0.3, 5.0), (2.0, 30.0), (7.0, 300.0))
            x = ω / (kB * T)
            @test qfi_kernel(ω, T) ≈ 4 * tanh(x / 2) * (1 - exp(-x))
        end
        x = 1e-4; T = 10.0; ω = x * kB * T
        @test qfi_kernel(ω, T) ≈ 2x^2 - x^3 rtol = 1e-6   # K = 2x² − x³ + O(x⁴)
        Ks = [qfi_kernel(ω, 20.0) for ω in range(0.01, 20; length=200)]
        @test all(0 .<= Ks .<= 4) && issorted(Ks)
        @test_throws ArgumentError qfi_kernel(1.0, -1.0)
    end

    @testset "spectral QFI equals the definition (Eq. S19 vs S25)" begin
        J, h = 1.0, 0.4                             # meV
        H, S1, S2 = _dimer(J, h)
        for O in (S1[1] - S2[1], S1[3] + S2[3], S1[1] + 0.5S2[2] + 0.5S2[2]')
            O = (O + O') / 2
            for T in (0.0, 2.0, 8.0, 40.0)
                β = T == 0 ? Inf : 1 / (kB * T)
                Fdef = T == 0 ? begin
                    E, V = eigen(H); g = V[:, 1]
                    4 * real(g' * O * O * g - (g' * O * g)^2)
                end : _qfi_definition(H, O, β)
                ωl, wl = _positive_lines(H, O, T == 0 ? 1e12 : β)
                isempty(ωl) && continue
                ord = sortperm(ωl)
                Sω = reshape(wl[ord], 1, :)        # one momentum row; weights as unit-width bins
                r = qfi(Sω, ωl[ord] .+ (1:length(ωl)) .* 1e-9, T; Δω=1.0)
                @test r.FQ[1] ≈ Fdef rtol = 1e-6 atol = 1e-10
            end
        end
    end

    @testset "bin widths, masks and NaN handling" begin
        @test energy_bin_widths([0.0, 1.0, 2.0, 3.0]) == [1.0, 1.0, 1.0, 1.0]
        @test energy_bin_widths([0.0, 1.0, 3.0]) == [1.0, 1.5, 2.0]
        @test_throws ArgumentError energy_bin_widths([0.0, 0.0, 1.0])

        ωs = collect(-1.0:0.5:5.0)
        S = ones(2, length(ωs))
        r = qfi(S, ωs, 0.0)
        @test r.FQ ≈ fill(4 * 0.5 * count(>(0), ωs), 2)  # T = 0: 4 Σ S Δω over ω > 0
        @test r.nbins == fill(count(>(0), ωs), 2)
        r2 = qfi(S, ωs, 0.0; ωrange=(1.0, 3.0))
        @test r2.FQ ≈ fill(4 * 0.5 * count(ω -> 1 < ω <= 3, ωs), 2)
        S[1, end] = NaN
        r3 = qfi(S, ωs, 0.0)
        @test r3.FQ[1] ≈ r.FQ[1] - 2.0 && r3.nbins[1] == r.nbins[1] - 1
        D = qfi_density(S, ωs, 0.0)
        @test isnan(D[1, end]) && all(D[:, ωs .<= 0] .== 0)
    end

    @testset "uncertainty propagation" begin
        ωs = collect(range(0.5, 6; length=12)); T = 15.0
        S = rand(MersenneTwister(1), 3, 12); σ = 0.1 .* (1 .+ S)
        r = qfi(S, ωs, T; σ=σ)
        w = energy_bin_widths(ωs); K = qfi_kernel.(ωs, T)
        @test r.σFQ ≈ [sqrt(sum((K .* w .* σ[i, :]) .^ 2)) for i in 1:3]
        @test r.FQ ≈ [sum(K .* w .* S[i, :]) for i in 1:3]
    end

    @testset "normalization and witness" begin
        Sspin = 5 / 2
        r = QFIResult([12Sspin^2 * 0.5, 12Sspin^2 * 2.5], [0.0, 12Sspin^2 * 0.4], [1, 1], 5.0)
        @test nqfi(r; S=Sspin) ≈ [0.5, 2.5]
        @test nqfi(r; S=Sspin, ncomponents=1) ≈ 3 .* [0.5, 2.5]
        @test nqfi_lower_bound(r; S=Sspin, nσ=2) ≈ [0.5, 2.5 - 0.8]
        w = qfi_witness(r; S=Sspin, m=1)
        @test w.certified == [false, true] && w.m == 1
        @test qfi_witness(r; S=Sspin, m=2).certified == [false, false]
        @test_throws ArgumentError nqfi(r; S=-1)
    end

    @testset "band decomposition" begin
        rng = MersenneTwister(2)
        nq, nω = 6, 40
        ωs = collect(range(-1, 10; length=nω)); T = 20.0
        S = rand(rng, nq, nω); σ = 0.05 .* ones(nq, nω)
        M1 = falses(nq, nω); M1[:, 10:18] .= true
        M2 = falses(nq, nω); M2[:, 25:33] .= true
        br = band_qfi(S, [M1, M2], ωs, T; σ=σ)
        @test br.total.FQ ≈ vec(sum(br.FQ; dims=1)) .+ br.FQres
        @test br.total.σFQ .^ 2 ≈ vec(sum(br.σFQ .^ 2; dims=1)) .+ br.σres .^ 2
        @test all(0 .<= br.ξ .<= 1) && all(vec(sum(br.ξ; dims=1)) .<= 1 + 1e-12)
        @test br.FQ[1, :] ≈ qfi(S, ωs, T; mask=M1).FQ

        Mrest = .!(M1 .| M2)
        brfull = band_qfi(S, [M1, M2, Mrest], ωs, T)
        @test vec(sum(brfull.ξ; dims=1)) ≈ ones(nq)
        @test all(abs.(brfull.FQres) .< 1e-12)

        Mbad = copy(M1); Mbad[1, 25] = true
        @test_throws ArgumentError band_qfi(S, [Mbad, M2], ωs, T)

        dep = band_dependence(br; S=0.5, m=1, nσ=0)
        @test dep.nqfi_lb_without ≈ dep.nqfi_lb' .- br.FQ ./ 3
        @test dep.certified_without == (dep.nqfi_lb_without .> 1)
    end

    @testset "supports from a segmentation result" begin
        qs = range(0, 1; length=64); ωs = range(0.1, 10; length=128)
        Iqω = [exp(-(ω - 4 - 1.5sinpi(q))^2 / 0.2) + 0.6exp(-(ω - 8)^2 / 0.2) for q in qs, ω in ωs]
        seg = adaptive_autosplitbands2d(Iqω, qs, ωs)
        for support in (:components, :domains)
            masks = band_supports(seg, size(Iqω)...; support=support)
            @test length(masks) == seg.qfirst.nparents
            @test all(size(M) == size(Iqω) for M in masks)
            @test all(sum(masks) .<= 1)
            br = band_qfi(Iqω, seg, ωs, 10.0; support=support)
            @test br.total.FQ ≈ vec(sum(br.FQ; dims=1)) .+ br.FQres
            @test all(x -> isnan(x) || 0 <= x <= 1, br.ξ)
        end
        @test_throws ArgumentError band_supports(seg, size(Iqω)...; support=:bogus)
    end

    @testset "form factor and conversion to S_mag" begin
        # Reference values: periodictable (Brown, ITC Vol. C) Mn2+ <j0>.
        @test magnetic_form_factor(0.0) ≈ 0.9992 atol = 1e-12
        @test magnetic_form_factor(1.0) ≈ 0.9323149968067334 rtol = 1e-10
        @test magnetic_form_factor(2.0) ≈ 0.7629713782168471 rtol = 1e-10
        @test magnetic_form_factor(4.0) ≈ 0.37669841061323533 rtol = 1e-10

        I = [1.0 2.0; 3.0 NaN]; σ = [0.1 0.2; 0.3 0.4]; Q = [0.0, 2.0]
        S, σS = to_smag(I, Q; σ=σ)
        F0, F2 = magnetic_form_factor(0.0), magnetic_form_factor(2.0)
        @test S[1, :] ≈ 1.5 .* I[1, :] ./ F0^2
        @test S[2, 1] ≈ 1.5 * 3.0 / F2^2 && isnan(S[2, 2])
        @test σS ≈ 1.5 .* σ ./ [F0^2 F0^2; F2^2 F2^2]
        S1, _ = to_smag(I, Q; polarization=1, g=2, scale=2, form_factor=_ -> 1.0)
        @test S1[1, :] ≈ 2 .* I[1, :]
    end
end
