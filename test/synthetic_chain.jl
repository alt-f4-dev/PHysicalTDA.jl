using Test, PHysicalTDA, LinearAlgebra, Statistics, Sunny
include(joinpath(@__DIR__, "fixtures", "synthetic_afm.jl"))

@testset "Synthetic antiferromagnet: reduction chain" begin
    dir = mktempdir()
    files = write_synthetic_afm(dir)
    truth = files.truth
    κtrue = afm_true_scale()
    hA = load_mdnorm(files.A.D, files.A.N)
    @test lattice(hA).a == truth.a

    # Indexing
    sites = [(x, y) for x in (-0.5, 0.5) for y in (-0.5, 0.5)]
    peaks = find_bragg_peaks(hA; positions=sites, energy_window=(-0.3, 0.3))
    fit = fit_index_correction(peaks)
    @test fit.M ≈ truth.M atol = 5e-3
    @test fit.t ≈ truth.t atol = 2e-3
    hkl = corrected_indices(hA, fit.M, fit.t)

    # Background
    disp(q) = afm_dispersion(q[1], q[2])
    bragg = [((0, 0), 0.15), ((0.5, 0.5), 0.15)]
    m_prof = exclusion_mask(hA; sites=bragg, dispersion=disp, window=1.0, qmax=2.5, indices=hkl)
    m_shell = m_prof .| exclusion_mask(hA; sites=[((0, 0), 0.2), ((0.5, 0.5), 0.2)], indices=hkl)
    b0 = background_profile(hA; exclude=m_prof, fill=true)
    # The background peaks at 6 meV, under the excluded top of the magnon band; bridging that gap along energy
    # (:omega_then_q) would cut through the peak, so missing cells are taken from neighbouring shells instead.
    b = background_shells(hA; exclude=m_shell, dq=0.3, qmax=2.5, min_bins=20, smooth=0.5, smooth_from=1.0, fallback=b0,
                          interpolate=:q)
    ω = hA.centres[3]
    high = findall(x -> 7.0 <= x <= 8.5, ω)
    @test maximum(abs(b(1.0, ω[j])[1] / afm_background(ω[j]) - 1) for j in high) < 0.05

    # Absolute scales
    inc = incoherent_scale(b0; formula="Rb2MnF4")
    @test inc.κ ≈ κtrue / (1 + truth.fAl) rtol = 0.05
    τs = [2π * norm(x) / truth.aAl for x in ((1, 1, 1), (2, 0, 0))]
    radial = radial_elastic_profile(files.B.D, files.B.N; dim=4, energy_window=inc.window, dq=0.02, qmax=3.5)
    lines = fit_powder_lines(radial.Q, radial.R, τs)
    al = powder_line_flat(lines, [powder_line_ratio(:Al, x; a=truth.aAl) for x in ((1, 1, 1), (2, 0, 0))])
    profB = background_profile(files.B.D, files.B.N; dim=4,
                               exclude=slice -> exclusion_mask(slice; qbands=[(τ, 0.3) for τ in τs]))
    incB = incoherent_scale(profB; formula="Rb2MnF4", fwhm=inc.fwhm)
    share = al.flat / incB.integral
    @test share ≈ truth.fAl / (1 + truth.fAl) atol = 0.05
    κAl = container_corrected_scale(incB.κ, share, al.σflat / incB.integral)
    @test κAl.κ ≈ κtrue rtol = 0.06

    # Conversion and cut
    hB = load_mdnorm(files.B.D, files.B.N)
    g = magnetic_conversion(hA; ion=:Mn2, l_average=(hB.centres[4], vec(sum(hB.N; dims=(1, 2, 3)))))
    gB = magnetic_conversion(hB; ion=:Mn2, average_over=4)
    onA = map_nearest_bins(gB, hB, hA)
    near_half = [abs(x[1] - 0.5) < 0.1 && abs(x[2] - 0.5) < 0.1 for x in hkl]
    @test maximum(abs(g[i, k, 40, 1] / onA.values[i, k, 40, 1] - 1) for i in axes(g, 1), k in axes(g, 2) if near_half[i, k]) < 0.05
    C4 = reciprocal_point_group(Sunny.Crystal(Sunny.lattice_vectors(truth.a, truth.a, truth.c, 90, 90, 90), [[0, 0, 0]], 139); dims=2)
    m_cut = exclusion_mask(hA; sites=[((0, 0), 0.05)], qmax=2.5, indices=hkl)
    cut = cut_path(hA, [(0.0, 0.0), (0.5, 0.5)]; width=0.05, step=0.05 * sqrt(2), fold=C4, indices=hkl,
                   conv=g, background=b, exclude=m_cut)
    hs = cut.s ./ sqrt(2)
    model = [afm_smag(x, x, e) for x in hs, e in ω]
    m1 = first_moment_scale(cut.values, cut.σ, ω, model; T=truth.T, window=(0.45, 7.5))
    @test m1.κ ≈ κtrue rtol = 0.03
    @test count(m1.used) >= length(hs) ÷ 2

    zone = fold_to_zone(magnetic_histogram(hA, g; background=b, exclude=m_cut), C4; grid=(0:0.05:1, 0:0.05:1), indices=hkl)
    xc = collect(0.025:0.05:0.975)
    model_zone = [afm_smag(x, y, e) for x in xc, y in xc, e in ω]
    m0 = total_moment_scale(intensity(zone), model_zone, ω; window=(0.45, 7.5))
    @test m0.κ ≈ κtrue rtol = 0.08
    @test m0.covered_fraction > 0.9

    # With the generating background, the scales isolate the conversion, cut and sum-rule machinery.
    exact = (q, e) -> (afm_background(e), 0.0)
    cutx = cut_path(hA, [(0.0, 0.0), (0.5, 0.5)]; width=0.05, step=0.05 * sqrt(2), fold=C4, indices=hkl,
                    conv=g, background=exact, exclude=m_cut)
    @test first_moment_scale(cutx.values, cutx.σ, ω, model; T=truth.T, window=(0.45, 7.5)).κ ≈ κtrue rtol = 0.03
    zonex = fold_to_zone(magnetic_histogram(hA, g; background=exact, exclude=m_cut), C4; grid=(0:0.05:1, 0:0.05:1), indices=hkl)
    @test total_moment_scale(intensity(zonex), model_zone, ω; window=(0.45, 7.5)).κ ≈ κtrue rtol = 0.03

    # QFI against the truth near the antiferromagnetic zone centre
    r = qfi(κtrue .* cut.values, ω, truth.T; σ=κtrue .* cut.σ, σ_sys=κtrue .* cut.σ_background, ωrange=(0.45, Inf))
    rtrue = qfi(model, ω, truth.T; ωrange=(0.45, Inf))
    i = argmin(abs.(hs .- 0.45))
    @test r.FQ[i] ≈ rtrue.FQ[i] rtol = 0.08
    @test all(r.σsys .>= 0) && any(r.σsys .> 0)
    hc = cut_histogram(cut)
    @test isequal(intensity(hc), cut.values)
end
