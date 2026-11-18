# Smoke test of the analysis layouts on the synthetic antiferromagnet. Runs when Makie is available in the
# test environment; the layouts live in the Makie extension.
using Test, PHysicalTDA, Statistics

if Base.get_extension(PHysicalTDA, :PHysicalTDAMakieExt) === nothing && Base.find_package("Makie") !== nothing
    @eval using Makie
end
const _layout_ext = Base.get_extension(PHysicalTDA, :PHysicalTDAMakieExt)
if _layout_ext === nothing
    @info "Makie not available; skipping the layout smoke test"
else
    include(joinpath(@__DIR__, "fixtures", "synthetic_afm.jl"))
    ext = _layout_ext
    @testset "Analysis layouts" begin
        files = write_synthetic_afm(mktempdir())
        h = load_mdnorm(files.A.D, files.A.N)
        hkl = corrected_indices(h, AFM.M, AFM.t)
        g = magnetic_conversion(h; l_average=(collect(-3.75:0.5:3.75), ones(16)))
        bk = (q, e) -> (afm_background(e), 0.0)
        cut = cut_path(h, [(0.0, 0.0), (0.5, 0.5)]; width=0.05, step=0.025 * sqrt(2), indices=hkl, conv=g, background=bk)
        κ = afm_true_scale()
        S, σ, ω = κ .* cut.values, κ .* cut.σ, cut.ω
        hs = cut.s ./ sqrt(2)
        r = qfi(S, ω, AFM.T; σ, σ_sys=κ .* cut.σ_background, ωrange=(0.45, Inf))
        low = falses(size(S)); low[:, ω .< 3] .= true
        br = band_qfi(S, [low, .!low], ω, AFM.T; σ, ωrange=(0.45, Inf))
        br2 = band_qfi(S, [low .& (ω' .< 2), .!(low .& (ω' .< 2))], ω, AFM.T; σ, ωrange=(0.45, Inf))
        sens = qfi_sensitivity(S, ω, AFM.T; σ, lower_limits=[0.35, 0.55], S=AFM.S)
        model = qfi([afm_smag(x, x, e) for x in hs, e in ω], ω, AFM.T; ωrange=(0.45, Inf)).FQ
        fig, ax = ext.plot_qfi_summary(r, hs; normalization=(; S=AFM.S), scale_σ=0.05, bands=br, sensitivity=sens,
                                       comparison=compare_band_qfi(br, br2), models=(; model), markers=(0.5,),
                                       band_labels=["low", "high"], xlabel="(h, h) (r.l.u.)")
        @test fig isa ext.Figure
        @test keys(ax) == (:FQ, :nqfi, :ξ, :sensitivity, :comparison)
        @test length(ext.plot_qfi_summary(r, hs; normalization=(; S=AFM.S))[2]) == 2
        @test_throws DimensionMismatch ext.plot_qfi_summary(r, hs[1:end-1]; normalization=(; S=AFM.S))

        I = replace(S, NaN => 0.0)
        I[:, ω .< 0.4] .= 0                    # the inelastic map only
        valid = cut.valid
        seg = adaptive_autosplitbands2d(I, hs, ω; verbose=false, max_slope=60.0)
        params = RidgeRecoveryParams(ω, hs; max_slope=60)
        cand = find_ridge_candidates(I, hs, ω; params, valid)
        children = recover_child_ridges(I, hs, ω; params, parents=seg, valid)
        refined = refine_supports(I, σ, valid, hs, ω, seg; resolution=(0.05, AFM.fwhm))
        traj = mask_dispersions(I, refined, hs, ω; valid)
        fig2, ax2 = ext.plot_ridge_workflow(I, hs, ω; valid, segmentation=seg, candidates=cand, children,
                                            supports=refined, trajectories=traj, markers=(0.5,), xlabel="(h, h) (r.l.u.)")
        @test fig2 isa ext.Figure
        @test keys(ax2) == (:data, :candidates, :parents, :children, :supports, :trajectories)
        @test keys(ext.plot_ridge_workflow(I, hs, ω)[2]) == (:data,)
        @test_throws DimensionMismatch ext.plot_ridge_workflow(I[1:end-1, :], hs, ω)
    end
end
