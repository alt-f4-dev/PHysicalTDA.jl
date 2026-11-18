using LinearAlgebra, Sunny

function rg_histogram(centres, D, N; D2=copy(N), names=("[H,0,0]", "[0,K,0]", "DeltaE"), cell=(4.0, 4.0, 10.0))
    coords = [collect(Float64, axis) for axis in centres]
    edges = map(coords) do x
        length(x) == 1 && return [x[1] - 0.5, x[1] + 0.5]
        [(3x[1] - x[2]) / 2; (x[1:end-1] .+ x[2:end]) ./ 2; (3x[end] - x[end-1]) / 2]
    end
    units = [n == "DeltaE" ? "DeltaE" : "r.l.u." for n in names]
    frames = [n == "DeltaE" ? "General Frame" : "HKL" for n in names]
    lat = (a=cell[1], b=cell[2], c=cell[3], alpha=90.0, beta=90.0, gamma=90.0)
    meta = (; runs=[(; name="experiment0", lattice=lat)])
    return NormalizedHistogram(D, D2, N, edges, coords, collect(names), units, frames, meta)
end

@testset "Cuts as histograms" begin
    hs = collect(0.0:0.05:0.5)
    ω = [1.0, 2.0, 3.0, 4.0]
    n = length(hs)
    N = [1.0 + 0.1 * (a + 2b + c) for a in 1:n, b in 1:n, c in 1:4]
    h = rg_histogram((hs, hs, ω), 2 .* N .+ [0.01 * a * c for a in 1:n, b in 1:n, c in 1:4], N)
    path = [(0.0, 0.0), (0.5, 0.5)]
    step = sqrt(2) / 20
    cut = cut_path(h, path; width=0.01, step)
    hc = cut_histogram(cut)
    @test size(hc.D) == size(cut.values) && energy_axis(hc) == 2
    @test isequal(intensity(hc), cut.values)
    @test isequal(uncertainty(hc), cut.σ)
    coarse = rebin(hc; factors=(2, 1))
    direct = cut_path(h, path; width=0.01, step=2step)
    @test intensity(coarse) ≈ direct.values
    @test uncertainty(coarse) ≈ direct.σ
    @test size(rebin_to_resolution(hc; factors=(5, 2)).D) == (2, 2)
    folded = cut_histogram(cut_path(h, [(-0.5, -0.5), (0.5, 0.5)]; width=0.01, step, fold=[[1 0; 0 1], [-1 0; 0 -1]]))
    @test !isempty(folded.meta.numerator_covariance)
    @test_throws ArgumentError cut_histogram((; values=1))
end

@testset "Cartesian cut widths" begin
    hs = collect(0.0:0.05:0.5)
    ω = [1.0, 2.0]
    n = length(hs)
    N = ones(n, n, 2)
    h = rg_histogram((hs, hs, ω), 3 .* N, N)
    path = [(0.0, 0.0), (0.5, 0.5)]
    k = 2π / 4
    rlu = cut_path(h, path; width=0.02, step=sqrt(2) / 10)
    cart = cut_path(h, path; width=0.02k, step=k * sqrt(2) / 10, metric=:cartesian)
    @test sum(cart.weight) == sum(rlu.weight)   # same cylinder; bin-edge ties may round differently
    @test cart.s ≈ k .* rlu.s
    @test cart.s_units == "Å^-1" && rlu.s_units == "r.l.u."
    ls = collect(0.0:0.1:1.0)
    h3 = rg_histogram(([0.0], [0.0], ω, ls), ones(1, 1, 2, length(ls)), ones(1, 1, 2, length(ls));
                      names=("[H,0,0]", "[0,K,0]", "DeltaE", "[0,0,L]"))
    c3 = cut_path(h3, [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]; width=0.01, step=2π / 10 / 5, metric=:cartesian)
    @test c3.s_edges[end] ≈ 2π / 10
    @test count(>(0), c3.weight[:, 1]) == 5
    @test_throws ArgumentError cut_path(h, path; width=0.02, step=0.1, metric=:polar)
end

@testset "Point groups from a Sunny crystal" begin
    cryst = Crystal(lattice_vectors(4.2, 4.2, 13.8, 90, 90, 90), [[0, 0, 0]], 139)
    g3 = reciprocal_point_group(cryst)
    g2 = reciprocal_point_group(cryst; dims=2)
    @test length(g3) == 16 && length(g2) == 8
    @test [0 -1; 1 0] in g2 && [0 1; 1 0] in g2
    @test all(R -> abs(det(R)) == 1, g3)
    low = Crystal(lattice_vectors(4.0, 5.0, 6.0, 90, 90, 90), [[0, 0, 0]], 1)
    @test reciprocal_point_group(low) == [Matrix{Int}(I, 3, 3)]
    @test length(reciprocal_point_group(low; inversion=true)) == 2
    hs = collect(-0.5:0.05:0.5)
    n = length(hs)
    D = [1.0 + 0.3x + 0.2y^2 for x in hs, y in hs, _ in 1:2]
    h = rg_histogram((hs, hs, [1.0, 2.0]), D, ones(n, n, 2))
    path = [(0.0, 0.0), (0.5, 0.5)]
    viaC = cut_path(h, path; width=0.02, step=sqrt(2) / 10, fold=cryst)
    viaM = cut_path(h, path; width=0.02, step=sqrt(2) / 10, fold=g2)
    @test isequal(viaC.values, viaM.values) && viaC.weight == viaM.weight
    @test_throws ArgumentError reciprocal_point_group(cryst; dims=1)
end

@testset "Ridge gaps across rows without peaks" begin
    qs = collect(range(0, 1; length=41))
    ωs = collect(range(0, 8; length=60))
    truth = 1 .+ 6 .* qs
    I = [exp(-(ω - ω0)^2 / (2 * 0.15^2)) for ω0 in truth, ω in ωs]
    I[20, :] .= 0
    plain = build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=6)
    bridged = build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=6, max_gap_rows=1)
    @test plain.beta0 == 2 && bridged.beta0 == 1
    @test bridged.diagnostics[:gap_edges] == 1 && plain.diagnostics[:gap_edges] == 0
    I[19:21, :] .= 0
    @test build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=6, max_gap_rows=1).beta0 == 2
    @test build_ridge_topology(I, qs, ωs; smooth_radius_ω=1, max_jump_bins=6, max_gap_rows=3).beta0 == 1
    @test_throws ArgumentError build_ridge_topology(I, qs, ωs; max_gap_rows=-1)
    r = adaptive_autosplitbands2d(I, qs, ωs; verbose=false, ridge_max_gap_rows=3)
    @test r.ridge_diagnostics[:gap_edges] >= 1
end
