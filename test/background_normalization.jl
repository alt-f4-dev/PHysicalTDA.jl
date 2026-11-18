using LinearAlgebra, Statistics

const BN_CELL = (a=4.0, b=4.0, c=10.0, alpha=90.0, beta=90.0, gamma=90.0)

function bn_histogram(centres, D, N; D2=copy(N), names=("H", "K", "DeltaE"))
    coords = [collect(Float64, axis) for axis in centres]
    edges = map(coords) do x
        length(x) == 1 && return [x[1] - 0.5, x[1] + 0.5]
        [(3x[1] - x[2]) / 2; (x[1:end-1] .+ x[2:end]) ./ 2; (3x[end] - x[end-1]) / 2]
    end
    units = [n == "DeltaE" ? "DeltaE" : "r.l.u." for n in names]
    frames = [n == "DeltaE" ? "General Frame" : "HKL" for n in names]
    meta = (; runs=[(; name="experiment0", lattice=BN_CELL)])
    return NormalizedHistogram(D, D2, N, edges, coords, collect(names), units, frames, meta)
end

@testset "Exclusion masks and |Q|" begin
    hs = collect(-1.0:0.05:1.0)
    ω = [0.0, 1.0, 2.0]
    n = length(hs)
    h = bn_histogram((hs, hs, ω), ones(n, n, 3), ones(n, n, 3))
    Q = momentum_magnitude(h)
    @test size(Q) == (n, n)
    i, j = findfirst(≈(0.5), hs), findfirst(≈(0.25), hs)
    @test Q[i, j] ≈ 2π * sqrt(0.5^2 + 0.25^2) / 4
    mask = exclusion_mask(h; sites=[((0, 0), 0.1), ((0.5, 0.5), 0.05)])
    near(x, y) = mask[findfirst(≈(x), hs), findfirst(≈(y), hs), 1]
    @test near(0.0, 0.0) && near(1.0, -0.1) && near(-0.5, 0.5) && near(0.55, 0.45)
    @test !near(0.15, 0.0) && !near(0.5, 0.0) && !near(0.6, 0.5)
    @test all(mask[:, :, 1] .== mask[:, :, 3])
    disp = exclusion_mask(h; dispersion=q -> 1.0 + q[1], window=0.3)
    @test disp[findfirst(≈(0.0), hs), 1, 2] && !disp[findfirst(≈(0.0), hs), 1, 1]
    @test disp[findfirst(≈(0.9), hs), 1, 3] && !disp[findfirst(≈(0.5), hs), 1, 3]
    modes = exclusion_mask(h; dispersion=q -> (0.0, 2.0), window=0.1)
    @test all(modes[:, :, 1]) && !any(modes[:, :, 2]) && all(modes[:, :, 3])
    qcut = exclusion_mask(h; qmax=1.0)
    @test qcut == repeat(Q .> 1.0, 1, 1, 3)
    band = exclusion_mask(h; qbands=[(1.0, 0.05)])
    @test band[:, :, 2] == (abs.(Q .- 1.0) .<= 0.05)
    shifted = [(x + 0.5, y, 0.0) for x in hs, y in hs]
    @test exclusion_mask(h; sites=[((0, 0), 0.01)], indices=shifted)[findfirst(≈(-0.5), hs), findfirst(≈(0.0), hs), 1]
    @test_throws ArgumentError exclusion_mask(h; sites=[((0, 0), -0.1)])
end

@testset "Background profiles and |Q| shells" begin
    hs = collect(0.0:0.025:1.0)
    ω = collect(0.0:0.5:5.0)
    n, m = length(hs), length(ω)
    Q = [2π * sqrt(x^2 + y^2) / 4 for x in hs, y in hs]
    dq = 0.25
    truth(q, e) = 1.0 + 0.1 * min(floor(Int, q / dq), 6) + 0.02e
    band(x, y) = 2.0 + 2 * (x + y)
    N = fill(2.0, n, n, m)
    D = [N[a, b, c] * (truth(Q[a, b], ω[c]) + (abs(ω[c] - band(hs[a], hs[b])) <= 0.5 ? 5.0 : 0.0))
         for a in 1:n, b in 1:n, c in 1:m]
    D2 = fill(0.04, n, n, m)
    h = bn_histogram((hs, hs, ω), D, N; D2)
    excl = exclusion_mask(h; dispersion=q -> band(q[1], q[2]), window=0.5)

    prof = background_profile(h; exclude=excl)
    @test size(prof.value) == (1, m)
    @test prof(0.0, 0.0) == prof(10.0, 0.0)
    sel = .!excl[:, :, 1]
    @test prof.value[1, 1] ≈ sum(D[:, :, 1][sel]) / sum(N[:, :, 1][sel])
    @test prof.σ[1, 1] ≈ sqrt(sum(D2[:, :, 1][sel])) / sum(N[:, :, 1][sel])
    @test isequal(prof(0.0, 99.0), (NaN, NaN))

    b = background_shells(h; exclude=excl, dq, min_bins=5)
    @test b.qedges[1] == 0 && b.qedges[end] ≈ maximum(Q)
    ok = 0
    for c in 1:size(b.value, 1), j in 1:m
        b.source[c, j] in (:measured, :omega) || continue
        q = (b.qedges[c] + b.qedges[c+1]) / 2
        @test b.value[c, j] ≈ truth(q, ω[j]) atol = 1e-12
        ok += 1
    end
    @test ok >= 0.8 * length(b.value)
    @test any(==(:omega), b.source)
    v, σ = b(Q[3, 3], ω[2])
    @test v ≈ truth(Q[3, 3], ω[2]) && σ > 0
    @test isequal(b(maximum(Q) + 1, 1.0), (NaN, NaN))
    noω = background_shells(h; exclude=excl, dq, min_bins=5, interpolate=:q)
    @test !any(==(:omega), noω.source)

    # |Q| filling, the fallback and min_bins
    hole = copy(N); hole[:, :, 4] .= 0
    hh = bn_histogram((hs, hs, ω), D .* (hole .> 0), hole; D2)
    fb = background_shells(hh; dq, min_bins=5, interpolate=:q, fallback=7.5)
    @test all(fb.value[:, 4] .== 7.5) && all(fb.source[:, 4] .== :fallback)
    @test all(isnan, background_shells(hh; dq, min_bins=5, interpolate=:q).value[:, 4])
    sparse = background_shells(h; dq, min_bins=10^6, interpolate=:q, fallback=background_profile(h))
    @test all(==(:fallback), sparse.source)
    @test sparse.value[3, 2] ≈ background_profile(h).value[1, 2]
    one = background_shells(h; dq, qmax=0.4, min_bins=1)
    @test isequal(one(0.5, 0.0), (NaN, NaN)) && isfinite(one(0.3, 0.0)[1])
    flat = bn_histogram((hs, hs, ω), 3 .* N, N)
    sm = background_shells(flat; dq, min_bins=1, smooth=2.0, smooth_from=1.0)
    @test all(sm.value .≈ 3.0)
    @test occursin("TabulatedBackground", sprint(show, b))
end

@testset "Magnetic histograms reproduce converted cuts" begin
    hs = collect(0.0:0.05:0.5)
    ω = [1.0, 2.0, 3.0]
    n = length(hs)
    N = [1.0 + 0.1 * (a + b + c) for a in 1:n, b in 1:n, c in 1:3]
    D = 3 .* N
    h = bn_histogram((hs, hs, ω), D, N; D2=0.01 .* N)
    g = fill(0.7, size(D))
    bk = (q, e) -> (0.5 + 0.1q, 0.05)
    hm = magnetic_histogram(h, g; background=bk)
    path = [(0.0, 0.0), (0.5, 0.5)]
    direct = cut_path(h, path; width=0.03, step=0.1, conv=g, background=bk)
    viaS = cut_path(hm, path; width=0.03, step=0.1)
    @test isequal(direct.valid, viaS.valid)
    @test direct.values[direct.valid] ≈ viaS.values[viaS.valid]
    @test direct.σ[direct.valid] ≈ viaS.σ[viaS.valid]
    S, σS = to_smag(h, g; background=bk)
    @test intensity(hm) ≈ S && uncertainty(hm) ≈ σS
    mask = falses(size(D)); mask[1, 1, :] .= true
    @test !any(measured(magnetic_histogram(h, g; background=bk, exclude=mask))[1, 1, :])
    @test count(measured(magnetic_histogram(h, g; background=(q, e) -> e > 2.5 ? NaN : 0.0))) == 2n^2
end

@testset "Neutron cross sections" begin
    @test neutron_xs(:Mn).σ_inc == 0.4
    @test neutron_xs("Al").b_coh == 3.449
    @test neutron_xs(:V).σ_inc == 5.08
    rmf = formula_xs("Rb2MnF4")
    @test rmf.σ_inc ≈ 2 * 0.5 + 0.4 + 4 * 0.0008
    @test rmf.composition == Dict(:Rb => 2.0, :Mn => 1.0, :F => 4.0)
    @test formula_xs("Cu(NO3)2").composition == Dict(:Cu => 1.0, :N => 2.0, :O => 6.0)
    @test formula_xs("Sr14Cu24O41").σ_abs ≈ 14 * 1.28 + 24 * 3.78 + 41 * 0.00019
    @test formula_xs("La1.85Sr0.15CuO4").composition[:Sr] ≈ 0.15
    @test_throws ArgumentError neutron_xs(:Gd)
    @test_throws ArgumentError formula_xs("Cu(NO3")
    @test_throws ArgumentError formula_xs("cu")
end

@testset "Incoherent elastic scale" begin
    ω = collect(-3.0:0.05:6.0)
    fwhm, A, floor0 = 0.25, 2.0e-4, 3.0e-5
    σE = fwhm / (2 * sqrt(2 * log(2)))
    line = @. floor0 + A * exp(-ω^2 / (2σE^2)) / (sqrt(2π) * σE)
    @test elastic_fwhm(ω, line) ≈ fwhm rtol = 0.02
    N = ones(3, 3, length(ω))
    D = reshape(repeat(line, inner=9), 3, 3, :) .* N
    h = bn_histogram(([0.1, 0.2, 0.3], [0.1, 0.2, 0.3], ω), D, N)
    r = incoherent_scale(h; formula="Rb2MnF4")
    @test r.integral ≈ A rtol = 1e-3
    @test r.floor ≈ floor0
    @test r.window == (-1.0, 1.0)
    @test r.κ ≈ formula_xs("Rb2MnF4").σ_inc / (4π * r.integral * GAMMA_R0_SQ)
    @test incoherent_scale(h; σ_inc=1.0, magnetic_ions=2).κ ≈ 1 / (8π * r.integral * GAMMA_R0_SQ)
    @test incoherent_scale(h; σ_inc=1.0, fwhm=0.5).window == (-2.0, 1.25)
    @test_throws ArgumentError incoherent_scale(h)
    @test_throws ArgumentError incoherent_scale(h; formula="Rb2MnF4", σ_inc=1.0)
    c = container_corrected_scale(100.0, 0.5, 0.05)
    @test c.κ ≈ 200 && c.κ_lo ≈ 100 / 0.6
    @test isnan(container_corrected_scale(100.0, 0.99, 0.0).κ)
end

@testset "Aluminium powder lines" begin
    a = 4.032
    r111, r200 = powder_line_ratio(:Al, (1, 1, 1); a), powder_line_ratio(:Al, (2, 0, 0); a)
    @test r111 ≈ 241 rtol = 2e-3
    @test r200 ≈ 136 rtol = 4e-3
    @test PHysicalTDA._cubic_multiplicity((2, 2, 0)) == 12
    @test_throws ArgumentError powder_line_ratio(:Al, (1, 0, 0); a)
    τs = [2π * sqrt(3) / a, 4π / a]
    Q = collect(2.0:0.02:3.8)
    areas = [0.03, 0.017]
    R = @. 1e-3 + 2e-4 * (Q - 3) + areas[1] * exp(-(Q - τs[1] * 1.004)^2 / (2 * 0.06^2)) / (sqrt(2π) * 0.06) +
           areas[2] * exp(-(Q - τs[2] * 1.004)^2 / (2 * 0.06^2)) / (sqrt(2π) * 0.06)
    fit = fit_powder_lines(Q, R, τs; widths=0.05:0.001:0.07)
    @test fit.ε ≈ 0.004 atol = 1e-9
    @test fit.σQ ≈ 0.06
    @test fit_powder_lines(Q, R, τs).area ≈ areas rtol = 0.05
    @test fit.area ≈ areas rtol = 0.02
    flat = powder_line_flat(fit, [r111, r200])
    @test all(flat.detected)
    @test flat.flat ≈ mean(areas ./ [r111, r200]) rtol = 0.05
    @test flat.σflat >= std(flat.flats)
    Rnan = copy(R); Rnan[1:2:end] .= NaN
    @test fit_powder_lines(Q, Rnan, τs) === nothing
    @test isnan(powder_line_flat(nothing, [r111]).flat)

    hs = collect(-0.95:0.1:0.95)
    ls = collect(-4.75:0.5:4.75)
    ω = [-0.1, 0.0, 0.1]
    shape = (length(hs), length(hs), length(ω), length(ls))
    N = ones(shape)
    hq = bn_histogram((hs, hs, ω, ls), N, N; names=("H", "K", "DeltaE", "L"))
    prof = radial_elastic_profile(hq; energy_window=(-0.15, 0.15), dq=0.1, qmax=4.0)
    @test all(x -> isnan(x) || x ≈ 0.3, prof.R)
    @test length(prof.Q) == 40
end

@testset "Sum-rule scales" begin
    ωs = collect(0.25:0.25:8.0)
    qs = collect(0.05:0.05:0.5)
    T = 5.0
    model = [exp(-(w - (1 + 8q))^2 / 0.2) for q in qs, w in ωs]
    κtrue = 3000.0
    S = model ./ κtrue
    σ = fill(1e-7, size(S))
    r = first_moment_scale(S, σ, ωs, model; T, window=(0.5, 7.5), correction=0.96, correction_σ=0.01)
    @test r.κ ≈ 0.96 * κtrue
    @test r.κ_lo ≈ 0.95 * κtrue
    @test all(r.used)
    S[3, 10] = NaN
    kT = KB_MEV_PER_K * T
    full = [sum(w * (1 - exp(-w / kT)) * x * 0.25 for (w, x) in zip(ωs, model[i, :])) for i in axes(model, 1)]
    r2 = first_moment_scale(S, σ, ωs, model; T, window=(0.5, 4.0), model_full=full)
    @test r2.used[1] && !r2.used[3] && isnan(r2.κ_q[3]) && !r2.used[end]
    @test r2.κ ≈ κtrue
    @test count(first_moment_scale(S, fill(1.0, size(S)), ωs, model; T, window=(0.5, 7.5)).used) == 0
    zone = cat([model .* (1 + 0.1k) for k in 1:3]...; dims=3)
    zone = permutedims(zone, (1, 3, 2))
    Sz = zone ./ κtrue
    Sz[1, 1, :] .= NaN
    t = total_moment_scale(Sz, zone, ωs; window=(0.5, 7.5), correction=0.985, correction_σ=0.015)
    @test t.κ ≈ 0.985 * κtrue
    @test t.κ_lo ≈ 0.97 * κtrue
    @test 0.9 < t.covered_fraction < 1
    table = compare_scales([(name="incoherent", κ=1700.0), (name="first moment", κ=r.κ, κ_lo=r.κ_lo),
                            (name="total moment", κ=t.κ, κ_lo=t.κ_lo)]; reference="total moment")
    @test table[3].ratio ≈ 1 && table[1].used == 1700
    @test select_scale(table).name == "incoherent"
    @test select_scale(table; among=("first moment", "total moment")).name == "first moment"
    @test_throws ArgumentError compare_scales([(name="a", κ=1.0)]; reference="b")
end
