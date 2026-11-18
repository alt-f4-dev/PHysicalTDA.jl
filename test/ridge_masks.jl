using Test, PHysicalTDA

@testset "Parent-aware and measured-bin ridge recovery" begin
    q = collect(0.0:0.1:1.0)
    ω = collect(0.0:0.1:8.0)
    first_track = 4 .* abs.(q .- 0.5)
    second_track = 4 .+ 0.4 .* q
    I = [exp(-((e - a) / 0.15)^2) + 0.7exp(-((e - b) / 0.15)^2) - 2 for (a,b) in zip(first_track, second_track), e in ω]
    measured = trues(size(I))
    measured[8, :] .= false
    I[.!measured] .= Inf
    lower = BitMatrix([abs(e-a) <= 0.36 for a in first_track, e in ω])
    upper = BitMatrix([abs(e-b) <= 0.36 for b in second_track, e in ω])
    p = RidgeRecoveryParams(ω, q; max_slope=5, peak_smooth_radius=0)
    @test p.support_touch ≈ 0.5
    @test p.fallback_window >= p.support_touch
    @test p.min_prominence == 0.1
    @test RidgeRecoveryParams(ω, q; max_slope=5, min_prominence=0.23).min_prominence == 0.23
    @test_throws ArgumentError RidgeRecoveryParams(ω, reverse(q); max_slope=5)
    @test_throws ArgumentError RidgeRecoveryParams(ω, q; max_slope=0)
    t = recover_child_ridges(I, q, ω; parents=[upper, lower, falses(size(I))], valid=measured, params=p)
    @test t.global_children == [1, 1, 0]
    @test isempty(t.parent_groups[3])
    @test size(t.tracks[3]) == (0, length(q))
    @test t.tracks[2][1, 6] == 0
    @test all(isnan(t.tracks[k][1, 8]) for k in 1:2)
    @test all(t.local_multiplicity[:, 8] .== 0)
    @test !any(t.support_mask .& .!measured)
    @test !any(t.support_mask .& .!(upper .| lower))
    for (k, truth) in enumerate((second_track, first_track))
        good = findall(vec(measured[:, 1]))
        @test maximum(abs.(t.tracks[k][1, good] .- truth[good])) <= 0.11
    end
    @test length(child_ridge_dispersions(t).segments) == 4
    @test !any(t.graph_nodes[a].iq < 8 < t.graph_nodes[b].iq for (a,b) in t.graph_edges)
    @test_throws ArgumentError recover_child_ridges(I, q, ω; parents=[upper, upper], valid=measured)
    @test_throws DimensionMismatch recover_child_ridges(I, q, ω; valid=trues(3, 4))
    @test all(measured[n.iq, n.iω] for n in find_ridge_candidates(I, q, ω; params=p, valid=measured).nodes)

    two = recover_child_ridges(I, q, ω; parents=upper .| lower, valid=measured, params=p)
    @test two.global_children == [2]
    observed = findall(vec(measured[:, 1]))
    @test all(two.local_multiplicity[1, observed] .== 2)
    @test maximum(abs.(vec(minimum(two.tracks[1][:, observed]; dims=1)) .- first_track[observed])) <= 0.11
    @test maximum(abs.(vec(maximum(two.tracks[1][:, observed]; dims=1)) .- second_track[observed])) <= 0.11

    # An invalid neighbour must not turn a negative monotone spectrum into a peak.
    ramp = [-10.0 + j / 100 for _ in q, j in eachindex(ω)]
    valid = trues(size(ramp))
    valid[:, end-2:end] .= false
    ramp[.!valid] .= 0.0
    @test isempty(find_ridge_candidates(ramp, q, ω; params=p, valid).nodes)
    @test isempty(recover_child_ridges(ramp, q, ω; parents=trues(size(ramp)), valid, params=p).tracks[1])

    # A real segmentation-shaped result uses component masks in physical parent order.
    components = [(parent_band=1, linear_indices=findall(vec(upper))), (parent_band=2, linear_indices=findall(vec(lower)))]
    segmentation = (qfirst=(nparents=2,), components=components)
    from_segmentation = recover_child_ridges(I, q, ω; params=p, parents=segmentation, valid=measured)
    @test isequal(from_segmentation.tracks, t.tracks[1:2])
end

@testset "Fine-grid support refinement and uncertainty" begin
    q = collect(0.0:0.1:1.0)
    ω = collect(0.0:0.1:6.0)
    track = 1 .+ 2 .* q
    I = [exp(-((e - a) / 0.07)^2) + 0.8exp(-((e - 4.5) / 0.1)^2) for a in track, e in ω]
    σ = fill(0.02, size(I))
    valid = trues(size(I))
    valid[4, end] = false
    I[4, end] = NaN
    seed = falses(size(I))
    seed[1, 11] = true
    upper = BitMatrix([abs(e - 4.5) < 0.12 for _ in q, e in ω])
    kwargs = (resolution=(0.0, 0.1), bridge=0.201, minpix=4, close=0)
    result = refine_supports(I, σ, valid, q, ω, [seed, upper, falses(size(I))]; kwargs...)
    @test length(result) == 3
    @test !any(result[3])
    @test all(result[1][i, argmin(abs.(ω .- track[i]))] for i in eachindex(q))
    @test all(result[2][:, 46])
    @test !any(result[1] .& result[2])
    @test !any((result[1] .| result[2]) .& .!valid)
    @test !any(refine_supports(I, σ, valid, q, ω, [seed]; resolution=(0, 0), bridge=0.05, minpix=4)[1])
    @test_throws DimensionMismatch refine_supports(I, σ, valid, q, ω, [trues(2, 2)]; kwargs...)

    # Check normalized convolution and propagation through two passes against direct kernel sums.
    X = [2.0 3.0; 5.0 101.0]
    E = [0.2 0.3; 0.4 5.0]
    good = Bool[1 1; 1 0]
    grid = [0.0, 1.0]
    kernel = PHysicalTDA._ridge_kernel(grid, 2.0)
    smoothed, error = PHysicalTDA._smooth_ridge_measurements(X, E, good, kernel, kernel)
    w = [1.0, 0.5, 0.5]
    @test smoothed[1, 1] ≈ sum(w .* X[good]) / sum(w)
    @test error[1, 1] ≈ sqrt(sum(abs2, w .* E[good])) / sum(w)
    @test isnan(smoothed[2, 2])

    # Closing a small energy hole never bridges an unmeasured bin.
    X = repeat([0.0 5.0 0.0 5.0 0.0], 3, 1)
    mask = X .> 0
    good = trues(size(X))
    good[2, 3] = false
    filled = refine_supports(X, ones(size(X)), good, 1:3, 1:5, mask; resolution=(0, 0), close=1, minpix=1)
    @test filled[1][:, 3] == Bool[1, 0, 1]

    rb = rebin_to_resolution(I, q, ω; factors=(2, 2), remainder=:merge)
    coarse = trues(size(rb.I))
    via_rebin = refine_supports(I, σ, valid, q, ω, coarse; kwargs..., rebin=rb)
    manual = refine_supports(I, σ, valid, q, ω, upsample_mask(coarse, rb); kwargs...)
    @test via_rebin == manual
end

@testset "Trajectories on arbitrary masks" begin
    q, ω = [0.0, 0.1, 0.3, 0.5, 1.0, 2.0], [1.0, 2.0, 3.0]
    I = [1.0 0 0; 1 3 0; 0 100 0; 0 0 2; -2 -1 -3; 0 0 1]
    M = trues(size(I))
    M[1, :] .= false
    M[end, :] .= false
    valid = trues(size(I)); valid[3, :] .= false
    t = mask_dispersions(I, [M, falses(size(M))], q, ω; valid)
    @test t.q == q
    @test t.ω[1, 2:5] ≈ [1.75, 2.375, 3.0, 2.0]
    @test t.imputed[1, :] == Bool[1, 0, 1, 0, 1, 1]
    @test all(isnan, t.ω[2, :])
    @test all(t.imputed[2, :])
    @test all(isnan, t.ω[1, [1, 6]])
    b = mask_dispersions(I, M, q, ω; valid, centre=:bin)
    @test b.ω[1, 2:5] ≈ [2.0, 2.5, 3.0, 2.0]
    @test !b.imputed[1, 5]
    @test_throws ArgumentError mask_dispersions(I, M, q, ω; centre=:invalid)
end

@testset "Distinct child assignment" begin
    costs = [1.0 2 5 3; 1 4 8 2; 3 2 2 1]
    assignment = PHysicalTDA._ridge_minimum_assignment(costs)
    @test allunique(assignment)
    @test sum(costs[b, k] for (b,k) in pairs(assignment)) ==
          minimum(costs[1,a] + costs[2,b] + costs[3,c] for a in 1:4, b in 1:4, c in 1:4 if allunique((a,b,c)))
    @test PHysicalTDA._ridge_minimum_assignment(costs[:, 1:2]) == [1, 1, 2]
    more = [abs(b - 2k) for b in 2:2:12, k in 1:12]
    @test PHysicalTDA._ridge_minimum_assignment(more) == collect(1:6)
end
