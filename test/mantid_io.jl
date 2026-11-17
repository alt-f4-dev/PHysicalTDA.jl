using HDF5
include("fixtures/mdnorm.jl")

@testset "Mantid MDHisto and MDNorm input" begin
    mktempdir() do dir
        D = reshape(collect(1.0:24.0), 2, 3, 4, 1)
        D2 = 4 .* D
        N = fill(2.0, size(D))
        N[1] = 0
        path_D, path_N, path_I = write_synthetic_mdnorm(dir, D, D2, N; metadata=true, norm_order=:reversed)
        raw = load_mdhisto(path_D)
        @test raw isa MDHistogram
        @test raw.signal == D
        @test raw.errors2 == D2
        @test raw.names == collect(FIXTURE_NAMES)
        @test raw.units == collect(FIXTURE_UNITS)
        @test raw.frames == collect(FIXTURE_FRAMES)
        @test raw.edges == collect(FIXTURE_GRID)
        @test raw.centres == [[1.0, 3.0], [-2.0, 0.0, 2.0], [0.5, 1.5, 2.5, 3.5], [0.5]]
        @test energy_axis(raw) == 3
        @test hkl_centres(raw) == [1 3 1 3 1 3; -2 -2 0 0 2 2; 0.5 0.5 0.5 0.5 0.5 0.5]
        @test axis_vector("[-H, 2.5H, 0]") == [-1.0, 2.5, 0.0]
        @test axis_vector("K") == [0, 1, 0]
        @test_throws ArgumentError axis_vector("[H,0,0,0]")

        @test lattice(raw) == (; a=4.2, b=4.2, c=13.7, alpha=90.0, beta=90.0, gamma=90.0)
        @test orientation_matrix(raw) == FIXTURE_UB
        @test projection(raw) == FIXTURE_W
        @test incident_energy(raw) == 31.0
        @test incident_energy(raw; run=2) == 29.0
        @test incident_energy(raw; run=1) === nothing
        @test run_title(raw; run=10) == "fixture run 10"
        @test [run.name for run in raw.meta.runs] == ["experiment0", "experiment2", "experiment10"]
        @test raw.meta.q_convention == "Inelastic"
        @test goniometer(raw).matrix == FIXTURE_W
        @test only(goniometer(raw).axes) == (; name="omega", direction=[0.0, 1.0, 0.0], angle=12.0, unit="deg", sense="CCW")
        @test lattice(load_mdhisto(path_D; metadata=false)) === nothing

        h = load_mdnorm(path_D, path_N; path_I)
        @test h isa NormalizedHistogram
        @test (h.D, h.D2, h.N) == (D, D2, N)
        @test measured(h) == (N .> 0)
        @test isnan(intensity(h)[1])
        @test isnan(uncertainty(h)[1])
        @test intensity(h)[2] == D[2] / 2
        @test uncertainty(h)[2] == sqrt(D2[2]) / 2
        valid = measured(h)
        @test sum(h.D[valid]) / sum(h.N[valid]) == sum(D[2:end]) / (2 * 23)
        @test sum(h.D2[valid]) / sum(h.N[valid])^2 == sum(D2[2:end]) / (2 * 23)^2

        bins = (2:2, :, 2:3, 1:1)
        subset = load_mdnorm(path_D, path_N; path_I, ranges=bins)
        @test size(subset.D) == (1, 3, 2, 1)
        @test subset.D == D[bins...]
        @test subset.D2 == D2[bins...]
        @test subset.N == N[bins...]
        @test subset.edges == [[2.0, 4.0], FIXTURE_GRID[2], [1.0, 2.0, 3.0], FIXTURE_GRID[4]]

        centres = Float64[]
        eachslice_mdnorm(path_D, path_N; dim=3, ranges=bins) do slab, centre
            push!(centres, centre)
            j = length(centres)
            @test slab.D == subset.D[:, :, j:j, :]
            @test slab.D2 == subset.D2[:, :, j:j, :]
            @test slab.N == subset.N[:, :, j:j, :]
            @test slab.edges[3] == subset.edges[3][j:j+1]
        end
        @test centres == [1.5, 2.5]
        @test_throws ArgumentError eachslice_mdnorm(identity, path_D, path_N; dim=5)
        @test_throws BoundsError load_mdhisto(path_D; ranges=(3:3, :, :, :))
        @test_throws ArgumentError load_mdhisto(path_D; ranges=(2:-1:1, :, :, :))
        @test_throws DimensionMismatch load_mdhisto(path_D; ranges=(:, :))

        without_errors = load_mdnorm(path_N, path_N)
        @test without_errors.D2 === nothing
        @test uncertainty(without_errors) === nothing
        eachslice_mdnorm(path_N, path_N; dim=4, metadata=false) do slab, centre
            @test slab.D == N
            @test slab.D2 === nothing
            @test centre == 0.5
        end
        bad_grid = ([1.0, 3, 5], FIXTURE_GRID[2:end]...)
        bad_N = write_mdhisto_fixture(joinpath(dir, "bad_N.nxs"), N; grid=bad_grid)
        @test_throws DimensionMismatch load_mdnorm(path_D, bad_N)
        bad_I = write_mdhisto_fixture(joinpath(dir, "bad_I.nxs"), ones(size(D)))
        @test_throws ArgumentError load_mdnorm(path_D, path_N; path_I=bad_I)
    end
end

@testset "Mantid axis layout and bin geometry" begin
    mktempdir() do dir
        path = joinpath(dir, "layout.nxs")
        signal = reshape(collect(1.0:8.0), 2, 4)
        centres = ([1.0, 3.0], [-2.0, 0.0, 2.0, 5.0])
        for order in (:forward, :reversed)
            write_mdhisto_fixture(path, signal; grid=centres, order, axis_label=false)
            raw = load_mdhisto(path)
            @test raw.signal == signal
            @test raw.centres == collect(centres)
            @test raw.edges == [[0.0, 2.0, 4.0], [-3.0, -1.0, 1.0, 3.5, 6.5]]
            @test load_mdhisto(path; ranges=(2:2, 3:4)).signal == signal[2:2, 3:4]
        end

        write_mdhisto_fixture(path, signal; grid=centres, axis_keys=("qy", "scanIndex"))
        named = load_mdhisto(path)
        @test named.signal == signal
        @test named.centres == collect(centres)

        square = reshape(collect(1.0:4.0), 2, 2)
        grid = (centres[1], centres[1])
        write_mdhisto_fixture(path, square; grid, axis_label=false)
        @test_throws ArgumentError load_mdhisto(path)
        write_mdhisto_fixture(path, square; grid, order=:reversed)
        @test load_mdhisto(path).signal == square
        write_mdhisto_fixture(path, square; grid=(centres[1], [1.0, 2, 3, 4, 5]))
        @test_throws DimensionMismatch load_mdhisto(path)
        write_mdhisto_fixture(path, square; grid, errors2=ones(1, 2))
        @test_throws DimensionMismatch load_mdhisto(path)

        write_mdhisto_fixture(path, [2.0]; grid=([5.0],))
        @test_throws ArgumentError load_mdhisto(path)
        write_mdhisto_fixture(path, [2.0]; grid=([5.0],), singleton_width=2.0)
        @test load_mdhisto(path).edges == [[4.0, 6.0]]
    end
end
