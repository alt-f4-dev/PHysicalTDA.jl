const FIXTURE_GRID = ([0.0, 2, 4], [-3.0, -1, 1, 3], [0.0, 1, 2, 3, 4], [0.25, 0.75])
const FIXTURE_NAMES = ("[H,0,0]", "[0,K,0]", "DeltaE", "[0,0,L]")
const FIXTURE_UNITS = ("r.l.u.", "r.l.u.", "DeltaE", "r.l.u.")
const FIXTURE_FRAMES = ("HKL", "HKL", "General Frame", "HKL")
const FIXTURE_UB = [1.0 2 3; 0 4 5; 0 0 6]
const FIXTURE_W = [1.0 0 0; 1 1 0; 0 0 1]

function write_mdhisto_fixture(path, signal; errors2=nothing, grid=FIXTURE_GRID, order=:forward,
                               axis_label=true, metadata=false, singleton_width=nothing, axis_keys=nothing)
    h5open(path, "w") do file
        entry = create_group(file, "MDHistoWorkspace")
        attributes(entry)["SaveMDVersion"] = Int32(2)
        attributes(entry)["QConvention"] = "Inelastic"
        group = create_group(entry, "data")
        dims = ntuple(identity, ndims(signal))
        perm = order === :forward ? dims : reverse(dims)
        keys = axis_keys === nothing ? map(i -> "D$(i - 1)", dims) : axis_keys
        group["signal"] = permutedims(signal, perm)
        if axis_label
            attributes(group["signal"])["axes"] = join((keys[i] for i in reverse(perm)), ':')
        end
        errors2 === nothing || (group["errors_squared"] = permutedims(errors2, perm))
        for i in dims
            key = keys[i]
            group[key] = grid[i]
            attrs = attributes(group[key])
            attrs["long_name"] = FIXTURE_NAMES[i]
            attrs["units"] = FIXTURE_UNITS[i]
            attrs["frame"] = FIXTURE_FRAMES[i]
            if length(grid[i]) == 1 && singleton_width !== nothing
                attrs["bin_width"] = singleton_width
            end
        end
        if metadata
            for i in (0, 2, 10)
                entry["experiment$i/logs/Ei/value"] = [31.0 - i]
                entry["experiment$i/logs/run_title/value"] = "fixture run $i"
            end
            run = entry["experiment0"]
            for (key, value) in zip(("a", "b", "c", "alpha", "beta", "gamma"), (4.2, 4.2, 13.7, 90.0, 90.0, 90.0))
                run["sample/oriented_lattice/unit_cell_$key"] = value
            end
            run["sample/oriented_lattice/orientation_matrix"] = permutedims(FIXTURE_UB)
            run["logs/W_MATRIX/value"] = vec(permutedims(FIXTURE_W))
            run["logs/goniometer/rotation_matrix"] = vec(permutedims(FIXTURE_W))
            run["logs/goniometer/num_axes"] = Int32(1)
            run["logs/goniometer/axis0/name"] = "omega"
            run["logs/goniometer/axis0/angle"] = 12.0
            run["logs/goniometer/axis0/rotationaxis"] = [0.0, 1.0, 0.0]
            attributes(run["logs/goniometer/axis0/angle"])["unit"] = "deg"
            attributes(run["logs/goniometer/axis0/angle"])["sense"] = "CCW"
            run["instrument/ignored"] = "not part of the metadata contract"
        end
    end
    return path
end

function write_synthetic_mdnorm(dir, D, D2, N; grid=FIXTURE_GRID, metadata=false, norm_order=:forward)
    path_D = write_mdhisto_fixture(joinpath(dir, "synthetic_D.nxs"), D; errors2=D2, grid, metadata)
    path_N = write_mdhisto_fixture(joinpath(dir, "synthetic_N.nxs"), N; grid, order=norm_order)
    I = map((d, n) -> n > 0 ? d / n : -1.0, D, N)
    path_I = write_mdhisto_fixture(joinpath(dir, "synthetic_I.nxs"), I; grid)
    return (; path_D, path_N, path_I)
end
