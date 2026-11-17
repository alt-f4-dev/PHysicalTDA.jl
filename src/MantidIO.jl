const _MDHISTO_ENTRY = "/MDHistoWorkspace"
const _MDHISTO_DATA = _MDHISTO_ENTRY * "/data"

"""
    MDHistogram

Mantid `SaveMD` histogram in `D0, D1, …` order, or on named grids in the axis order presented by HDF5.jl. The raw `signal` and optional `errors2` arrays retain all axes, including single-bin integrated dimensions. Coordinate vectors, names, units, frames and experiment records are available as `edges`, `centres`, `names`, `units`, `frames` and `meta`.
"""
struct MDHistogram{N,M<:NamedTuple}
    signal::Array{Float64,N}
    errors2::Union{Nothing,Array{Float64,N}}
    edges::Vector{Vector{Float64}}
    centres::Vector{Vector{Float64}}
    names::Vector{String}
    units::Vector{String}
    frames::Vector{String}
    meta::M
end

"""
    NormalizedHistogram

An MDNorm histogram with numerator `D`, optional numerator variance `D2`, and normalization `N` on the same grid. Use `intensity`, `uncertainty`, and `measured` for derived arrays.
"""
struct NormalizedHistogram{R,M<:NamedTuple}
    D::Array{Float64,R}
    D2::Union{Nothing,Array{Float64,R}}
    N::Array{Float64,R}
    edges::Vector{Vector{Float64}}
    centres::Vector{Vector{Float64}}
    names::Vector{String}
    units::Vector{String}
    frames::Vector{String}
    meta::M
end

const _MantidHistogram = Union{MDHistogram,NormalizedHistogram}

function _read_h5_if_present(parent, name::AbstractString)
    haskey(parent, name) || return nothing
    return _with_h5_object(parent, name) do obj
        obj isa HDF5.Dataset ? read(obj) : nothing
    end
end

_h5_attr(obj, name, default) = haskey(HDF5.attributes(obj), name) ? HDF5.read_attribute(obj, name) : default

function _md_axis_matches(count::Int, coordinates::Vector{Float64})
    return length(coordinates) == count || length(coordinates) == count + 1
end

function _md_axis_keys(rank::Int, label)
    canonical = ["D$(i - 1)" for i in 1:rank]
    label === nothing && return canonical, nothing
    # NeXus lists axes slowest first, while HDF5.jl presents dimensions fastest first.
    stored = reverse(strip.(split(String(label), ':')))
    length(stored) == rank && all(!isempty, stored) && allunique(stored) ||
        throw(ArgumentError("invalid MDHisto signal axes: $label"))
    if all(key -> occursin(r"^D[0-9]+$", key), stored)
        issetequal(stored, canonical) || throw(ArgumentError("invalid MDHisto signal axes: $label"))
        return canonical, stored
    end
    # Older SaveMD v2 files use dimension IDs as dataset names.
    return stored, stored
end

function _md_axis_order(shape::NTuple{R,Int}, axes, axis_keys, stored) where {R}
    forward = ntuple(identity, Val(R))
    backward = reverse(forward)
    matches(order) = all(i -> _md_axis_matches(shape[order[i]], axes[i].coordinates), forward)
    if stored !== nothing
        order = ntuple(i -> findfirst(==(axis_keys[i]), stored)::Int, Val(R))
        matches(order) || throw(DimensionMismatch("MDHisto signal shape $shape disagrees with its axes attribute"))
        return order
    end
    direct, reversed = matches(forward), matches(backward)
    direct && reversed && R > 1 && throw(ArgumentError("ambiguous MDHisto axis order; signal needs an axes attribute"))
    direct && return forward
    reversed && return backward
    throw(DimensionMismatch("MDHisto signal shape $shape matches neither D0-first nor reversed axis order"))
end

function _md_axis_geometry(axis, n::Int)
    coordinates = axis.coordinates
    if length(coordinates) == n + 1
        return coordinates, [(coordinates[i] + coordinates[i + 1]) / 2 for i in 1:n]
    elseif n == 1
        width = axis.bin_width
        width isa Real && isfinite(width) && width > 0 ||
            throw(ArgumentError("singleton centre axis $(axis.name) requires a positive bin_width attribute"))
        centre = only(coordinates)
        return [centre - width / 2, centre + width / 2], coordinates
    end
    edges = Vector{Float64}(undef, n + 1)
    for i in 2:n
        edges[i] = (coordinates[i - 1] + coordinates[i]) / 2
    end
    edges[1] = 2 * first(coordinates) - edges[2]
    edges[end] = 2 * last(coordinates) - edges[end - 1]
    return edges, coordinates
end

function _md_layout(group::HDF5.Group)
    shape, label = _with_h5_object(group, "signal") do dataset
        size(dataset), _h5_attr(dataset, "axes", nothing)
    end
    return _md_layout(group, shape, label)
end

function _md_layout(group::HDF5.Group, shape::NTuple{R,Int}, label) where {R}
    R > 0 && all(>(0), shape) || throw(ArgumentError("MDHisto signal must be a nonempty array"))
    axis_keys, stored = _md_axis_keys(R, label)
    axes = map(axis_keys) do key
        haskey(group, key) || throw(ArgumentError("MDHisto axis dataset '$key' is missing"))
        _with_h5_object(group, key) do dataset
            coordinates = read(dataset)
            coordinates isa AbstractVector{<:Real} || throw(ArgumentError("MDHisto axis $key must be a numeric vector"))
            all(isfinite, coordinates) && all(i -> coordinates[i] > coordinates[i - 1], 2:length(coordinates)) ||
                throw(ArgumentError("MDHisto axis $key must be finite and strictly increasing"))
            (; coordinates=convert(Vector{Float64}, coordinates), name=String(_h5_attr(dataset, "long_name", key)),
               unit=String(_h5_attr(dataset, "units", "")), frame=String(_h5_attr(dataset, "frame", "")),
               bin_width=_h5_attr(dataset, "bin_width", nothing))
        end
    end
    order = _md_axis_order(shape, axes, axis_keys, stored)
    counts = ntuple(i -> shape[order[i]], Val(R))
    geometry = map(i -> _md_axis_geometry(axes[i], counts[i]), 1:R)
    errors2 = haskey(group, "errors_squared")
    if errors2
        _with_h5_object(group, "errors_squared") do dataset
            size(dataset) == shape || throw(DimensionMismatch("errors_squared shape does not match signal"))
        end
    end
    return (; counts, order, errors2, edges=first.(geometry), centres=last.(geometry),
            names=getproperty.(axes, :name), units=getproperty.(axes, :unit), frames=getproperty.(axes, :frame))
end

_md_ranges(layout, ::Nothing) = nothing

function _md_ranges(layout, ranges)
    length(ranges) == length(layout.counts) || throw(DimensionMismatch("ranges must have one entry per workspace axis"))
    return ntuple(length(layout.counts)) do i
        range = ranges[i]
        range isa Colon && return 1:layout.counts[i]
        range isa AbstractUnitRange{<:Integer} && !isempty(range) ||
            throw(ArgumentError("ranges[$i] must be : or a nonempty unit range"))
        1 <= first(range) <= last(range) <= layout.counts[i] || throw(BoundsError(1:layout.counts[i], range))
        Int(first(range)):Int(last(range))
    end
end

function _md_read(group::HDF5.Group, key, order, ranges)
    return _with_h5_object(group, key) do dataset
        _md_read(dataset, order, ranges)
    end
end

function _md_read(dataset::HDF5.Dataset, order::NTuple{R,Int}, ranges) where {R}
    forward = ntuple(identity, Val(R))
    stored_ranges = ranges === nothing ? nothing : ntuple(i -> ranges[findfirst(==(i), order)::Int], Val(R))
    values = (stored_ranges === nothing ? read(dataset) : dataset[stored_ranges...])::Array{<:Real,R}
    return convert(Array{Float64,R}, order == forward ? values : permutedims(values, order))
end

function _md_coordinates(layout, ranges)
    ranges === nothing && return layout.edges, layout.centres
    edges = [layout.edges[i][first(r):last(r)+1] for (i, r) in enumerate(ranges)]
    centres = [layout.centres[i][r] for (i, r) in enumerate(ranges)]
    return edges, centres
end

function _md_validate_grid(a, b, label)
    for field in (:edges, :centres, :names, :units, :frames)
        getproperty(a, field) == getproperty(b, field) || throw(DimensionMismatch("MDNorm $label grid differs in $field"))
    end
end

function _md_pack_normalized(D, D2, N, layout, ranges, meta)
    edges, centres = _md_coordinates(layout, ranges)
    return NormalizedHistogram(D, D2, N, edges, centres, layout.names, layout.units, layout.frames, meta)
end

function _md_read_normalized(data, norm, d, n, ranges, meta)
    D = _md_read(data, "signal", d.order, ranges)
    D2 = d.errors2 ? _md_read(data, "errors_squared", d.order, ranges) : nothing
    N = _md_read(norm, "signal", n.order, ranges)
    return _md_pack_normalized(D, D2, N, d, ranges, meta)
end

"""
    load_mdhisto(path; metadata=true, ranges=nothing) -> MDHistogram

Read the `/MDHistoWorkspace/data` layout written by Mantid `SaveMD` v2. In workspace axis order, each entry of `ranges` is a contiguous bin range or `Colon()` for the full axis. Singleton ranges retain their dimensions. Use `metadata=false` to skip experiment records.
"""
function load_mdhisto(path::AbstractString; metadata::Bool=true, ranges=nothing)
    h5open(path, "r") do file
        _with_h5_object(file, _MDHISTO_DATA) do group
            layout = _md_layout(group)
            bins = _md_ranges(layout, ranges)
            signal = _md_read(group, "signal", layout.order, bins)
            errors2 = layout.errors2 ? _md_read(group, "errors_squared", layout.order, bins) : nothing
            edges, centres = _md_coordinates(layout, bins)
            meta = metadata ? _md_metadata(file) : (; runs=NamedTuple[], q_convention=nothing)
            return MDHistogram(signal, errors2, edges, centres, layout.names, layout.units, layout.frames, meta)
        end
    end
end

function _with_mdnorm(f, path_D, path_N, metadata)
    h5open(path_D, "r") do data_file
        h5open(path_N, "r") do norm_file
            _with_h5_object(data_file, _MDHISTO_DATA) do data
                _with_h5_object(norm_file, _MDHISTO_DATA) do norm
                    d, n = _md_layout(data), _md_layout(norm)
                    _md_validate_grid(d, n, "_D/_N")
                    meta = metadata ? _md_metadata(data_file) : (; runs=NamedTuple[], q_convention=nothing)
                    return f(data, norm, d, n, meta)
                end
            end
        end
    end
end

"""
    load_mdnorm(path_D, path_N; path_I=nothing, ranges=nothing, metadata=true) -> NormalizedHistogram

Read matching MDNorm numerator and normalization workspaces. If `path_I` is given, also check the saved intensity against `D/N` on the selected measured bins. Experiment metadata comes from the numerator file.
"""
function load_mdnorm(path_D::AbstractString, path_N::AbstractString; path_I::Union{Nothing,AbstractString}=nothing,
                     ranges=nothing, metadata::Bool=true)
    _with_mdnorm(path_D, path_N, metadata) do data, norm, d, n, meta
        bins = _md_ranges(d, ranges)
        h = _md_read_normalized(data, norm, d, n, bins, meta)
        if path_I !== nothing
            h5open(path_I, "r") do file
                _with_h5_object(file, _MDHISTO_DATA) do group
                    saved = _md_layout(group)
                    _md_validate_grid(d, saved, "_D/_I")
                    I = _md_read(group, "signal", saved.order, bins)
                    for i in eachindex(h.D, h.N, I)
                        if h.N[i] > 0 && !isapprox(h.D[i] / h.N[i], I[i]; rtol=1e-6, atol=0)
                            throw(ArgumentError("MDNorm _I differs from D/N at selected bin $i"))
                        end
                    end
                end
            end
        end
        return h
    end
end

"""
    eachslice_mdnorm(f, path_D, path_N; dim, ranges=nothing, metadata=true)

Call `f(slice, centre)` once per selected bin of axis `dim`, in increasing bin order. Each slice retains `dim` as a singleton dimension. Files stay open during iteration, and only one slice is needed in memory unless the callback retains the results.
"""
function eachslice_mdnorm(f::F, path_D::AbstractString, path_N::AbstractString; dim::Integer, ranges=nothing,
                          metadata::Bool=true) where {F}
    _with_mdnorm(path_D, path_N, metadata) do data, norm, d, n, meta
        1 <= dim <= length(d.counts) || throw(ArgumentError("dim must index a workspace axis"))
        bins = _md_ranges(d, ranges)
        bins === nothing && (bins = map(i -> 1:i, d.counts))
        _with_h5_object(data, "signal") do D
            _with_h5_object(norm, "signal") do N
                if d.errors2
                    _with_h5_object(data, "errors_squared") do D2
                        _md_foreach_slice(f, D, D2, N, d, n, bins, meta, Int(dim))
                    end
                else
                    _md_foreach_slice(f, D, nothing, N, d, n, bins, meta, Int(dim))
                end
            end
        end
    end
    return nothing
end

function _md_foreach_slice(f::F, data, errors2, norm, d, n, bins::NTuple{R,UnitRange{Int}}, meta, dim::Int) where {F,R}
    for j in bins[dim]
        slice_bins = ntuple(i -> i == dim ? (j:j) : bins[i], Val(R))
        D = _md_read(data, d.order, slice_bins)
        D2 = errors2 === nothing ? nothing : _md_read(errors2, d.order, slice_bins)
        N = _md_read(norm, n.order, slice_bins)
        f(_md_pack_normalized(D, D2, N, d, slice_bins, meta), d.centres[dim][j])
    end
    return nothing
end

"""
    intensity(h::NormalizedHistogram)

Return `D/N` where `N > 0`, and `NaN` elsewhere.
"""
intensity(h::NormalizedHistogram) = map((D, N) -> N > 0 ? D / N : NaN, h.D, h.N)

"""
    uncertainty(h::NormalizedHistogram)

Return `sqrt(D2)/N` where `N > 0`, and `NaN` elsewhere. Return `nothing` when the numerator variance is unavailable.
"""
function uncertainty(h::NormalizedHistogram)
    h.D2 === nothing && return nothing
    return map((D2, N) -> N > 0 ? sqrt(D2) / N : NaN, h.D2, h.N)
end

"""
    measured(h::NormalizedHistogram)

Return the MDNorm coverage mask, `N .> 0`.
"""
measured(h::NormalizedHistogram) = h.N .> 0

_md_scalar(value) = value isa AbstractArray && length(value) == 1 ? only(value) : value

function _md_matrix(value, path)
    value === nothing && return nothing
    value isa AbstractArray{<:Real} && (size(value) == (3, 3) || size(value) == (9,)) ||
        throw(ArgumentError("$path must hold a 3×3 matrix or nine row-major entries"))
    # Mantid serializes matrix rows consecutively and HDF5.jl reads the first dimension fastest.
    return permutedims(reshape(Float64.(value), 3, 3))
end

function _md_goniometer(run)
    path = "logs/goniometer"
    haskey(run, path) || return nothing
    return _with_h5_object(run, path) do group
        matrix = _md_matrix(_read_h5_if_present(group, "rotation_matrix"), "goniometer/rotation_matrix")
        n = _md_scalar(_read_h5_if_present(group, "num_axes"))
        axes = map(0:(something(n, 0) - 1)) do i
            _with_h5_object(group, "axis$i") do axis
                name = _md_scalar(_read_h5_if_present(axis, "name"))
                direction = _read_h5_if_present(axis, "rotationaxis")
                angle, unit, sense = _with_h5_object(axis, "angle") do dataset
                    _md_scalar(read(dataset)), _h5_attr(dataset, "unit", nothing), _h5_attr(dataset, "sense", nothing)
                end
                (; name, direction, angle, unit, sense)
            end
        end
        (; matrix, axes)
    end
end

function _md_log(run, name)
    path = "logs/$name"
    haskey(run, path) || return nothing
    return _with_h5_object(run, path) do obj
        if obj isa HDF5.Dataset
            return _md_scalar(read(obj))
        end
        return _md_scalar(_read_h5_if_present(obj, "value"))
    end
end

const _MD_CELL_KEYS = (:a, :b, :c, :alpha, :beta, :gamma)

function _md_run(entry, name::String)
    return _with_h5_object(entry, name) do run
        cell = map(_MD_CELL_KEYS) do key
            _md_scalar(_read_h5_if_present(run, "sample/oriented_lattice/unit_cell_$key"))
        end
        lattice = all(isnothing, cell) ? nothing : NamedTuple{_MD_CELL_KEYS}(cell)
        ub = _md_matrix(_read_h5_if_present(run, "sample/oriented_lattice/orientation_matrix"), "orientation_matrix")
        projection = _md_matrix(_md_log(run, "W_MATRIX"), "W_MATRIX")
        (; name, lattice, orientation_matrix=ub, incident_energy=_md_log(run, "Ei"), projection,
           goniometer=_md_goniometer(run), run_title=_md_log(run, "run_title"))
    end
end

function _md_metadata(file::HDF5.File)
    return _with_h5_object(file, _MDHISTO_ENTRY) do entry
        names = filter(name -> occursin(r"^experiment\d+$", name), keys(entry))
        sort!(names; by=name -> parse(Int, SubString(name, 11)))
        (; runs=[_md_run(entry, name) for name in names], q_convention=_h5_attr(entry, "QConvention", nothing))
    end
end

function _md_run_field(h::_MantidHistogram, field::Symbol, run::Integer)
    run >= 0 || throw(ArgumentError("run must be a nonnegative experiment index"))
    name = "experiment$run"
    index = findfirst(record -> record.name == name, h.meta.runs)
    return index === nothing ? nothing : getproperty(h.meta.runs[index], field)
end

"""
    lattice(h; run=0)

Return Mantid unit-cell parameters `(a, b, c, alpha, beta, gamma)` for `experiment<run>`, or `nothing` if absent. Lengths are in ångströms and angles in degrees.
"""
lattice(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :lattice, run)

"""
    orientation_matrix(h; run=0)

Return the stored Mantid UB matrix for `experiment<run>` in row/column order, or `nothing` if absent.
"""
orientation_matrix(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :orientation_matrix, run)

"""
    incident_energy(h; run=0)

Return the incident energy from the Mantid `Ei` log of `experiment<run>`, or `nothing` if absent. Multi-value logs are retained as arrays.
"""
incident_energy(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :incident_energy, run)

"""
    projection(h; run=0)

Return the `W_MATRIX` projection for `experiment<run>` in row/column order, or `nothing` if absent.
"""
projection(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :projection, run)

"""
    goniometer(h; run=0)

Return the stored rotation matrix and ordered motor records (angle, unit, rotation axis, and sense) for `experiment<run>`, or `nothing` if absent.
"""
goniometer(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :goniometer, run)

"""
    run_title(h; run=0)

Return the title from the `run_title` log of `experiment<run>`, or `nothing` if absent.
"""
run_title(h::_MantidHistogram; run::Integer=0) = _md_run_field(h, :run_title, run)

"""
    energy_axis(h)

Return the index of the unique axis labelled `DeltaE` (by name or unit).
"""
function energy_axis(h::_MantidHistogram)
    axes = findall(i -> any(s -> lowercase(strip(s)) == "deltae", (h.names[i], h.units[i])), eachindex(h.names))
    length(axes) == 1 || throw(ArgumentError("expected exactly one DeltaE axis"))
    return only(axes)
end

"""
    axis_vector(name) -> Vector{Float64}

Parse a Mantid HKL direction such as `"[H,H,0]"` or `"[-H,2H,0]"` into three coefficients. The single-axis labels `H`, `K`, and `L` are also supported.
"""
function axis_vector(name::AbstractString)
    label = uppercase(replace(name, r"\s+" => ""))
    label in ("H", "K", "L") && return [label == axis ? 1.0 : 0.0 for axis in ("H", "K", "L")]
    startswith(label, '[') && endswith(label, ']') || throw(ArgumentError("invalid Mantid HKL direction: $name"))
    parts = split(label[2:end-1], ',')
    length(parts) == 3 || throw(ArgumentError("expected three components in HKL direction: $name"))
    return map(parts) do part
        literal = tryparse(Float64, part)
        literal === 0.0 && return 0.0
        term = match(r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)?)\*?[HKL]$", part)
        term === nothing && throw(ArgumentError("invalid component '$part' in HKL direction: $name"))
        coefficient = term.captures[1]
        coefficient in ("", "+") && return 1.0
        coefficient == "-" && return -1.0
        return parse(Float64, coefficient)
    end
end

"""
    hkl_centres(h) -> Matrix{Float64}

Return a `3 × npoints` matrix of stored HKL bin centres, with the energy axis omitted. Momentum bins follow Julia's column-major order on the remaining histogram axes, including integrated momentum axes.
"""
function hkl_centres(h::_MantidHistogram)
    e = energy_axis(h)
    axes = filter(!=(e), eachindex(h.names))
    directions = [axis_vector(h.names[i]) for i in axes]
    return _md_hkl_points(Tuple(h.centres[i] for i in axes), directions)
end

function _md_hkl_points(centres::Tuple, directions)
    points = zeros(3, prod(length, centres; init=1))
    for (j, coords) in enumerate(Iterators.product(centres...))
        for i in eachindex(directions)
            direction = directions[i]
            points[1, j] += direction[1] * coords[i]
            points[2, j] += direction[2] * coords[i]
            points[3, j] += direction[3] * coords[i]
        end
    end
    return points
end
