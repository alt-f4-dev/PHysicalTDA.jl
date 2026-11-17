#DataUtils.jl (NOT a submodule)

#Generic HDF5 save (dispatches by type)
function saveHDF5 end
function loadHDF5 end

#---------------------#
# LSWTResults - HDF5  #
#---------------------#
"""
    saveHDF5(crystal, sweepParam, sweepValue, qnames, results::LSWTResults; fname, root)

Save LSWT results to HDF5 using the canonical hierarchy:

    root/crystal/LSWT/sweepParam/sweepValue/direction/{qpoints, energies, intensities}

This enforces a structured scientific database layout where:
- `crystal` identifies the material,
- `LSWT` identifies the computational method (determined by dispatch),
- `sweepParam` defines the swept physical parameter (e.g. "Temp"),
- `sweepValue` defines the specific parameter value (e.g. "0K"),
- `direction` corresponds to q-space paths provided by `qnames`.

# Arguments
- `crystal::String`  
  Crystal or material identifier (e.g., "La2CuO4").

- `sweepParam::String`  
  Name of the swept parameter (e.g., "Temp", "Field", "Pressure").

- `sweepValue::String`  
  Label for the specific parameter value (e.g., "0K", "10K", "B=5T").

- `qnames::Vector{String}`  
  Direction labels corresponding to q-space paths (typically from the `@name` macro).

- `results::LSWTResults`  
  LSWT computation results container.

# Keyword Arguments
- `fname::String="crystal-data"`  
  Output filename (without `.h5` extension).

- `root::String="crystals"`  
  Top-level group name in the HDF5 hierarchy.

# Notes
- Data is appended to the file, allowing multiple sweeps and methods to coexist.
- Metadata stored in `results.meta` is written to:

      root/crystal/LSWT/sweepParam/sweepValue/meta/

- Physical state variables (e.g., temperature, field, doping) should be encoded via `sweepParam` and `sweepValue`, not stored as metadata.
"""
function saveHDF5(crystal::String, sweepParam::String, sweepValue::String,
                  qnames::Vector{String}, results::LSWTResults; 
                  fname::String="crystal-data", root::String="crystals")
    method = "LSWT"
    h5open(fname*".h5", "cw") do h5
        basepath = "$root/$crystal/$method/$sweepParam/$sweepValue"
        for (idx, name) in enumerate(qnames)
            path = "$basepath/$name"
            qpoints = reduce(hcat, results.qpaths[idx].qs)
            h5["$path/qpoints"] = qpoints
            h5["$path/energies"] = collect(results.energies)
            h5["$path/intensities"] = results.intensities[idx].data
        end
        if !isempty(results.meta)
            metapath = "$basepath/meta"
            for (key,val) in results.meta
                try 
                    h5["$metapath/$key"] = val
                catch e
                    @warn "Could not save metadata key `$key` : $e"
                end
            end
        end
    end
    @info "Saved LSWT results to $(fname).h5"
end
"""
    loadHDF5(fname; keypath=nothing) -> data

Load data from an HDF5 file.

# Arguments
- `fname::String`: Filename without .h5 extension
- `keypath::Union{String,Nothing}`: HDF5 internal path to specific dataset or group.
  If `nothing` (default), loads entire file structure as nested Dict.

# Returns
- `Dict{String,Any}`: Nested dictionary mirroring HDF5 hierarchy, or the specific dataset/group if `keypath` is provided

# Examples
```julia
# Load entire file structure
data = loadHDF5("crystal-data")
data["crystals"]["La2CuO4"]["LSWT"]["Temp"]["0K"]["H"]["intensities"]  # Access nested data

#Load specific group
LCO_data = loadHDF5("crystal-data"; keypath="crystals/La2CuO4/LSWT")
LCO_data["Temp"]["0K"]["H"]["intensities"]  # Only La2CuO4 LSWT data loaded

#Load specific direction
H_data = loadHDF5("crystal-data"; keypath="crystals/La2CuO4/LSWT/Temp/0K/H")
H_data["intensities"]  # Only H-direction data

#Load specific dataset (array)
intensities = loadHDF5("crystal-data"; keypath="crystals/La2CuO4/LSWT/Temp/0K/H/intensities")
#intensities is now a Matrix{Float64}

#Load q-points
qpoints = loadHDF5("crystal-data"; keypath="crystals/La2CuO4/LSWT/Temp/0K/H/qpoints")
```

# Notes
- HDF5 hierarchy uses "/" as path separator (like filesystem paths)
- Use `keypath` to load only what you need (faster for large files)
- Returned structure depends on what's stored at `keypath`:
  - Group → Dict with nested structure
  - Dataset → Array or scalar value
"""
function loadHDF5(fname::String; keypath::Union{String,Nothing}=nothing)
    h5open(fname*".h5","r") do h5 
        if keypath == nothing; return read(h5); else; return read(h5[keypath]); end
    end
end

#----------------#
# HDF5 Structure #
#----------------#
function _with_h5_object(f, parent, name::AbstractString)
    obj = parent[name]
    try
        return f(obj)
    finally
        close(obj)
    end
end

function _tree_preview(value, maxstring::Int)
    if value isa AbstractString
        prefix = first(value, maxstring)
        return repr(length(prefix) == length(value) ? prefix : prefix * "…")
    elseif value isa AbstractArray
        return "[" * join((_tree_preview(item, maxstring) for item in value), ", ") * "]"
    end
    return repr(value)
end

function _tree_has_preview(item)
    dims = size(item)
    dims === nothing && return false
    # Cap previews by element count and fixed-width storage size before doing any data I/O.
    count = prod(dims)
    0 < count <= 8 || return false
    dtype = HDF5.datatype(item)
    try
        return sizeof(dtype) <= 65_536 ÷ count
    finally
        close(dtype)
    end
end

function _print_tree_attributes(node, prefix, options)
    options.attributes || return
    attrs = HDF5.attributes(node)
    names = keys(attrs)
    for name in Iterators.take(names, options.maxchildren)
        item = attrs[name]
        try
            preview = _tree_has_preview(item) ? _tree_preview(read(item), options.maxstring) : "(preview omitted)"
            println(options.io, prefix, "@", name, " = ", preview)
        finally
            close(item)
        end
    end
    length(names) > options.maxchildren && println(options.io, prefix, "… ", length(names) - options.maxchildren,
                                                  " more attributes")
end

function _print_tree_dataset(item, options)
    dims = size(item)
    dims_str = dims === nothing ? "empty" : isempty(dims) ? "scalar" : join(dims, "×")
    print(options.io, "  [", eltype(item), ": ", dims_str, "]")
    if options.values && _tree_has_preview(item)
        print(options.io, " = ", _tree_preview(read(item), options.maxstring))
    end
    println(options.io)
end

function _print_tree_item(item, child_prefix, depth, options)
    if item isa HDF5.Group
        println(options.io, "/")
        _print_tree_attributes(item, child_prefix, options)
        _print_tree(item, child_prefix, depth, options)
    elseif item isa HDF5.Dataset
        _print_tree_dataset(item, options)
        _print_tree_attributes(item, child_prefix, options)
    else
        println(options.io, "  [", typeof(item), "]")
    end
end

function _print_tree(node, prefix, depth, options)
    depth >= options.maxdepth && return
    names = keys(node)
    for (i, name) in enumerate(Iterators.take(names, options.maxchildren))
        last = i == length(names)
        print(options.io, prefix, last ? "└── " : "├── ", name)
        child_prefix = prefix * (last ? "    " : "│   ")
        _with_h5_object(node, name) do item
            _print_tree_item(item, child_prefix, depth + 1, options)
        end
    end
    length(names) > options.maxchildren && println(options.io, prefix, "└── … ", length(names) - options.maxchildren,
                                                  " more children")
end

"""
    h5tree(fname; keypath=nothing, maxdepth=typemax(Int), attributes=false, values=false, maxchildren=typemax(Int), maxstring=200, io=stdout)

Print an HDF5 or NeXus tree to `io`. Open an existing path directly and otherwise try the legacy `.h5` suffix. Set `keypath` to inspect one group or dataset.

Enable `attributes` to display attributes and `values` to preview scalar and short datasets. `maxchildren` limits the children and attributes displayed per object, and `maxstring` limits each string preview in characters. Large fixed-width values are omitted.
"""
function h5tree(fname::AbstractString; keypath::Union{Nothing,AbstractString}=nothing, maxdepth::Int=typemax(Int),
                attributes::Bool=false, values::Bool=false, maxchildren::Int=typemax(Int), maxstring::Int=200, io::IO=stdout)
    min(maxdepth, maxchildren, maxstring) >= 0 || throw(ArgumentError("tree limits must be nonnegative"))
    path = isfile(fname) ? fname : fname * ".h5"
    options = (; maxdepth, attributes, values, maxchildren, maxstring, io)
    h5open(path, "r") do file
        if keypath === nothing || keypath == "/"
            println(io, path)
            _print_tree_attributes(file, "", options)
            _print_tree(file, "", 0, options)
        else
            clean_path = rstrip(keypath, '/')
            parts = rsplit(clean_path, '/'; limit=2)
            name = last(parts)
            parent_path = length(parts) == 2 ? first(parts) : ""
            isempty(parent_path) || println(io, parent_path, "/")
            print(io, isempty(parent_path) ? "" : "└── ", name)
            _with_h5_object(file, keypath) do item
                _print_tree_item(item, "    ", 0, options)
            end
        end
    end
    return nothing
end
