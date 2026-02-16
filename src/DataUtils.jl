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
- Data is appended to the file (`"a"` mode), allowing multiple sweeps and methods to coexist.
- Metadata stored in `results.meta` is written to:

      root/crystal/LSWT/sweepParam/sweepValue/meta/

- Physical state variables (e.g., temperature, field, doping) should be encoded via `sweepParam` and `sweepValue`, not stored as metadata.
"""
function saveHDF5(crystal::String, sweepParam::String, sweepValue::String,
                  qnames::Vector{String}, results::LSWTResults; 
                  fname::String="crystal-data", root::String="crystals")
    method = "LSWT"
    h5open(fname*".h5", "a") do h5
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
function _print_tree(node, prefix::String, is_last::Bool, depth::Int, maxdepth::Int)
    depth >= maxdepth && return
    
    # Get all keys
    keys_list = keys(node)
    n = length(keys_list)
    
    for (i, key) in enumerate(keys_list)
        is_last_item = (i == n)
        
        # Print current item
        branch = is_last_item ? "└── " : "├── "
        item = node[key]
        
        if isa(item, HDF5.Group)
            println(prefix * branch * key * "/")
            # Recurse into group
            new_prefix = prefix * (is_last_item ? "    " : "│   ")
            _print_tree(item, new_prefix, is_last_item, depth + 1, maxdepth)
        else
            # Dataset - show type and shape
            dims = size(item)
            dtype = eltype(item)
            dims_str = isempty(dims) ? "scalar" : join(dims, "×")
            println(prefix * branch * key * "  [$dtype: $dims_str]")
        end
    end
end
"""
    h5tree(fname; keypath=nothing, maxdepth=Inf)

Display HDF5 file structure in tree format.

# Arguments
- `fname::String`: HDF5 filename (without .h5 extension)
- `keypath::Union{String,Nothing}`: Start from specific path (default: root)
- `maxdepth::Int`: Maximum depth to display (default: Inf = all levels)

# Examples
```julia
#Show entire file
h5tree("crystal-data")

#Show specific subtree
h5tree("crystal-data"; keypath="crystals/La2CuO4")

#Limit depth
h5tree("crystal-data"; keypath="crystals/La2CuO4/LSWT", maxdepth=1)
```
"""
function h5tree(fname::String; keypath::Union{String,Nothing}=nothing, maxdepth::Int=typemax(Int))
    h5open(fname * ".h5", "r") do h5
        root = keypath === nothing ? h5 : h5[keypath]
        
        if keypath === nothing
            # At root - show filename as root
            println(fname * ".h5")
            _print_tree(root, "", true, 0, maxdepth)
        else
            # Strip trailing slashes
            clean_path = rstrip(keypath, '/')
            
            # At keypath - show as root with connecting branch
            parts = split(clean_path, '/')
            parent_path = join(parts[1:end-1], '/')
            node_name = parts[end]
            
            if !isempty(parent_path)
                println(parent_path * "/")
                println("└── " * node_name * "/")
                _print_tree(root, "    ", true, 0, maxdepth)
            else
                println(node_name * "/")
                _print_tree(root, "    ", true, 0, maxdepth)
            end
        end
    end
end
