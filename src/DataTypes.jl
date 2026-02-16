#DataTypes.jl (NOT a submodule)
abstract type PHysicalTDAResult end

"""
Result container for Linear Spin Wave Theory calculations obtained from Sunny.

Stores broadened intensity data along multiple q-space paths,
compatible with Sunny.jl output format.

# Fields
- `intensities::Vector{Sunny.Intensities{T, Sunny.QPath, 2}}`: Broadened intensity objects
- `qpaths::Vector{Sunny.QPath}`: Q-space paths
- `energies::E`: Energy grid (typically `StepRangeLen` or similar)
- `meta::Dict{String,Any}`: Metadata (temperature, parameters, etc.)

# Constructors
```julia
# From LSWT() output:
LSWTResults(nt::NamedTuple, meta=Dict())

# Full construction:
LSWTResults(intensities, qpaths, energies, meta)
```
"""
struct LSWTResults{T<:Real, E<:AbstractVector{T}} <: PHysicalTDAResult
    intensities::Vector{Sunny.Intensities{T, Sunny.QPath, 2}}
    qpaths::Vector{Sunny.QPath}
    energies::E
    meta::Dict{String,Any}
end
function LSWTResults(intensities, qpaths, energies; meta=Dict{String,Any}())
    T = eltype(energies)
    LSWTResults{T,typeof(energies)}(intensities, qpaths, energies, meta)
end
LSWTResults(nt::NamedTuple; meta::Dict{String,Any}=Dict{String,Any}()) = LSWTResults(nt.intensities, 
                                                                                     nt.qpaths, 
                                                                                     nt.energies; 
                                                                                     meta=meta)
