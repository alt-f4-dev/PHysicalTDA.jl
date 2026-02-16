module PHysicalTDA

#Required Dependencies
using Sunny 
using LinearAlgebra, StaticArrays 
using Ripserer, PersistenceDiagrams 
using Statistics, Base.Threads
using HDF5

#Internal Core
include("CoreUtils.jl")
include("DataTypes.jl")
include("DataUtils.jl")

#Result containers
export LSWTResults

#File IO
export saveHDF5, loadHDF5, h5tree
export @name


include("CrystalDB.jl")
#Crystals
export La2CuO4, TMMC 

include("LSWTtools.jl")
#LSWT Tools 
export LSWT, convert4D, collapse

include("TDAPH.jl")
#Persistence Homology
export MAD, pd_array_intensities, pd_sunny_intensities

#Topological Invariants
export persistence_entropy, persistence_entropy_curve, betti_curve, betti_curvature

#Topologically-Guided Physical Observables
export autosplitbands, splitbands, merge_energy_bands!


#Lazy Plotting
export enable_visuals, plot_persistence_diagram, plot_lifetime_diagram, plot_persistence_entropy, plot_betti_curvature, plot_betti_surface, plot_dataset_overview


const _viz_loaded = Ref(false)
function enable_visuals(; backend::Symbol = :gl)
    _viz_loaded[] && return nothing

    function _lazy_import(modsym::Symbol)
        try
            if !isdefined(Main, modsym)
                Base.require(Main, modsym)
            end
            return getfield(Main, modsym)
        catch err
            @warn "Failed to load $modsym - $(err)"
        end
    end

    if backend === :gl
        @info "Activating GLMakie backend..."
        mod = _lazy_import(:GLMakie)
        mod !== nothing && getfield(mod, :activate!)()
    elseif backend === :cairo
        @info "Activating CairoMakie backend..."
        mod = _lazy_import(:CairoMakie)
        mod !== nothing && getfield(mod, :activate!)()
    else
        error("Unknown backend: $backend (use :gl or :cairo)")
    end 
    Base.invokelatest(include, joinpath(@__DIR__,"TopoViz.jl"))
    _viz_loaded[] = true
    return nothing
end

#Confirm plotting utilities are enabled 
@inline _ensure_viz(; backend::Symbol=:gl) = (_viz_loaded[] || enable_visuals(backend=backend))

#Lazy Plot Entry Point
function plot_persistence_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_persistence_diagram, args...; kwargs...)
end
function plot_lifetime_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_lifetime_diagram, args...; kwargs...)
end
function plot_persistence_entropy(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_persistence_entropy, args...; kwargs...)
end
function plot_betti_curvature(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_betti_curvature, args...; kwargs...)
end
function plot_betti_surface(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_betti_surface, args...; kwargs...)
end
function plot_dataset_overview(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(TopoViz.plot_dataset_overview, args...; kwargs...)
end

end #module
