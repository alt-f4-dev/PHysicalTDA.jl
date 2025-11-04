module PHysicalTDA

#Required Dependencies
using Sunny, LinearAlgebra, StaticArrays, Ripserer, PersistenceDiagrams, Statistics, Base.Threads

#Public API
export La2CuO4, LSWT, convert4D, collapse, 
       MAD, pd_array_intensity, pd_sunny_intensities, 
       persistence_entropy, betti_curvature

#Internal Submodules 
include("CrystalDB.jl"); include("LSWTtools.jl"); include("TDAPH.jl")

#Lazy Plotting
export enable_visuals, 
       plot_persistence_diagram, lifetime_diagram,
       plot_persistence_entropy, plot_betti_curvature, 
       plot_betti_surface, plot_dataset_overview
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
        #@eval begin
            #import GLMakie
            #GLMakie.activate!()
            #mod = Base.require(Base, :GLMakie)
            #getfield(mod, :activate!)()
        #end
        #Core.eval(@__MODULE__, :(import GLMakie))
        #GLMakie.activate!()
        mod = _lazy_import(:GLMakie)
        mod !=== nothing && getfield(mod, :activate!)()
    elseif backend === :cairo
        @info "Activating CairoMakie backend..."
        #@eval begin
            #import CairoMakie
            #CairoMakie.activate!()
            #mod = Base.require(Base, :CairoMakie)
            #getfield(mod, :activate!)()
        #end
        #Core.eval(@__MODULE__, :(import CairoMakie))
        #CairoMakie.activate!()
        mod = _lazy_import(:CairoMakie)
        mod !=== nothing && getfield(mod, :activate!)()
    else
        error("Unknown backend: $backend (use :gl or :cairo)")
    end 
    @eval include(joinpath(@__DIR__, "TopoViz.jl"))
    _viz_loaded[] = true
    return nothing
end

@inline _ensure_viz(; backend::Symbol=:gl) = (_viz_loaded[] || enable_visuals(backend=backend))

#Lazy Plot Entry
function plot_persistence_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoViz.plot_persistence_diagram(args...; kwargs...)
end
function lifetime_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoViz.lifetime_diagram(args...; kwargs...)
end
function plot_persistence_entropy(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoViz.plot_persistence_entropy(args...; kwargs...)
end
function plot_betti_curvature(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoViz.plot_betti_curvature(args...; kwargs...)
end
function plot_betti_surface(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoViz.plot_betti_curvature(args...; kwargs...)
end
function plot_data_set_overview(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return TopoVize.plot_betti_curvature(args...; kwargs...)
end

end #module
