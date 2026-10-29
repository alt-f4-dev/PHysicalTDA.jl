module PHysicalTDA

#Required Dependencies
using Sunny 
using LinearAlgebra, StaticArrays 
using Ripserer, PersistenceDiagrams 
using Statistics, Base.Threads
using HDF5

#Re-exports
birth = PersistenceDiagrams.birth
death = PersistenceDiagrams.death
PersistenceImage = PersistenceDiagrams.PersistenceImage
export birth, death, PersistenceImage

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

#Topologically-Guided Physical Observables (Naive)
export autosplitbands, naivesplitbands, merge_energy_bands!

include("BandSegmentation.jl") #Physical Bands
export ComponentResult, RidgeNode, RidgeCycleCandidate, RidgeTopologyResult
export adaptive_autosplitbands2d, build_ridge_topology, ridge_dispersions
export components_to_band_marginals, splitbands

include("BandTopology.jl")
export restrict_to_domain
export domain_marginals
export robust_branch_count
export analyze_domain_topologies
export summarize_domain_branches


#Lazy Plotting
# Implemented in ext/PHysicalTDAMakieExt.jl, which Julia loads automatically
# when Makie is loaded (e.g. `using CairoMakie` or `using GLMakie`).
export enable_visuals, plot_persistence_diagram, plot_lifetime_diagram, plot_persistence_entropy, plot_betti_curvature, plot_betti_surface, plot_dataset_overview

const _viz_loaded = Ref(false)
_viz_ext() = Base.get_extension(@__MODULE__, :PHysicalTDAMakieExt)

"""
    enable_visuals(; backend=:gl)

Activate a Makie backend for the `plot_*` functions: `:gl` (GLMakie) or
`:cairo` (CairoMakie). The backend package must be installed in the active
environment; it is loaded if needed. Plotting is provided by a package
extension, so Makie is not loaded by `using PHysicalTDA` alone.
"""
function enable_visuals(; backend::Symbol = :gl)
    _viz_loaded[] && _viz_ext() !== nothing && return nothing
    modsym = backend === :gl    ? :GLMakie :
             backend === :cairo ? :CairoMakie :
             error("Unknown backend: $backend (use :gl or :cairo)")
    @info "Activating $modsym backend..."
    mod = try
        isdefined(Main, modsym) ? getfield(Main, modsym) : Base.require(Main, modsym)
    catch err
        @warn "Failed to load $modsym - $(err)"
        nothing
    end
    mod !== nothing && Base.invokelatest(getfield(mod, :activate!))
    _viz_ext() === nothing && error("PHysicalTDA plotting requires Makie. Add CairoMakie or GLMakie " *
                                    "to your environment and load it, e.g. `using CairoMakie`.")
    _viz_loaded[] = true
    return nothing
end

#Confirm plotting utilities are enabled 
@inline _ensure_viz(; backend::Symbol=:gl) = (_viz_loaded[] || enable_visuals(backend=backend))

#Lazy Plot Entry Point
# The extension may be loaded after the caller's world age, so resolve it
# inside invokelatest.
_viz_call(name::Symbol, args...; kwargs...) = getglobal(_viz_ext(), name)(args...; kwargs...)
function plot_persistence_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_persistence_diagram, args...; kwargs...)
end
function plot_lifetime_diagram(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_lifetime_diagram, args...; kwargs...)
end
function plot_persistence_entropy(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_persistence_entropy, args...; kwargs...)
end
function plot_betti_curvature(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_betti_curvature, args...; kwargs...)
end
function plot_betti_surface(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_betti_surface, args...; kwargs...)
end
function plot_dataset_overview(args...; backend::Symbol=:gl, kwargs...)
    _ensure_viz(backend=backend)
    return Base.invokelatest(_viz_call, :plot_dataset_overview, args...; kwargs...)
end

end #module
