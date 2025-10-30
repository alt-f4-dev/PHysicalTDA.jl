module PHysicalTDA

#Required Dependencies
using Sunny, GLMakie, LinearAlgebra, StaticArrays, Ripserer, PersistenceDiagrams, Statistics, Base.Threads

#Public API
export La2CuO4, LSWT, convert4D, collapse, 
       MAD, lifetime_diagram, pd_array_intensity, pd_sunny_intensities, persistence_entropy, betti_curvature,
       plot_persistence_entropy, plot_betti_curvature, plot_betti_surface, plot_dataset_overview

#Internal Submodules 
include("CrystalDB.jl"); include("LSWTtools.jl"); include("TDAPH.jl"); include("TopoViz.jl")

end #module
