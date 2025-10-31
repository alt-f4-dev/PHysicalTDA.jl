module TopoViz

using Makie #Lazy Plot

#import functions
using ..PHysicalTDA: betti_curve
using Ripserer: birth, death

#export plotting utils
export lifetime_diagram,
       plot_betti_curvature, plot_persistence_entropy, 
       plot_betti_surface, plot_dataset_overview

#-----------------#
# Low-Level Plots #
#-----------------#

function plot_persistence_diagram(PD; label::AbstractString="intensity", 
                                      filtration::AbstractString="superlevel",
                                      ndims::Int=2, maxdim::Int=1, 
                                      title::AbstractString=nothing)
    T = isnothing(title) ? "Persistence Diagram: ($filtration filtration on cubical complex of $(ndims)D array" : title
    fig = Figure(size=(720, 720))
    ax = Axis(fig[1,1], xlabel="birth ($label)", ylabel="death ($label)", title=T)
    for d in 0:maxdim
        bd = PD[d+1]; isempty(bd) && continue
        bth = birth.(bd); dth = death.(bd)
        scatter!(ax, bth, dth; label="Δ$d", marker=:circle)
    end
    Legend(fig[1,2], ax)
    return fig
end

"""
Plot lifetimes (|birth-death|) vs birth from a persistance diagram.
Iterates over homology dimensions [0,1,...,maxdim] with labels.
Assumes input is a Ripserer PD grouped by dimension. Returns Makie `Figure`.
"""
function plot_lifetime_diagram(pd ;maxdim::Int=1)
	living_fig = Figure(size = (720, 720))
    	living_ax = Axis(living_fig[1, 1],
			 xlabel = "Birth (% Intensity)", 
			 ylabel = "Lifetime",
			 title = "Feature Lifetime = |birth - death|")
	label=["⁰","¹","²","³","⁴","⁵","⁶"]
	for d in 0:maxdim
		births = [birth(p) for p in pd[d + 1]]
    		deaths = [death(p) for p in pd[d + 1]]
		living = abs.( (births .- deaths) )
		scatter!(living_ax, births, living; label = "Δ$(label[d+1])", marker = :circle)
	end
	# Add a legend in a separate grid cell
    	Legend(living_fig[1, 2], living_ax)
	return living_fig
end

function plot_betti_curvature(τ::AbstractVector, κ::Dict{Int, <:AbstractVector}; 
                              dims=0:1, title="Betti Curvature")
    fig = Figure(size=(750,520))
    ax = Axis(fig[1,1], xlabel="τ", ylabel="κₚ(τ)", title=title)
    for p in dims
        haskey(κ,p) || continue
        lines!(ax, τ, κ[p], label="κₚ")
    end
    Legend(fig[1,2], ax)
    return fig
end


#---------------------#
# Dataset-Level Plots #
#---------------------#

function plot_persistence_entropy(dataset; p::Int=0, cvals=nothing, title="Entropy-Persistence Map")
    S = [get(r.S, p, 0.0) for r in dataset]
    E = [get(r.E, p, 0.0) for r in dataset]
    C = cvals === nothing ? [r.idx for r in dataset] : cvals
    fig = Figure(size=(680,520))
    ax = Axis(fig[1,1], xlabel="Sₚ", ylabel="Eₚ", title=title)
    sc = scatter!(ax, S, E; color=C, markersize=9)
    Colorbar(fig[1,2], sc; label=cvals===nothing ? "path index" : "color")
    return fig
end


function plot_betti_surface(dataset; p::Int=0, nτ::Int=256, title="Betti Surface βₚ(q,τ)")
    tmins = Float64[]; tmaxs = Float64[]
    for r in dataset
        for Δ in r.PD
            isempty(Δ) && continue
            append!(tmins, filter(isfinite, (birth(x) for x in Δ)))
            append!(tmaxs, filter(isfinite, (death(x) for x in Δ)))
        end
    end
    if isempty(tmins) || isempty(tmaxs)
        error("No finite births/deaths in dataset")
    end
    τmin, τmax = minimum(tmins), maximum(tmaxs)
    τ = range(τmin + eps(Float64), τmax - eps(Float64); length=nτ)
    #β-Matrix: rows = paths, cols = τ
    βmatrix = zeros(Float64, length(dataset), nτ)
    for (i, r) in enumerate(dataset)
        β = betti_curve(r.PD, τ; dims=(p,))
        βmatrix[i,:] .= get(β, p, zeros(Int, nτ))
    end
    fig = Figure(size=(760, 520))
    ax = Axis(fig[1,1], xlabel="τ", ylabel="path index", title=title)
    hm = heatmap!(ax, collect(τ), 1:length(dataset), βmatrix)
    Colorbar(fig[1,2], hm; label="βₚ")
    return fig
end

#------------------#
# Dataset Overview #
#------------------#

# Dataset-level dashboard
function plot_dataset_overview(dataset; p::Int=0, nτ::Int=256, savepath::Union{Nothing,AbstractString}=nothing)
    τinfo = "β surface uses common τ-grid"
    fig = Figure(size=(1200,820))
    # Entropy–Persistence map
    S = [get(r.S,p,0.0) for r in dataset]; E = [get(r.E,p,0.0) for r in dataset]; C = [r.idx for r in dataset]
    ax1 = Axis(fig[1,1], xlabel="S_$p", ylabel="E_$p", title="Entropy–Persistence")
    sc = scatter!(ax1, S, E; color=C, markersize=10)
    Colorbar(fig[1,2], sc; label="path")
    # Betti surface
    # Build common τ and βmat
    tmins = Float64[]; tmaxs = Float64[]
    for r in dataset
        for Δ in r.PD
            isempty(Δ) && continue
            append!(tmins, filter(isfinite, (birth(x) for x in Δ)))
            append!(tmaxs, filter(isfinite, (death(x) for x in Δ)))
        end
    end
    if isempty(tmins) || isempty(tmaxs)
        Label(fig[2,1], "No finite births/deaths in dataset"); savepath===nothing || save(fig, savepath); return fig
    end
    τ = range(minimum(tmins)+eps(), maximum(tmaxs)-eps(); length=nτ)
    βmat = zeros(Float64, length(dataset), nτ)
    for (i, r) in enumerate(dataset)
        β = betti_curve(r.PD, τ; dims=(p,))
        βmat[i, :] .= get(β, p, zeros(Int, nτ))
    end
    ax2 = Axis(fig[2,1], xlabel="τ", ylabel="path", title="Betti Surface β_$p(q,τ)")
    hm = heatmap!(ax2, collect(τ), 1:length(dataset), βmat)
    Colorbar(fig[2,2], hm; label="β_$p")
    Label(fig[2,3], τinfo)
    savepath === nothing || save(fig, savepath)
    return fig
end

end #module
