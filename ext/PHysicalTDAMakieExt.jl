module PHysicalTDAMakieExt
# Plotting for PHysicalTDA, loaded automatically when Makie is loaded.
# Called through the lazy wrappers `PHysicalTDA.plot_*`.

using Makie
using PHysicalTDA: betti_curve, birth, death
using PHysicalTDA: QFIResult, BandQFIResult, ChildRidgeResult, RidgeCandidates, qfi_uncertainty, qfi_bound,
                   band_supports


#export plotting utils (mutating)

#-----------------#
# Low-Level Plots #
#-----------------#
function plot_persistence_diagram(PD; label::AbstractString="intensity", 
                                      filtration::AbstractString="superlevel",
                                      ndims::Int=2, maxdim::Int=1, 
                                      title::Union{AbstractString,Nothing}=nothing)
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
            append!(tmins, filter(isfinite, [birth(x) for x in Δ]))
            append!(tmaxs, filter(isfinite, [death(x) for x in Δ]))
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
            append!(tmins, filter(isfinite, [birth(x) for x in Δ]))
            append!(tmaxs, filter(isfinite, [death(x) for x in Δ]))
        end
    end
    if isempty(tmins) || isempty(tmaxs)
        Label(fig[2,1], "No finite births/deaths in dataset"); savepath===nothing || save(savepath, fig); return fig
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
    savepath === nothing || save(savepath, fig)
    return fig
end


#--------------------#
# Analysis Layouts   #
#--------------------#
_letter(k) = "(" * Char('a' + k - 1) * ")"
_colour(k) = Makie.wong_colors()[mod1(k, 7)]

function _decorate!(ax, markers, shade)
    for (lo, hi) in shade
        vspan!(ax, lo, hi; color=(:gray, 0.15))
    end
    for x in markers
        vlines!(ax, x; color=:gray, linestyle=:dash, linewidth=1)
    end
end

function plot_qfi_summary(r::QFIResult, qs; normalization, m::Integer=1, nσ::Real=2, scale_σ::Real=0,
                          systematic::Bool=true, bands=nothing, sensitivity=nothing, comparison=nothing,
                          models=NamedTuple(), markers=(), shade=(), xlabel="q", band_labels=nothing,
                          title=nothing, size=nothing)
    length(qs) == length(r.FQ) || throw(DimensionMismatch("qs must have one entry per momentum of r"))
    q = collect(Float64, qs)
    B = qfi_bound(; normalization...)
    σtot = qfi_uncertainty(r; scale_σ, systematic)
    lb = (r.FQ .- nσ .* σtot) ./ B
    panels = [:FQ, :nqfi]
    bands === nothing || push!(panels, :ξ)
    sensitivity === nothing || push!(panels, :sensitivity)
    comparison === nothing || push!(panels, :comparison)
    ncol = length(panels) <= 2 ? length(panels) : 2
    nrow = cld(length(panels), ncol)
    fig = Figure(size=size === nothing ? (520ncol, 380nrow) : size)
    title === nothing || Label(fig[0, 1:ncol], title; fontsize=18, font=:bold)
    panes = Dict{Symbol,Any}()
    for (k, panel) in enumerate(panels)
        pos = fig[cld(k, ncol), mod1(k, ncol)]
        if panel === :FQ
            ax = Axis(pos; xlabel, ylabel="F_Q", title=_letter(k) * " F_Q(q)")
            _decorate!(ax, markers, shade)
            band!(ax, q, r.FQ .- σtot, r.FQ .+ σtot; color=(:steelblue, 0.3))
            lines!(ax, q, r.FQ; color=:steelblue, label="measured")
            for (name, curve) in pairs(models)
                lines!(ax, q, collect(Float64, curve); linestyle=:dash, label=string(name))
            end
            axislegend(ax; position=:lt, framevisible=false)
        elseif panel === :nqfi
            ax = Axis(pos; xlabel, ylabel="nQFI lower bound", title=_letter(k) * " nQFI_lb(q), $(nσ)σ")
            _decorate!(ax, markers, shade)
            hlines!(ax, m; color=:black, linestyle=:dot)
            lines!(ax, q, lb; color=:black)
            cert = lb .> m
            any(cert) && scatter!(ax, q[cert], lb[cert]; color=:crimson, markersize=7, label="certified (m = $m)")
            any(cert) && axislegend(ax; position=:lt, framevisible=false)
        elseif panel === :ξ
            ax = Axis(pos; xlabel, ylabel="ξ", title=_letter(k) * " band fractions")
            _decorate!(ax, markers, shade)
            for b in axes(bands.ξ, 1)
                name = band_labels === nothing ? "band $b" : band_labels[b]
                ok = isfinite.(bands.ξ[b, :]) .& isfinite.(bands.σξ[b, :])
                band!(ax, q[ok], bands.ξ[b, ok] .- bands.σξ[b, ok], bands.ξ[b, ok] .+ bands.σξ[b, ok]; color=(_colour(b), 0.25))
                lines!(ax, q, bands.ξ[b, :]; color=_colour(b), label=name)
            end
            axislegend(ax; position=:lb, framevisible=false)
        elseif panel === :sensitivity
            ax = Axis(pos; xlabel, ylabel="nQFI lower bound", title=_letter(k) * " lower energy limit")
            _decorate!(ax, markers, shade)
            hlines!(ax, m; color=:black, linestyle=:dot)
            for (j, e) in enumerate(sensitivity.lower_limits)
                lines!(ax, q, sensitivity.nqfi_lb[j, :]; label="ω > $(round(e; digits=3))")
            end
            axislegend(ax; position=:lt, framevisible=false)
        else
            ax = Axis(pos; xlabel, ylabel="Δξ", title=_letter(k) * " support comparison")
            _decorate!(ax, markers, shade)
            hlines!(ax, 0; color=:black, linestyle=:dot)
            for b in axes(comparison.Δξ, 1)
                name = band_labels === nothing ? "band $b" : band_labels[b]
                ok = isfinite.(comparison.Δξ[b, :]) .& isfinite.(comparison.σΔξ[b, :])
                errorbars!(ax, q[ok], comparison.Δξ[b, ok], nσ .* comparison.σΔξ[b, ok]; color=_colour(b))
                scatter!(ax, q[ok], comparison.Δξ[b, ok]; color=_colour(b), markersize=6, label=name)
            end
            axislegend(ax; position=:lt, framevisible=false)
        end
        panes[panel] = ax
    end
    return fig, NamedTuple{Tuple(panels)}(Tuple(panes[p] for p in panels))
end

function _overlay_masks!(ax, q, ω, masks)
    labels = fill(NaN, length(q), length(ω))
    for (b, M) in enumerate(masks)
        labels[M] .= b
    end
    any(isfinite, labels) || return
    heatmap!(ax, q, ω, labels; colormap=Makie.Categorical(:tab10), alpha=0.45)
end

function plot_ridge_workflow(I::AbstractMatrix{<:Real}, qs, ωs; valid=nothing, segmentation=nothing,
                             candidates=nothing, children=nothing, supports=nothing, trajectories=nothing,
                             markers=(), shade=(), xlabel="q", ylabel="ω (meV)", colormap=:viridis,
                             colorrange=nothing, title=nothing, size=nothing)
    q, ω = collect(Float64, qs), collect(Float64, ωs)
    Base.size(I) == (length(q), length(ω)) || throw(DimensionMismatch("I must be length(qs) × length(ωs)"))
    X = Float64.(I)
    valid === nothing || (X[.!valid] .= NaN)
    range = colorrange === nothing ? (0.0, max(eps(), quantile_finite(X, 0.99))) : colorrange
    panels = [:data]
    candidates === nothing || push!(panels, :candidates)
    segmentation === nothing || push!(panels, :parents)
    children === nothing || push!(panels, :children)
    supports === nothing || push!(panels, :supports)
    trajectories === nothing || push!(panels, :trajectories)
    ncol = min(3, length(panels))
    nrow = cld(length(panels), ncol)
    fig = Figure(size=size === nothing ? (430ncol + 80, 360nrow) : size)
    title === nothing || Label(fig[0, 1:ncol], title; fontsize=18, font=:bold)
    panes = Dict{Symbol,Any}()
    names = Dict(:data => "intensity", :candidates => "ridge candidates", :parents => "parent supports",
                 :children => "child ridges", :supports => "refined supports", :trajectories => "trajectories")
    hm = nothing
    for (k, panel) in enumerate(panels)
        ax = Axis(fig[cld(k, ncol), mod1(k, ncol)]; xlabel, ylabel, title=_letter(k) * " " * names[panel])
        hm = heatmap!(ax, q, ω, X; colormap, colorrange=range, nan_color=:white)
        _decorate!(ax, markers, shade)
        if panel === :candidates
            c = candidates isa ChildRidgeResult ? candidates.candidates : candidates
            for gp in c.groups
                lines!(ax, [c.q[gp.iq], c.q[gp.iq]], [c.ω[gp.lo], c.ω[gp.hi]]; color=(:white, 0.6), linewidth=1)
            end
            scatter!(ax, [n.q for n in c.nodes], [n.ω for n in c.nodes]; color=:red, markersize=4)
        elseif panel === :parents
            masks = segmentation isa AbstractVector ? segmentation : band_supports(segmentation, length(q), length(ω))
            _overlay_masks!(ax, q, ω, masks)
        elseif panel === :children
            for (ip, T) in enumerate(children.tracks), b in axes(T, 1)
                lines!(ax, q, Float64.(T[b, :]); color=_colour(ip), linewidth=2)
            end
            nodes = children.graph_nodes
            isempty(children.split_nodes) || scatter!(ax, [nodes[i].q for i in children.split_nodes],
                                                      [nodes[i].ω for i in children.split_nodes]; color=:white, marker=:utriangle)
            isempty(children.merge_nodes) || scatter!(ax, [nodes[i].q for i in children.merge_nodes],
                                                      [nodes[i].ω for i in children.merge_nodes]; color=:white, marker=:dtriangle)
        elseif panel === :supports
            _overlay_masks!(ax, q, ω, supports)
        elseif panel === :trajectories
            for b in axes(trajectories.ω, 1)
                lines!(ax, trajectories.q, trajectories.ω[b, :]; color=_colour(b), linewidth=2)
                imp = trajectories.imputed[b, :] .& isfinite.(trajectories.ω[b, :])
                any(imp) && scatter!(ax, trajectories.q[imp], trajectories.ω[b, imp]; color=:white, strokecolor=_colour(b),
                                     strokewidth=1, markersize=6)
            end
        end
        panes[panel] = ax
    end
    Colorbar(fig[1:nrow, ncol+1], hm; label="intensity")
    return fig, NamedTuple{Tuple(panels)}(Tuple(panes[p] for p in panels))
end

function quantile_finite(X, p)
    v = sort!(filter(isfinite, vec(X)))
    isempty(v) && return 1.0
    return v[clamp(ceil(Int, p * length(v)), 1, length(v))]
end

end #module
