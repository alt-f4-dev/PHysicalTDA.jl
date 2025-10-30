using Sunny
using .PHysicalTDA            # exports: pd_sunny_intensities, persistence_entropy, betti_curve, betti_curvature
using Ripserer
import PersistenceDiagrams: birth, death
# --- per-path PD + fingerprints ---
struct PathTopo
    idx::Int
    PD                     # Ripserer PD (vector by degree)
    S::Dict{Int,Float64}
    E::Dict{Int,Float64}
    τ::Vector{Float64}
    β::Dict{Int,Vector{Int}}
    κ::Dict{Int,Vector{Float64}}
end
function main()
    # --- τ-grid from finite births/deaths (robust) ---
    function τ_span_from_PD(PD; n::Int=256)
        b = Float64[]; d = Float64[]
        for Δ in PD
            isempty(Δ) && continue
            append!(b, filter(isfinite, birth.(collect(Δ))))
            append!(d, filter(isfinite, death.(collect(Δ))))
        end
        if isempty(b) || isempty(d)
	    return nothing  # <- signal "no finite support"
        end
        tmin, tmax = minimum(b), maximum(d)
        if !(tmax > tmin)
           return nothing
        end
        @assert !isempty(b) "empty PD"
        @assert tmax > tmin "degenerate PD range"
        δ = 1e-12 * max(1.0, abs(tmax-tmin))
        return range(tmin+δ, tmax-δ; length=n)
    end

    # --- finite-lifetime view (aligns with S_p, E_p) ---
    finite_PD(PD) = [isempty(PD[d]) ? PD[d] :
                    PD[d][map(x -> isfinite(birth(x)) && isfinite(death(x)) && death(x)>birth(x), PD[d])]
                    for d in 1:length(PD)]


    function pd_fingerprints_for_paths(cryst, sys; maxdim::Int=1, superlevel::Bool=true, normalize::Bool=true,
                                   dims=0:1, τpts::Int=512)
        result2D, results_paths, energies = LSWT(cryst, sys)
        results, paths = results_paths
        out = PathTopo[]
        for (i, res) in enumerate(results)
	    PD, _ = pd_sunny_intensities(res; maxdim=maxdim, superlevel=superlevel, normalize=normalize)
            S, E  = persistence_entropy(PD; dims=dims, tol=0.0)
            # use finite-only support for τ and β/κ (matches S,E)
            PDf   = finite_PD(PD)
            τ     = τ_span_from_PD(PDf; n=τpts)
            β, κ  = betti_curvature(PDf, τ; dims=dims, scheme=:forward)
            push!(out, PathTopo(i, PD, S, E, collect(τ), β, κ))
        end
        return out, paths, energies
    end

    # --- run ---
    cryst, sys = La2CuO4()
    dataset, paths, energies = pd_fingerprints_for_paths(cryst, sys; maxdim=1, superlevel=true, normalize=true, dims=0:1, τpts=512)

    # minimal report
    for rec in dataset
    @info "path=$(rec.idx)  S₀=$(round(get(rec.S,0,0.0),digits=4))  E₀=$(round(get(rec.E,0,0.0),digits=4))"
    end


    # exact ∫β_p dτ vs E_p (event-sweep)
    function exact_E_from_PD(PDfin, p)::Float64
        D = (p+1 <= length(PDfin)) ? PDfin[p+1] : ()
        isempty(D) && return 0.0
        b   = [birth(x) for x in D]; d = [death(x) for x in D]
        t   = vcat(b,d); s = vcat(ones(Int,length(b)),-ones(Int,length(d)))
        ord = sortperm(t); t = t[ord]; s = s[ord]
        area = 0.0; β = 0
        @inbounds for j in 1:length(t)-1
            β += s[j]; area += β*(t[j+1]-t[j])
        end
        area
    end


    for rec in dataset
        PDf = finite_PD(rec.PD)
        Ep0 = get(rec.E,0,0.0); Eint0 = exact_E_from_PD(PDf,0)
        @info "path=$(rec.idx)  exact∫β₀=$(Eint0)  E₀=$(Ep0)  rel.err=$(Ep0≈0 ? 0.0 : abs(Eint0-Ep0)/max(Ep0,1e-12))"
    end
    
    τaxis = range(0.0, 1.0; length=512)
    κₚ = betti_curvature(finite_PD(dataset[1].PD),τaxis; dims=0:1, scheme=:forward)
    

    plot_persistence_entropy(dataset; p = 0, title="La₂CuO₄: Persistence Entropy")
    plot_betti_surface(dataset; p = 0, nτ = 256, title="Betti Surface: β₀(q,τ)")
    plot_dataset_overview(dataset; p = 0, nτ=256, savepath=nothing)

end
main()
