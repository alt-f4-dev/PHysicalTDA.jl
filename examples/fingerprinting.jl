# Topological fingerprints (persistence entropy, total persistence, Betti
# curves) for LSWT spectra of La₂CuO₄ along three q-paths.
#
# Run from an environment containing PHysicalTDA, Sunny and CairoMakie:
#     julia --project examples/fingerprinting.jl
using PHysicalTDA, Sunny, CairoMakie

struct PathTopo
    idx::Int
    PD                       # persistence diagrams, grouped by degree
    S::Dict{Int,Float64}     # persistence entropy per degree
    E::Dict{Int,Float64}     # total persistence per degree
    τ::Vector{Float64}
    β::Dict{Int,Vector{Int}}
    κ::Dict{Int,Vector{Float64}}
end

# Keep only finite, non-degenerate bars (matches the lifetimes used by S, E).
finite_PD(PD) = [isempty(D) ? D :
                 D[map(x -> isfinite(birth(x)) && isfinite(death(x)) && death(x) > birth(x), D)]
                 for D in PD]

# Uniform τ grid spanning the finite births and deaths.
function τ_span(PD; n::Int=256)
    b = Float64[]; d = Float64[]
    for D in PD
        isempty(D) && continue
        append!(b, birth.(D)); append!(d, death.(D))
    end
    (isempty(b) || isempty(d)) && return nothing
    tmin, tmax = minimum(b), maximum(d)
    tmax > tmin || return nothing
    δ = 1e-12 * max(1.0, tmax - tmin)
    return collect(range(tmin + δ, tmax - δ; length=n))
end

# Exact ∫β_p dτ by an event sweep; equals total persistence E_p.
function exact_E(PDfin, p)
    D = p + 1 <= length(PDfin) ? PDfin[p+1] : ()
    isempty(D) && return 0.0
    t = vcat(birth.(D), death.(D)); s = vcat(ones(Int, length(D)), -ones(Int, length(D)))
    ord = sortperm(t); t = t[ord]; s = s[ord]
    area = 0.0; β = 0
    for j in 1:length(t)-1
        β += s[j]; area += β * (t[j+1] - t[j])
    end
    return area
end

function path_fingerprints(cryst, sys, paths; dims=0:1, τpts::Int=512)
    lswt = LSWT(cryst, sys, paths)
    out = PathTopo[]
    for (i, res) in enumerate(lswt.intensities)
        PD = pd_sunny_intensities(res; maxdim=maximum(dims), superlevel=true, normalize=true)
        S, E = persistence_entropy(PD; dims=dims)
        PDf = finite_PD(PD)
        τ = τ_span(PDf; n=τpts)
        τ === nothing && continue
        β, κ = betti_curvature(PDf, τ; dims=dims, scheme=:forward)
        push!(out, PathTopo(i, PD, S, E, τ, β, κ))
    end
    return out
end

function main(; outdir=mktempdir())
    cryst, sys = La2CuO4()
    paths = [q_space_path(cryst, [[-1, 0, 0], [0, 0, 0], [1, 0, 0]], 400),
             q_space_path(cryst, [[0, -1, 0], [0, 0, 0], [0, 1, 0]], 400),
             q_space_path(cryst, [[0, 0, -1], [0, 0, 0], [0, 0, 1]], 400)]
    dataset = path_fingerprints(cryst, sys, paths)

    for rec in dataset
        Eint = exact_E(finite_PD(rec.PD), 0)
        @info "path=$(rec.idx)  S₀=$(round(rec.S[0]; digits=4))  E₀=$(round(rec.E[0]; digits=4))  ∫β₀dτ=$(round(Eint; digits=4))"
    end

    f1 = plot_persistence_entropy(dataset; p=0, title="La₂CuO₄: persistence entropy", backend=:cairo)
    f2 = plot_betti_surface(dataset; p=0, nτ=256, title="Betti surface β₀(path, τ)", backend=:cairo)
    f3 = plot_dataset_overview(dataset; p=0, nτ=256, savepath=joinpath(outdir, "overview.png"), backend=:cairo)
    save(joinpath(outdir, "entropy.png"), f1)
    save(joinpath(outdir, "betti_surface.png"), f2)
    @info "figures written to $outdir"
    return dataset
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
