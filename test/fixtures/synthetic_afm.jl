# Synthetic MDNorm workspaces of a square-lattice Heisenberg antiferromagnet with a known answer. Every
# ingredient of a measured single-crystal data set that the reduction must undo is present: an indexing
# error, nuclear and magnetic Bragg peaks, an incoherent elastic line of known cross section, an aluminium
# flat level and powder lines, a momentum-independent background, and a one-magnon S_mag with the
# L-averaged dipole polarization. Set A is (H, K, E) with a single L bin; set B is coarser and resolves L.
using HDF5, LinearAlgebra, Random, Statistics

const AFM = (a=4.2, c=13.7, J=0.65, Δ=0.005, S=2.5, scale=1.0e3, fAl=1.0, σ_inc=2 * 0.5 + 0.4 + 4 * 0.0008,
             γr0²=0.2905, fwhm=0.25, T=5.0, Ei=10.0, aAl=4.032, bAl=3.449, σincAl=0.0082,
             M=[1.02 -0.005; 0.004 1.025], t=[-0.0025, -0.0013])

afm_γ(h, k) = (cospi(2h) + cospi(2k)) / 2
afm_dispersion(h, k) = 4AFM.J * AFM.S * sqrt(max((1 + AFM.Δ)^2 - afm_γ(h, k)^2, 0.0))
afm_weight(h, k) = AFM.S * sqrt((1 + AFM.Δ - afm_γ(h, k)) / (1 + AFM.Δ + afm_γ(h, k)))
_afm_gauss(x, fwhm) = (σ = fwhm / (2 * sqrt(2 * log(2))); exp(-x^2 / (2σ^2)) / (sqrt(2π) * σ))

"""
    afm_smag(h, k, ω)

True one-magnon `S_mag` per magnetic ion (component sum), resolution-broadened, at true `(h, k)`.
"""
afm_smag(h, k, ω) = ω > 0 ? afm_weight(h, k) * _afm_gauss(ω - afm_dispersion(h, k), AFM.fwhm) : 0.0
afm_background(ω) = 3e-5 + 2e-4 * exp(-(ω - 6)^2 / 4)

"""Scale converting `(D/N)/g` to `S_mag` per ion, `scale/(γ r₀)²`."""
afm_true_scale() = AFM.scale / AFM.γr0²

# Aluminium ring-to-flat ratio for an ideal fcc powder, written out independently of the package.
function _afm_ring_ratio(hkl, m)
    τ = 2π * norm(hkl) / AFM.aAl
    return (2π)^3 * m * 16 * (AFM.bAl^2 / 100) / (4 * AFM.aAl^3 * τ^2 * AFM.σincAl)
end

function _afm_write(path, edges, X; errors2=nothing, metadata=true)
    h5open(path, "w") do file
        entry = create_group(file, "MDHistoWorkspace")
        HDF5.attributes(entry)["SaveMDVersion"] = Int32(2)
        group = create_group(entry, "data")
        group["signal"] = X
        HDF5.attributes(group["signal"])["axes"] = "D3:D2:D1:D0"
        errors2 === nothing || (group["errors_squared"] = errors2)
        for (i, (e, name, unit, frame)) in enumerate(zip(edges, ("[H,0,0]", "[0,K,0]", "DeltaE", "[0,0,L]"),
                                                         ("r.l.u.", "r.l.u.", "DeltaE", "r.l.u."),
                                                         ("HKL", "HKL", "General Frame", "HKL")))
            key = "D$(i - 1)"
            group[key] = e
            HDF5.attributes(group[key])["long_name"] = name
            HDF5.attributes(group[key])["units"] = unit
            HDF5.attributes(group[key])["frame"] = frame
        end
        if metadata
            run = create_group(entry, "experiment0")
            for (key, value) in zip(("a", "b", "c", "alpha", "beta", "gamma"), (AFM.a, AFM.a, AFM.c, 90.0, 90.0, 90.0))
                run["sample/oriented_lattice/unit_cell_$key"] = value
            end
            run["logs/Ei/value"] = [AFM.Ei]
        end
    end
    return path
end

_centres(e) = (e[1:end-1] .+ e[2:end]) ./ 2

"""
    write_synthetic_afm(dir; seed=1) -> (A=(D, N), B=(D, N), truth)

Write sets A and B into `dir` and return their paths with the generating parameters.
"""
function write_synthetic_afm(dir; seed::Integer=1)
    rng = MersenneTwister(seed)
    Minv = inv(AFM.M)
    elastic(ω) = _afm_gauss(ω, AFM.fwhm)
    qnorm(h, k, l) = 2π * sqrt((h^2 + k^2) / AFM.a^2 + (l / AFM.c)^2)
    function polarization(h, k, ls)
        mean(begin
            Q = qnorm(h, k, l)
            F = magnetic_form_factor(Q, MN2_J0)
            F^2 * (1 + (2π * l / AFM.c / Q)^2) / 2
        end for l in ls)
    end

    # Set B: L-resolved, elastic window with aluminium powder lines; its exposure sets the L average.
    hB, ωB, lB = collect(-1.5:0.1:1.5), collect(-3.0:0.25:3.0), collect(-4.0:0.5:4.0)
    hcB, ωcB, lcB = _centres(hB), _centres(ωB), _centres(lB)
    lines = [(2π * sqrt(3) / AFM.aAl, _afm_ring_ratio((1, 1, 1), 8)), (4π / AFM.aAl, _afm_ring_ratio((2, 0, 0), 6))]
    NB = zeros(length(hcB), length(hcB), length(ωcB), length(lcB)); DB = similar(NB)
    for l in eachindex(lcB), j in eachindex(ωcB), k in eachindex(hcB), i in eachindex(hcB)
        q = qnorm(hcB[i], hcB[k], lcB[l])
        e = mean(elastic(ωB[j] + (ωB[j+1] - ωB[j]) * (u - 0.5) / 20) for u in 1:20)
        ring = sum(r * _afm_gauss(q - τ, 2.355 * 0.06) for (τ, r) in lines)
        I = afm_background(ωcB[j]) + AFM.σ_inc / (4π) * e * (1 + AFM.fAl + AFM.fAl * ring) / AFM.scale
        NB[i, k, j, l] = 1e-3
        DB[i, k, j, l] = NB[i, k, j, l] * I * (1 + 0.005randn(rng))
    end

    # Set A: (H, K, E), observed indices h_obs = M h + t, Bragg peaks at true integer and half-integer points.
    hA, ωA, lA = collect(-1.0:0.05:1.0), collect(-3.0:0.1:9.0), [-2.0, 2.0]
    hcA, ωcA = _centres(hA), _centres(ωA)
    shape = (length(hcA), length(hcA), length(ωcA), 1)
    DA, NA = zeros(shape), zeros(shape)
    for k in eachindex(hcA), i in eachindex(hcA)
        ht, kt = Minv * ([hcA[i], hcA[k]] .- AFM.t)
        q = qnorm(hcA[i], hcA[k], 0.0)
        g = polarization(hcA[i], hcA[k], lcB)   # measured momentum: stored indices with the file lattice
        bragg = sum(exp(-((ht - n)^2 + (kt - m)^2) / (2 * 0.025^2)) for n in -2:2, m in -2:2 if (n, m) != (0, 0)) +
                0.3sum(exp(-((ht - n - 0.5)^2 + (kt - m - 0.5)^2) / (2 * 0.025^2)) for n in -2:1, m in -2:1)
        for j in eachindex(ωcA)
            ωcA[j] < AFM.Ei - 1 && q < 3.5 || continue
            ω = ωcA[j]
            I = afm_background(ω) + (AFM.σ_inc / (4π) * elastic(ω) * (1 + AFM.fAl + 1e3bragg) +
                                     AFM.γr0² * g * afm_smag(ht, kt, ω)) / AFM.scale
            NA[i, k, j, 1] = 1e-3 * (1 + 0.3rand(rng))
            DA[i, k, j, 1] = NA[i, k, j, 1] * I
        end
    end
    D2A = (0.01 .* DA .+ 1e-9) .^ 2
    DA .+= sqrt.(D2A) .* randn(rng, size(DA))

    to_file(X) = X   # HDF5.jl writes column-major; "axes" lists D3:…:D0 in C order
    A = (D=_afm_write(joinpath(dir, "afm_A_D.nxs"), [hA, hA, ωA, lA], to_file(DA); errors2=to_file(D2A)),
         N=_afm_write(joinpath(dir, "afm_A_N.nxs"), [hA, hA, ωA, lA], to_file(NA)))
    B = (D=_afm_write(joinpath(dir, "afm_B_D.nxs"), [hB, hB, ωB, lB], to_file(DB); errors2=to_file((0.005 .* DB) .^ 2)),
         N=_afm_write(joinpath(dir, "afm_B_N.nxs"), [hB, hB, ωB, lB], to_file(NB)))
    return (; A, B, truth=AFM)
end
