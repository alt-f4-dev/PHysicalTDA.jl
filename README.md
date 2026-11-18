<p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/dark-mode-banner.svg">
      <img src="docs/src/assets/light-mode-banner.svg" alt="PHysicalTDA.jl" width="600"/>
    </picture>
</p>


[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alt-f4-dev.github.io/PHysicalTDA.jl/dev/)
[![Build Status](https://github.com/alt-f4-dev/PHysicalTDA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alt-f4-dev/PHysicalTDA.jl/actions/workflows/CI.yml?query=branch%3Amain)

`PHysicalTDA.jl` performs Topological Data Analysis (TDA) by utilizing Persistent Homology (PH) on scattering datasets.
It is written to handle both synthetic and raw datasets. The module currently supports Julia native and `Sunny.jl` types.
Mantid `SaveMD` v2 histogram files (`.nxs`) and MDNorm workspace pairs can be read directly.

## Installation

`PHysicalTDA.jl` is not yet registered. It supports Julia 1.10 (LTS) and every later 1.x release.

```julia
using Pkg
Pkg.add(url="https://github.com/alt-f4-dev/PHysicalTDA.jl")
```

## Quick start

```julia
using PHysicalTDA

qs = range(0, 1; length=128)                  # momentum transfer (r.l.u.)
ωs = range(0, 10; length=256)                 # energy transfer (meV)
Iqω = [exp(-(ω - 4 - 2sinpi(q))^2 / 0.3) + 0.6exp(-(ω - 8)^2 / 0.3) for q in qs, ω in ωs]

# Static fingerprint: superlevel persistence of the full intensity map
PD   = pd_array_intensities(Iqω; maxdim=1, superlevel=true, normalize=true)
S, E = persistence_entropy(PD; dims=0:1)                    # entropy and total persistence per degree
β    = betti_curve(PD, range(0, 1; length=64); dims=0:1)    # β₀(τ), β₁(τ)

# Dynamic fingerprint (naive): partition-resolved marginals of the 2D spectral weight
marg = partition_marginals(Iqω, qs, ωs)
marg.splitE.bandIq                                          # I(q) for each energy partition

# Dynamic fingerprint (physical): parent bands, their marginals, and ridge topology
seg = adaptive_autosplitbands2d(Iqω, qs, ωs)
seg.qfirst.nparents, seg.ridge_beta0, seg.ridge_beta1
dm  = domain_marginals(Iqω, seg, qs, ωs)
```

## Overview

| Task | Functions |
|---|---|
| Persistence diagrams | `pd_array_intensities`, `pd_sunny_intensities` |
| Topological invariants | `persistence_entropy`, `persistence_entropy_curve`, `betti_curve`, `betti_curvature` |
| Naive (slice-wise) partitions of spectral weight | `partition_marginals`, `slice_partition_marginals`, `merge_energy_bands!` |
| Physical band segmentation | `adaptive_autosplitbands2d`, `components_to_band_marginals`, `splitbands` |
| Resolution-matched rebinning | `rebin_to_resolution`, `upsample_mask` |
| Ridge topology | `build_ridge_topology`, `ridge_dispersions` |
| Child-ridge recovery | `find_ridge_candidates`, `recover_child_ridges`, `child_ridge_dispersions`, `RidgeRecoveryParams` |
| Parent-band domains | `restrict_to_domain`, `domain_marginals`, `analyze_domain_topologies`, `summarize_domain_branches`, `robust_branch_count` |
| Quantum Fisher information | `qfi_kernel`, `qfi`, `qfi_bound`, `nqfi`, `nqfi_lower_bound`, `qfi_witness`, `band_qfi`, `band_supports`, `band_dependence`, `symmetrize_q`, `to_smag`, `magnetic_form_factor` |
| Linear spin-wave theory with `Sunny.jl` | `LSWT`, `convert4D`, `collapse`, `La2CuO4`, `TMMC` |
| Robust statistics | `MAD` |
| HDF5 and NeXus inspection | `saveHDF5`, `loadHDF5`, `h5tree` |
| Mantid MDHistograms and MDNorm | `load_mdhisto`, `load_mdnorm`, `eachslice_mdnorm`, `intensity`, `uncertainty`, `measured` |

The cubical complex and the filtrations built on intensity data are described on the
[Conventions](https://alt-f4-dev.github.io/PHysicalTDA.jl/dev/conventions/) page; each function's
docstring documents its keywords.

## Mantid histograms

```julia
using PHysicalTDA

path_D, path_N = "sample_D.nxs", "sample_N.nxs"
h5tree(path_D; keypath="/MDHistoWorkspace/data", attributes=true, values=true, maxchildren=16)
h = load_mdnorm(path_D, path_N; path_I="sample_I.nxs")
I, σ, valid = intensity(h), uncertainty(h), measured(h)
energy_dim = energy_axis(h)
Ei = incident_energy(h)
```

Use `load_mdhisto` for a single workspace. The histogram carries `edges`, `centres`, `names`, `units`, and `frames` in `D0, D1, …` axis order. Older files with named axis datasets use the order presented by HDF5.jl. `hkl_centres(h)` returns a matrix with one stored HKL column per momentum bin. Integrated axes retain a single bin. Edge arrays are used directly, while centre arrays are converted to edges using adjacent midpoints and outer half-widths. A singleton centre requires a positive `bin_width` attribute. When an axis order is ambiguous, provide a `signal` dataset `axes` attribute, as Mantid does.

For MDNorm, keep `h.D`, `h.D2`, and `h.N` when combining bins. `intensity(h)` gives `D/N` and `uncertainty(h)` gives `sqrt(D2)/N` where `N > 0`. Both return `NaN` outside measured bins. If the numerator file has no variance, `h.D2` and `uncertainty(h)` are `nothing`. The optional `path_I` checks the saved ratio on measured bins.

Use `lattice`, `orientation_matrix`, `incident_energy`, `projection`, `goniometer`, and `run_title` for experiment metadata. Each accessor takes a zero-based `run` keyword identifying the Mantid `experimentN` group and returns `nothing` if its metadata is absent. Histograms retain the records under `h.meta.runs`; use `metadata=false` to skip them.

For a large `(H, K, E, L)` workspace, pass contiguous bin ranges or `:` in axis order. A slice callback receives a histogram with a singleton sliced dimension and its physical centre.

```julia
subset = load_mdnorm(path_D, path_N; ranges=(:, :, 5:20, 3:3))
eachslice_mdnorm(path_D, path_N; dim=4, ranges=(:, :, 5:20, :), metadata=false) do slice, L
    @info "Momentum slice" L measured_bins=count(measured(slice))
end
```

## MDNorm reductions

Fit Bragg positions in an elastic energy window and use the inverse affine map to place bins in cuts and symmetry folds. Magnetic conversion uses the measured momentum and the Cartesian moment axis. The reciprocal basis includes `2π`.

```julia
peaks = find_bragg_peaks(h; positions=bragg_hkl, energy_window=(-0.5, 0.5))
fit = fit_index_correction(peaks)
hkl = corrected_indices(h, fit.M, fit.t)
g = magnetic_conversion(h; ion=:Mn2, channel=:transverse)
cut = cut_path(h, [(0.0, 0.0), (0.5, 0.5)]; width=0.04, step=0.02, indices=hkl, conv=g)
```

`rebin`, `integrate`, and the histogram overload of `rebin_to_resolution` combine numerator, variance, and exposure before forming intensity. Use explicit factors for multidimensional grids. `weighted_average` and `magnetic_conversion(h; average_over=l_dim)` retain singleton dimensions for subsequent mapping with `map_nearest_bins`.

`cut_path` and `fold_to_zone` accept reciprocal-space point-operation matrices. Symmetry images sharing a destination count once per source bin. The folded histogram records shared-source covariance, and later reductions account for it. A 2D cut sums over any omitted momentum coordinate. A diagonal selection `abs(H-K) ≤ w` corresponds to `width=w/sqrt(2)`.

With conversion and background, the cut evaluates `sum(D - N*b)/sum(N*g)`. A background can be a scalar, a histogram-shaped array, or a callable of measured `(|Q|, ω)`. Background uncertainties returned by a callable are treated as independent between source bins.

## Quantum Fisher information

The QFI functions implement Sec. S2 of the Supplemental Material:
`F_Q(q, T) = ∫₀^∞ K_T(ω) S_mag(q, ω) dω` with `K_T(ω) = 4 tanh(βω/2)(1 − e^{−βω})`,
`nQFI = F_Q / (12 S² N)`, band-resolved contributions over disjoint parent-band supports, and the
witness criterion `G_m`. Energies are in meV and temperatures in K.

```julia
T  = 10.0                                          # K
r  = qfi(Iqω, ωs, T; ωrange=(0.5, Inf))            # F_Q(q, T), elastic line excluded
br = band_qfi(Iqω, seg, ωs, T; ωrange=(0.5, Inf))  # F_Q^(b), residual F_Q^(res), ξ_b
br.ξ                                               # nparents × nq relative contributions

# Certification needs S_mag in absolute units per magnetic ion, with uncertainties.
# The form factor is always given explicitly: a function of |Q| or ⟨j₀⟩ coefficients.
Smag, σS = to_smag(Iqω, Qmag; σ=σI, form_factor=j0_coeffs, scale=absolute_scale)
Ssym, σsym, _ = symmetrize_q(Smag, qpath; σ=σS)   # needed without inversion symmetry
w = qfi_witness(qfi(Ssym, ωs, T; σ=σsym, ωrange=(0.5, Inf)); S=5/2, m=1)
# Several spin sizes per normalization volume: qfi_witness(r; spins=[S₁, S₂, ...], m=1)
```

Nothing in the QFI functions is specific to a material:

- the magnetic form factor is a required argument of `to_smag` (`MN2_J0` is provided for Mn²⁺);
- the separable bound is set by a single spin `S` with `N` sites, by per-site `spins` for
  several magnetic species, or by an explicit `bound` (`qfi_bound`);
- `symmetrize_q` averages `S(q, ω)` and `S(−q, ω)`, which gives the QFI without assuming inversion
  symmetry (Eq. S25 assumes it).

The relative contributions `ξ_b` do not depend on the absolute scale. `nQFI` and `G_m` require an
isolated, absolutely normalized magnetic response; `to_smag` performs the conversion only when called.

## Rebinning and child-ridge recovery

Measured maps are usually rebinned to about one resolution width before segmentation. The QFI is
then computed on the original grid with the parent-band supports mapped back to it.

```julia
rb  = rebin_to_resolution(Iqω, qs, ωs; σ=σI, remainder=:merge)  # NaN bins not averaged; no bins dropped
seg = adaptive_autosplitbands2d(rb.I, rb.q, rb.ω; resolution_q=rb.resolution_q, resolution_ω=rb.resolution_ω)
disp = ridge_dispersions(seg, rb.q, rb.ω; centre=:weighted)  # parent-band dispersions, no re-segmentation
br  = band_qfi(Smag, band_supports(seg, rb), ωs, T; σ=σS)    # supports on the original grid

# Individual ridges inside each parent, with split, merge and crossing points
cr = recover_child_ridges(rb.I, rb.q, rb.ω)                  # RidgeRecoveryParams() = manuscript settings
child_ridge_dispersions(cr).segments
```

`rebin_to_resolution(...; average=:all)` and the default `RidgeRecoveryParams()` reproduce the
manuscript workflow exactly. The default ridge settings connect a ridge across momentum rows only
when it moves by about one energy bin or less between rows; for steep dispersions on coarse grids,
`RidgeRecoveryParams(ωs; max_slope)` scales the support widths and tolerances to the energy bin.

## Optional conventions

Each of these keywords defaults to the established behaviour, so existing results do not change unless an
alternative is requested.

| Keyword | Functions | Default | Alternative |
|---|---|---|---|
| `essential_cap` | `pd_array_intensities`, `pd_sunny_intensities` (with `make_death_finite=true`) | `:max_death`: each infinite death is capped at the largest finite death in that degree | `:threshold`: cap at the filtration value where the complex is complete (for a superlevel filtration, the background intensity) |
| `interval` | `betti_curve`, `betti_curvature` | `:closed`: a class is alive on [b, d] | `:halfopen`: alive on [b, d), the rank H_p convention |
| `peak_selection` | `partition_marginals`, `slice_partition_marginals` | `:height`: partition centres at the highest local maxima of each slice | `:persistence`: partition centres at the births of the most persistent H₀ classes (same number of partitions) |

## Unreleased

- New: `rebin_to_resolution`, `upsample_mask` and `band_supports(result, rebin)` (resolution-matched
  rebinning and supports on the original grid).
- New: `ridge_dispersions(result, qs, ωs; centre)` uses an existing segmentation; `centre=:weighted`
  gives intensity-weighted parent-band energies. The existing method is unchanged.
- New: child-ridge recovery (`find_ridge_candidates`, `recover_child_ridges`, `child_ridge_dispersions`,
  `RidgeRecoveryParams`), ported from the manuscript scripts.
- `rebin_to_resolution(...; remainder=:merge)` adds trailing bins that do not fill a block to the
  last block, so every original bin keeps a parent-band support. The default, `:drop`, is unchanged.
- `RidgeRecoveryParams(; include_energy_edges=true)` lets `find_ridge_candidates` and
  `recover_child_ridges` accept maxima in the first or last energy bin, for dispersion extrema at
  the edge of a clipped energy window. The default, `false`, is unchanged.
- No existing function changed behaviour.

## Changes in 0.1.2

**Behaviour changes**

- `periodic_q` now defaults to `false` throughout physical band segmentation (`adaptive_autosplitbands2d`
  and the routines it calls), so momentum axes are treated as open unless periodicity is requested.
  Pass `periodic_q=true` to recover the previous behaviour.
- `MAD(t)` returns the scalar median absolute deviation. `MAD(t; elementwise=true)` returns the
  scaled elementwise deviations that earlier versions returned.
- Plotting moved to a package extension (see [Plotting](#plotting)). Makie is no longer a hard
  dependency, and `PHysicalTDA.TopoViz` no longer exists.

**Renamed, with the old names kept as aliases**

- `naivesplitbands` → `slice_partition_marginals`
- `autosplitbands` → `partition_marginals`

The old names refer to the same functions, so existing code and results are unchanged.

**Other**

- Segmentation accepts range axes for `qs` and `ωs`.
- Threaded segmentation uses task-local scratch buffers, which are safe under Julia's task migration.
- `threshold` accepts any `Real`.
- Repaired plotting helpers, examples, and documentation.
- New tests for segmentation and the optional conventions.
- Supports Julia 1.10 (LTS) onward; CI tests 1.10, 1.11, 1.12 and the latest release.

## Plotting

Plotting helpers (`plot_persistence_diagram`, `plot_betti_surface`, ...) are provided by a package extension that loads when Makie does. Add a backend to your environment and load it before plotting:

```julia
using PHysicalTDA, CairoMakie      # or GLMakie
fig = plot_persistence_entropy(dataset; backend=:cairo)
```
