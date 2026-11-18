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
| Background | `exclusion_mask`, `momentum_magnitude`, `background_profile`, `background_shells`, `TabulatedBackground`, `magnetic_histogram` |
| Absolute normalization | `neutron_xs`, `formula_xs`, `elastic_fwhm`, `incoherent_scale`, `radial_elastic_profile`, `fit_powder_lines`, `powder_line_ratio`, `powder_line_flat`, `container_corrected_scale`, `first_moment_scale`, `total_moment_scale`, `compare_scales`, `select_scale` |

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

`cut_path` and `fold_to_zone` accept reciprocal-space point-operation matrices. Symmetry images sharing a destination count once per source bin. A folded cut includes every lattice translation of an image that reaches the path, so paths may end on zone boundaries such as `(1/2, 1/2)`. When the path length is not a multiple of `step`, `remainder=:keep`, `:merge`, or `:drop` controls the short final bin. The folded histogram records shared-source covariance, and later reductions account for it. A 2D cut sums over any omitted momentum coordinate. A diagonal selection `abs(H-K) ≤ w` corresponds to `width=w/sqrt(2)`.

`cut_path(...; metric=:cartesian)` takes `width` and `step` in Å⁻¹ and returns `s` in Å⁻¹, which suits paths that include `L`. `fold` also accepts a `Sunny.Crystal`; `reciprocal_point_group(cryst; dims=2)` lists its operations on `(H, K)`. `cut_histogram(cut)` turns a cut into a two-axis `NormalizedHistogram`, so `rebin` and `rebin_to_resolution` combine its bins as `Σ(D - N b)/Σ(N g)`.

With conversion and background, the cut evaluates `sum(D - N*b)/sum(N*g)`. A background can be a scalar, a histogram-shaped array, or a callable of measured `(|Q|, ω)`. Background uncertainties returned by a callable are treated as independent between source bins in `σ`; `σ_background` gives the bound for a background error fully correlated within each cut bin.

## Background

Build an exclusion mask from Bragg positions, zone centres, a model dispersion and `|Q|` limits, then estimate the background from the remaining bins. `background_profile` gives a momentum-independent `b(ω)`; `background_shells` gives `b(|Q|, ω)` in `|Q|` shells, filling cells with too few bins along energy, then along `|Q|`, then from a fallback profile. Both return a `TabulatedBackground`, which `cut_path`, `to_smag` and `magnetic_histogram` accept as a background.

```julia
bragg = [((0, 0), 0.10), ((0.5, 0.5), 0.10)]
ω_model(q) = 4J * S * sqrt(max(1 - ((cospi(2q[1]) + cospi(2q[2])) / 2)^2, 0))
profile_mask = exclusion_mask(h; sites=bragg, dispersion=ω_model, window=1.0, qmax=2.5, indices=hkl)
shell_mask = profile_mask .| exclusion_mask(h; sites=[((0, 0), 0.20), ((0.5, 0.5), 0.20)], indices=hkl)
b0 = background_profile(h; exclude=profile_mask, fill=true)
b = background_shells(h; exclude=shell_mask, dq=0.15, qmax=2.5, min_bins=40, smooth=0.5, smooth_from=1.0, fallback=b0)
cut = cut_path(h, [(0.0, 0.0), (0.5, 0.5)]; width=0.025/sqrt(2), step=0.0125*sqrt(2), indices=hkl, conv=g, background=b)
```

`magnetic_histogram(h, g; background)` returns a `NormalizedHistogram` with numerator `D - N*b` and normalization `N*g`, so `fold_to_zone`, `rebin` and `integrate` reduce `S_mag` directly. Masks are `true` where a bin is excluded; combine them with `.|`.

## Absolute normalization

Cross sections come from Sears (1992). The incoherent elastic line of the sample sets one scale; powder lines of a container measure its share of the same line; sum rules set two more. `compare_scales` and `select_scale` tabulate the estimates and choose the most conservative.

```julia
inc = incoherent_scale(h; exclude=profile_mask, formula="Rb2MnF4", magnetic_ions=1)
radial = radial_elastic_profile(path_D, path_N; dim=4, energy_window=inc.window, dq=0.02, qmax=5.0)
τs = [2π * norm(hkl) / 4.032 for hkl in ((1, 1, 1), (2, 0, 0))]
fit = fit_powder_lines(radial.Q, radial.R, τs)
al = powder_line_flat(fit, [powder_line_ratio(:Al, hkl; a=4.032) for hkl in ((1, 1, 1), (2, 0, 0))])
κ_al = container_corrected_scale(inc.κ, al.flat / inc.integral, al.σflat / inc.integral)
m1 = first_moment_scale(cut.values, cut.σ, cut.ω, S_model; T=5.0, window=(0.45, 7.5), correction=r₁, correction_σ=0.01)
zone = fold_to_zone(magnetic_histogram(h, g; background=b), group; grid=(0:0.025:1, 0:0.025:1), indices=hkl)
m0 = total_moment_scale(intensity(zone), S_model_zone, cut.ω; window=(0.45, 7.5), correction=0.985, correction_σ=0.015)
table = compare_scales([(name="incoherent", κ=inc.κ), (; name="Al-corrected", κ_al...),
                        (name="first moment", κ=m1.κ, κ_lo=m1.κ_lo), (name="total moment", κ=m0.κ, κ_lo=m0.κ_lo)])
κ = select_scale(table; among=("first moment", "total moment")).used
```

`κ` converts `(D/N)/g` into `S_mag` per magnetic ion. The incoherent scale is `σ_inc/(4π I_el (γr₀)² n_mag)`, with `I_el` the elastic line above the lower of its two flank levels. A powder line's `|Q|`-integrated area, divided by `powder_line_ratio`, is that material's flat elastic level. The sum rules compare measured and model moments; `correction` carries the ratio of the true to the model moment. The low end `κ_lo` includes `correction_σ` and, for the first moment, the spread over momentum. Absorption, multiple scattering and Debye–Waller corrections are not applied.

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

Systematic uncertainties enter the witness alongside the statistical one:

```julia
r = qfi(κ .* cut.values, cut.ω, T; σ=κ .* cut.σ, σ_sys=κ .* cut.σ_background, ωrange=(0.45, Inf))
qfi_uncertainty(r; scale_σ=σκ/κ)                   # statistical, correlated background and scale, in quadrature
w = qfi_witness(r; S=5/2, m=1, scale_σ=σκ/κ)       # nQFI_lb = (F_Q − nσ σ_tot)/B
sens = qfi_sensitivity(Smag, ωs, T; σ=σS, lower_limits=[0.35, 0.45, 0.55, 0.65], S=5/2)  # rows: E_min
compare_band_qfi(band_qfi(Smag, coarse_masks, ωs, T; σ=σS), band_qfi(Smag, refined_masks, ωs, T; σ=σS))
```

`σ_sys` is fully correlated across the energies of a momentum row, so it adds linearly over the
energy integral (`r.σsys`); `scale_σ` is the relative uncertainty of the absolute scale. With neither,
`nqfi_lower_bound` and `qfi_witness` are unchanged. `qfi_sensitivity` includes energies above each
lower limit, `ω > E_min`. `compare_band_qfi` returns `Δξ`, a conservative `σΔξ` and whether the two
sets of supports agree within `nσ`.

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

`rebin_to_resolution(...; average=:all)` and `RidgeRecoveryParams(; include_energy_edges=false)` reproduce
the manuscript workflow exactly; its ridge scripts must also pass `include_energy_edges=false` to
`build_ridge_topology`, whose default is now `true`. The default ridge settings connect a ridge across momentum rows only
when it moves by about one energy bin or less between rows; for steep dispersions on coarse grids,
`RidgeRecoveryParams(ωs; max_slope)` scales the support widths and tolerances to the energy bin.

Refine a parent support on the original grid and recover its child ridges using the measured-bin mask.

```julia
parents = refine_supports(Iqω, σI, valid, qs, ωs, seg; rebin=rb,
                          resolution=(rb.resolution_q, rb.resolution_ω), bridge=0.15)
params = RidgeRecoveryParams(ωs, qs; max_slope=12)
children = recover_child_ridges(Iqω, qs, ωs; params, parents, valid)
trajectories = mask_dispersions(Iqω, parents, qs, ωs; valid)
```

The resolution tuple gives Gaussian FWHM along (q) and energy. Set `bridge` to the maximum expected energy displacement between neighbouring momentum rows in the support. Refinement returns one mask per input parent, and recovery retains that order. The trajectory result flags estimates without positive intensity weight and rows filled by interpolation.

The segmentation can also run on the unrebinned grid. Give the instrumental resolution (a number or a function of energy), the steepest expected slope, and `physical_units=true`, which counts smoothing radii, the prominence noise estimate and the ridge tube in resolution elements rather than pixels.

```julia
seg = adaptive_autosplitbands2d(Iqω, qs, ωs; resolution_q=0.025, resolution_ω=ω -> 0.22 + 0.01ω,
                                max_slope=45.0, physical_units=true)
resolution_report(seg)   # Γq, Γω, the apparent Γω, their sources, and where a function was evaluated
seg.qfirst.track_lengths # rows occupied by each peak track; short tracks do not add parent bands
```

`ridge_max_gap_rows` lets the ridge graph bridge rows where a ridge has no peak (for example a masked column).
`max_slope` (energy per unit of `qs`) links peaks of neighbouring rows one to one, closest first, within `max_slope·Δq`. Parent bands are counted from peaks on tracks that occupy at least `min_track_fraction` (default 0.25) of the rows, so short secondary features do not add parents. It also joins support components and ridge-graph nodes across that displacement. Without a supplied width, `resolution_estimate=:flat` estimates it from the columns where the strongest feature moves least, which on a steep dispersion is closer to the instrumental width than the default median half-height width. `ridge_long_beta0` counts the ridge-graph components that occupy at least `min_track_fraction` of the rows with peaks.

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
- New: background estimation (`exclusion_mask`, `momentum_magnitude`, `background_profile`,
  `background_shells`, `TabulatedBackground`, `magnetic_histogram`) and absolute normalization
  (`neutron_xs`, `formula_xs`, `elastic_fwhm`, `incoherent_scale`, `radial_elastic_profile`,
  `fit_powder_lines`, `powder_line_ratio`, `powder_line_flat`, `container_corrected_scale`,
  `first_moment_scale`, `total_moment_scale`, `compare_scales`, `select_scale`).
- New: `adaptive_autosplitbands2d(...; max_slope, min_track_fraction, physical_units, resolution_estimate)`,
  a callable `resolution_ω`, `resolution_report`, and the result fields `apparent_resolution_ω`,
  `resolution_source`, `resolution_energy` and `ridge_long_beta0`. All are off by default.
- Changed default: `include_energy_edges = true` in `ridge_peak_slices`, `build_ridge_topology`,
  `analyze_domain_topologies` and `RidgeRecoveryParams()`, matching `adaptive_autosplitbands2d` and
  `qfirst_parent_domains`. `adaptive_autosplitbands2d` now passes its setting to its ridge topology, and
  the ridge graph drops isolated maxima in the first or last energy bin. To reproduce earlier results,
  pass `include_energy_edges=false` to those functions and `ridge_include_energy_edges=false` to
  `adaptive_autosplitbands2d`.

- New: `cut_path(...; metric=:cartesian)`, `cut_histogram`, `reciprocal_point_group` and a
  `Sunny.Crystal` as `fold`; cut results also carry `s_edges`, `ω_edges`, `numerator_covariance` and
  `s_units`. `build_ridge_topology(...; max_gap_rows)` and `adaptive_autosplitbands2d(...;
  ridge_max_gap_rows)` bridge rows without a peak (default 0, unchanged).
- New: `plot_qfi_summary` and `plot_ridge_workflow` (Makie extension), and
  `magnetic_conversion(...; l_average=(values, weights))`, which averages the polarization factor over `L` at
  each bin's own `(H, K)`.
- New: `qfi(...; σ_sys)`, `QFIResult.σsys`, `BandQFIResult.σsys`, `qfi_uncertainty`, and `scale_σ`
  and `systematic` in `nqfi_lower_bound`, `qfi_witness` and `band_dependence`; `qfi_sensitivity` and
  `compare_band_qfi`. Results without these options are unchanged.

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

Two layouts assemble the analysis results into figures. Material details are keywords: the separable bound
(`normalization`), model curves, marked momenta and shaded intervals, labels.

```julia
fig, ax = plot_qfi_summary(r, hs; normalization=(; S=5/2), scale_σ=σκ/κ, bands=br, sensitivity=sens,
                           comparison=compare_band_qfi(br, br_refined), models=(; LSWT=FQ_lswt),
                           markers=(0.5,), xlabel="(h, h) (r.l.u.)", backend=:cairo)
fig, ax = plot_ridge_workflow(I, hs, ωs; valid, segmentation=seg, candidates, children, supports=refined,
                              trajectories, markers=(0.5,), backend=:cairo)
```

Each returns the figure and its axes by panel name; panels appear only for the results given.
