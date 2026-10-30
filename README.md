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
It will support NXS, JLD2, and ASCII datasets for analysis of raw measurements in the future.

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
| Ridge topology | `build_ridge_topology`, `ridge_dispersions` |
| Parent-band domains | `restrict_to_domain`, `domain_marginals`, `analyze_domain_topologies`, `summarize_domain_branches`, `robust_branch_count` |
| Linear spin-wave theory with `Sunny.jl` | `LSWT`, `convert4D`, `collapse`, `La2CuO4`, `TMMC` |
| Robust statistics | `MAD` |
| File I/O | `saveHDF5`, `loadHDF5`, `h5tree` |

The cubical complex and the filtrations built on intensity data are described on the
[Conventions](https://alt-f4-dev.github.io/PHysicalTDA.jl/dev/conventions/) page; each function's
docstring documents its keywords.

## Optional conventions

Each of these keywords defaults to the established behaviour, so existing results do not change unless an
alternative is requested.

| Keyword | Functions | Default | Alternative |
|---|---|---|---|
| `essential_cap` | `pd_array_intensities`, `pd_sunny_intensities` (with `make_death_finite=true`) | `:max_death`: each infinite death is capped at the largest finite death in that degree | `:threshold`: cap at the filtration value where the complex is complete (for a superlevel filtration, the background intensity) |
| `interval` | `betti_curve`, `betti_curvature` | `:closed`: a class is alive on [b, d] | `:halfopen`: alive on [b, d), the rank H_p convention |
| `peak_selection` | `partition_marginals`, `slice_partition_marginals` | `:height`: partition centres at the highest local maxima of each slice | `:persistence`: partition centres at the births of the most persistent H₀ classes (same number of partitions) |

## Changes in 0.1.1

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
