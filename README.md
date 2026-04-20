<p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/dark-mode-banner.svg">
      <img src="docs/src/assets/light-mode-banner.svg" alt="PHysicalTDA.jl" width="600"/>
    </picture>
</p>


[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alt-f4-dev.github.io/PHysicalTDA.jl/dev/)
[![Build Status](https://github.com/alt-f4-dev/PHysicalTDA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/alt-f4-dev/PHysicalTDA.jl/actions/workflows/CI.yml?query=branch%3Amain)

`PHysicalTDA.jl` performs Topological Data Analysis (TDA) by utilizing Persistence Homology (PH) on scattering datasets.
It is written to be handle both synthetic and raw datasets. The module currently supports Julia native and `Sunny.jl` types.
It will support NXS, JLD2, and ASCII datasets for analysis of raw measurements in the future. 
