```@meta
CurrentModule = PHysicalTDA
```

# Spectroscopic Data as a Discretized Manifold
PHyicalTDA operates on datasets of the general form $I(\mathbf{q},\omega)$ where $\hbar\mathbf{q}$ denotes momentum transfer and $\hbar\omega$ denotes energy transfer.

In neutron scattering experiments, this intensity represents a measured spectral density which encodes the allowed collective excitations of the system. From the perspective of PHysicalTDA, the dataset is treated as:
- A *scalar field* defined over a discretized parameter space
- Independent of the microscopic oigin of the signal (phonons, magnons, continua, etc.)
- A sampled representation of an underlying *excitation manifold* embedded in $(\mathbf{q},\omega)$ space

# The Cubical Complex
Experimental datasets are naturally sampled on rectilinear grids due to detector geometry, binning, or interpolation. Accordingly, PHysicalTDA represents intensity data using a 
**cubical complex**:
- Each voxel corresponds to a sampled vaue of $I(\mathbf{q},\omega)$
- The grid defines adjacency relations used to construct the topological complex
- The scalar field values are used to define a *filtration* over the complex

The cubical complex is preferred over simplicial meshes (Vietoris-Rips, alpha, etc.) because it preserves the native structure of experimental data and avoids unnecessary triangulation artifacts.

# Filtrations on Intensity Data
Persistent homology is computed by applying a filtration to the cubical complex covering the intensity data. In PHysicalTDA:
- Filtrations are either *superlevel* or *sublevel* 
- Persistence tracks when a topological feature appears and vanishes
- Filtration parameters are chosen to be invariant under monotonic rescaling of intensity

This construction allows robust identification of topological features corresponding to excitation branches, gaps, continua, and intersections.

# Excitation Manifolds
From a geometric perspective, dispersive excitations trace out lower-dimensional manifolds embedded in the higher-dimensional $(\mathbf{q},\omega)$ space. Persistent homology provides a way to characterize these structures *without fitting dispersion relations* or assuming symmetries. Topological features can describe connectivity of excitation branches, band merging/splitting, and the stability of features under noise/resolution effects. 
