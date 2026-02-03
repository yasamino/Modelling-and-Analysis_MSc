# Cell Shape Measurements and Alignment

This repository contains analysis scripts for quantifying **cell shape anisotropy (elongation)** and **collective alignment** in multicellular systems.  
The implemented measurements are based on established methods from the literature and are intended for use with simulation data or segmented cell geometries.

---

## Overview

The scripts in this repository allow you to:

- Measure **single-cell shape elongation** from cell boundary geometry.
- Quantify **collective alignment** of cells across a tissue using a triangle-based order parameter.
- Visualize spatial and temporal evolution of cell shape and alignment.

These tools are suitable for epithelial tissues, cell aggregates, and other densely packed cellular systems.

---

## Example Outputs

Below are example animations generated using the analysis scripts in this repository.

### Cell Elongations
In this case the cells are colored based on the magnitude of elongation
![Cell shape alignment example](elong_colorma_XY_100_cell_mixinp.gif)

### Collective Elongation Directions Differenciated by Cell Type
![Cell alignment example](Type_diff_XY_100_cell_mixinp.gif)

---

## Implemented Measurements

### Cell Shape Elongation

Cell elongation is quantified using a tensor-based measure derived from the second moments of the cell boundary.

For each cell:
- The centroid is computed from boundary vertices.
- A shape tensor is constructed from vertex displacements.
- Eigenvalues of this tensor define the magnitude and orientation of cell elongation.

This approach follows the formulation described in **Appendix 1** of:

> Jordi Comelles et al. (2021).  
> *Epithelial colonies in vitro elongate through collective effects*.  
> **eLife**, 10:e57730.  
> https://elifesciences.org/articles/57730

This method provides a size-independent, geometry-based measure of single-cell anisotropy.

---

### Collective Cell Alignment

Collective alignment is quantified using a **triangle-based order parameter** that captures coherent orientation across the tissue.

The procedure:
- Cells are associated with a triangulation of cell centers.
- Each triangle is mapped to a reference equilateral triangle.
- A symmetric, traceless tensor is computed for each triangle.
- Area-weighted averaging yields a global alignment tensor.

This method separates:
- Local shape anisotropy
- Global orientational alignment

The implementation is based on:

> Merkel, M. & Manning, M. L. (2020).  
> *A general framework for epithelial tissue alignment*.  
> **Proceedings of the National Academy of Sciences (PNAS)**,  
> 117 (39) 24051–24061.  
> https://www.pnas.org/doi/full/10.1073/pnas.1916418117

When explicit vertex connectivity is not available, a **Delaunay triangulation of cell centers** is used as a proxy for the triangular tiling.

---

## Repository Contents

This repository includes:

- Scripts for **cell shape elongation analysis**
- Scripts for **collective alignment analysis**
- Visualization and plots
- Example animations


---

## References


1. Merkel, M. et al. (2020).  
   *Connecting cell geometries to tissue mechanics*.  
   **eLife**, 9:e57730.  
   https://elifesciences.org/articles/57730

2. Merkel, M. & Manning, M. L. (2020).  
   *A general framework for epithelial tissue alignment*.  
   **PNAS**, 117(39), 24051–24061.  
   https://www.pnas.org/doi/full/10.1073/pnas.1916418117

---

## Notes

- Alignment measurements depend on the chosen triangulation method.
- Scripts may require adaptation for specific datasets.
