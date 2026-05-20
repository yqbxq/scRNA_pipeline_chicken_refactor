# Spatial Module 06 Neighborhood Overview

`06d_spatial_neighborhood.R` consumes the annotated or subannotated ST panorama and exports a compact spot table with coordinates and region labels. The Python sidecar `workflow/04python/spatial_neighborhood.py` checks the `py_spatial` environment for Squidpy, then computes a radius-based interaction matrix, permutation z-score matrix, and mean-distance co-occurrence table.

| Output | Meaning |
| --- | --- |
| `neighborhood_interaction_matrix.tsv` | Observed label-label neighbor counts within `SPATIAL_NEIGHBORHOOD_RADIUS`. |
| `neighborhood_nhood_enrichment_zscore.tsv` | Permutation z-scores using `SPATIAL_NEIGHBORHOOD_PERMS`. |
| `neighborhood_co_occurrence.tsv` | Mean neighbor distance for each label pair. |
| `neighborhood_manifest.tsv` | Single-row status contract for downstream niche derivation. |

The module prefers `sub_region` labels when present, then falls back to `region` and `seurat_clusters`. Coordinate columns are detected from common ST names such as `x/y`, `imagecol/imagerow`, or Visium full-resolution pixel columns.

If the annotated panorama, Seurat, Python environment, Squidpy, or enough labeled spots are unavailable, the manifest records a skip status. The pipeline sets the `spatial_neighborhood` gate after this stage.
