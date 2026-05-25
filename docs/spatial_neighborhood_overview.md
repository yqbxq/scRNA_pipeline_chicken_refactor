# Spatial Module 06 Neighborhood Overview

`06d_spatial_neighborhood.R` consumes the annotated or subannotated ST panorama, exports a transient AnnData `.h5ad` via `sceasy::convertFormat()` or `MuDataSeurat::WriteH5AD()`, and deletes that temporary file after the Python sidecar finishes. The sidecar `workflow/04python/spatial_neighborhood.py` requires AnnData `obsm["spatial"]` plus the selected `obs` label column, then runs Squidpy graph functions directly.

| Output | Meaning |
| --- | --- |
| `neighborhood_interaction_matrix.tsv` | `sq.gr.interaction_matrix()` label-label neighbor counts. |
| `neighborhood_nhood_enrichment_zscore.tsv` | `sq.gr.nhood_enrichment()` permutation z-scores using `SPATIAL_NEIGHBORHOOD_PERMS`. |
| `neighborhood_co_occurrence.tsv` | `sq.gr.co_occurrence()` distance-binned co-occurrence rows. |
| `neighborhood_manifest.tsv` | Single-row status contract for downstream niche derivation. |

The module prefers `sub_region` labels when present, then falls back to `region` and `seurat_clusters`. Spatial coordinates are expected to arrive through the AnnData spatial bridge rather than through a hand-written coordinates TSV.

If the annotated panorama, Seurat, an AnnData export bridge, Python environment, Squidpy, or enough labeled spots are unavailable, the manifest records a skip status. The pipeline sets the `spatial_neighborhood` gate after this stage.
