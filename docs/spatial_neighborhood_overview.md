# Spatial Module 06 Neighborhood Overview

`06d_spatial_neighborhood.R` is H5AD-first. It first resolves a spatial H5AD mirror from `results/90a_export_h5ad/${SPATIAL_NEIGHBORHOOD_H5AD_MODULE:-spatial_03_region}` and runs Squidpy with `${SPATIAL_NEIGHBORHOOD_GROUP_BY:-region_label}`. The legacy annotated/subannotated ST panorama RDS conversion path remains as a fallback only when no H5AD mirror is available. The sidecar `workflow/04python/spatial_neighborhood.py` requires AnnData `obsm["spatial"]` plus the selected `obs` label column, then runs Squidpy graph functions directly.

| Output | Meaning |
| --- | --- |
| `neighborhood_interaction_matrix.tsv` | `sq.gr.interaction_matrix()` label-label neighbor counts. |
| `neighborhood_nhood_enrichment_zscore.tsv` | `sq.gr.nhood_enrichment()` permutation z-scores using `SPATIAL_NEIGHBORHOOD_PERMS`. |
| `neighborhood_co_occurrence.tsv` | `sq.gr.co_occurrence()` distance-binned co-occurrence rows. |
| `neighborhood_manifest.tsv` | Single-row status contract for downstream niche derivation. |

The module prefers `sub_region` labels when present, then falls back to `region` and `seurat_clusters`. Spatial coordinates are expected to arrive through the AnnData spatial bridge rather than through a hand-written coordinates TSV.

If the annotated panorama, Seurat, an AnnData export bridge, Python environment, Squidpy, or enough labeled spots are unavailable, the manifest records a skip status. The pipeline sets the `spatial_neighborhood` gate after this stage.
