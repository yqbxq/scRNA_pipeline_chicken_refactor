# Spatial Module 05 Marker / DE Overview

Spatial module 05 consumes the annotated panorama from module 03 and, when present, the subannotated panorama from module 04. Each script writes region outputs, and repeats the same contract for `sub_region` when that column exists.

| Script | Manifest | Purpose |
| --- | --- | --- |
| `05_region_marker_discovery.R` | `spatial_05_region_marker` | Exploratory region or sub-region marker discovery with `FindAllMarkers`. |
| `05a_region_pseudobulk.R` | `spatial_05a_region_pseudobulk` | Formal region pseudobulk DE when the biological replicate gate passes. |
| `05b_region_composition.R` | `spatial_05b_region_composition` | Region composition proportions, with formal propeller output only when the replicate gate passes. |
| `05_spatial_de.R` | `spatial_05_spatial_de` | Spot-level cross-condition Wilcoxon DE used as the exploratory fallback for N=2 studies. |

The default ST project has one SYF and one F5 sample, so formal replicate-aware paths usually write `skipped_replicate_gate`. The spot-level rescue remains available and is linked from the pseudobulk manifest through `exploratory_fallback`.

Outputs are written under `results/spatial/tables/spatial_05_*` and `results/spatial/figures/spatial_05_*`. The review report lands in `reports/eda/spatial_marker_de/`, and the pipeline holds at gate `spatial_marker_de` before enrichment or deconvolution.
