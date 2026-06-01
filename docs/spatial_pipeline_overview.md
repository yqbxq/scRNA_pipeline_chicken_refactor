# Spatial Pipeline Overview

The main ST workflow currently runs:

| stage | manifest | purpose |
| --- | --- | --- |
| `01` | `spatial_01_build_objects` | Build per-section spatial Seurat objects. |
| `01a` / `01b` / `01c` | pre/post QC manifests | Review, filter, and re-check spatial spots. |
| `02` / `02a` / `02b` | normalization, integration, clustering manifests | Compare normalization, integration, and clustering choices. |
| `03` / `03a` | `spatial_03_region_annotation`, `spatial_03a_region_annotation_eda` | Assign whole-panorama regions and review diagnostics. |
| `04a` / `04b` / `04c` | `spatial_04a_subcluster_build`, `spatial_04b_subcluster_annotate`, `spatial_04c_subcluster_eda` | Optional region-subset subclustering and sub-region review. |
| `05` | `spatial_05_region_marker` | Region and optional sub-region marker discovery. |
| `05a` | `spatial_05a_region_pseudobulk` | Formal region pseudobulk DE when the replicate gate passes, with spot-level fallback links otherwise. |
| `05b` | `spatial_05b_region_composition` | Descriptive region composition and formal propeller output when the replicate gate passes. |
| `05_spatial_de` | `spatial_05_spatial_de` | Exploratory spot-level spatial DE for N=2 rescue. |
| `06a` / `06b` / `06c` | `spatial_06a_region_go`, `spatial_06b_region_kegg`, `spatial_06c_region_enrichment_eda` | Region GO/KEGG enrichment and review summaries. |
| `06d` | `spatial_06d_neighborhood` | Squidpy-gated radius neighborhood enrichment and co-occurrence tables. |
| `06e` | `spatial_06e_niche` | Post-07 spatial niche derivation from neighborhood and deconvolution features. |
| `07a` / `07b` / `07c` / `07d` / `07e` / `07f` | RCTD, transfer, CARD, cell2location, compare, validation manifests | Default multi-method deconvolution, method recommendation, and validation contracts. |

The `spatial_region_annotation` gate is held twice: once after `03a` for region review, then reset and held again after `04c` for sub-region review. ST 05 adds `spatial_marker_de` after marker, DE, and composition outputs are written. ST 06 adds `spatial_enrichment` after `06c`, `spatial_neighborhood` after `06d`, and `spatial_niche` after post-07 `06e`. ST 07 holds `spatial_deconv` after method comparison and optional validation. Module 08 regulatory extensions remain pending as an independent follow-up module.

Additional ST06 docs:

- `docs/spatial_enrichment_overview.md`
- `docs/spatial_neighborhood_overview.md`
- `docs/spatial_niche_overview.md`
- `docs/spatial_deconv_overview.md`
- `docs/deconv_pairs_schema.md`
- `docs/spatial_deconv_compare.md`
- `docs/spatial_deconv_validation.md`
