# Spatial Module 06 Niche Overview

Spatial niche derivation is implemented as a post-07 gated ST06 step. `06e_niche_derivation.R` runs after deconvolution comparison, consumes the 06d neighborhood manifest plus the recommended 07 deconvolution output, and records an explicit skip status until a deconvolution method manifest has `status=ok`.

The planned 06e contract will consume:

| Input | Purpose |
| --- | --- |
| `spatial_06d_neighborhood` | Region-neighborhood features. |
| `spatial_07a_deconvolution_rctd` or fallback deconvolution manifest | Per-spot cell-type proportion features. |
| `spatial_07e_deconvolution_compare` | Recommended deconvolution method selection. |

When real deconvolution proportions are present, 06e derives `deconv_proportion_kmeans` niche labels, writes `niche_scores.tsv` and `niche_labels.tsv`, and saves `spatial_panorama_niched.rds`. It adds `spatial_niche*` metadata without overwriting `region`, `sub_region`, or `condition`. The `spatial_niche` gate is held after 06e before joint analyses.
