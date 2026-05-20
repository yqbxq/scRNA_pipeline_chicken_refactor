# Spatial Module 06 Niche Overview

Spatial niche derivation is intentionally deferred until module 07 deconvolution manifests exist. The ST06 branch reserves configuration and manifest paths for `spatial_06e_niche`, but does not invoke `06e_niche_derivation.R` in `60_run_spatial_extensions.sh`.

The planned 06e contract will consume:

| Input | Purpose |
| --- | --- |
| `spatial_06d_neighborhood` | Region-neighborhood features. |
| `spatial_07a_deconvolution_rctd` or fallback deconvolution manifest | Per-spot cell-type proportion features. |

The planned output is `spatial_panorama_niched.rds` plus niche score and label tables. It must add niche metadata without overwriting `region`, `sub_region`, or `condition`. The `spatial_niche` gate is already present in the template for the post-07 implementation.
