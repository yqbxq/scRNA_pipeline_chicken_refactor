# Spatial Module 07 Deconvolution Overview

Spatial module 07 runs deconvolution from a frozen single-cell reference into ST spots.

| Script | Manifest | Purpose |
| --- | --- | --- |
| `07a_deconvolution_rctd.R` | `spatial_07a_deconvolution_rctd` | Primary RCTD/spacexr slot. |
| `07b_deconvolution_transfer.R` | `spatial_07b_deconvolution_transfer` | Seurat TransferData fallback slot. |
| `07c_deconvolution_card.R` | `spatial_07c_deconvolution_card` | CARD alternative method slot. |
| `07d_deconvolution_cell2location.R` | `spatial_07d_deconvolution_cell2location` | cell2location Python sidecar slot. |
| `07e_deconvolution_compare.R` | `spatial_07e_deconvolution_compare` | Multi-method status/ranking and `recommended_method.txt`. |
| `07f_deconvolution_validation.R` | `spatial_07f_deconvolution_validation` | Optional scDesign3 synthetic validation. |

`60_run_spatial_extensions.sh` runs 07a and 07e by default. `--with-deconv-extra` adds 07b/07c/07d, and `--with-deconv-validation` adds 07f. The `spatial_deconv` gate is held after 07e or optional 07f.

All method scripts write manifests even when `deconv_pairs.tsv` is empty or method packages are unavailable. This keeps downstream review and 06e niche wiring deterministic while server environments are prepared.
