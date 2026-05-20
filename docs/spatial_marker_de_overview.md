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

## Implementation Notes

- **Pseudobulk DE backend (05a)**: uses `edgeR` directly (`DGEList -> filterByExpr -> calcNormFactors -> glmQLFit -> glmQLFTest`) instead of `muscat::pbDS`. `muscat::pbDS` internally wraps the same edgeR pipeline for its default `method="edgeR"`, so this keeps the statistical path equivalent while avoiding a hard `muscat` dependency. The design matrix is `~ group` by default and `~ batch + group` when `comparisons.tsv` `batch_var` resolves to at least two non-confounded levels in retained pseudobulk samples.
- **Composition backend (05b)**: uses `limma::lmFit` plus `eBayes` on logit-transformed proportions, `qlogis(pmin(pmax(x, 1e-5), 1 - 1e-5))`, instead of `speckle::propeller`. This matches propeller's default logit-transformed proportion model without requiring `speckle`. The manifest schema is backend-agnostic, so `run_formal_propeller_st()` can be swapped to `speckle::propeller` later if needed.
- **Spot-level rescue (05_spatial_de.R)**: provides the N=2 exploratory fallback when the replicate gate is closed. It uses `Seurat::FindMarkers(test.use = "wilcox")` and relies on Seurat's internal `p_val_adj` field for multiple-testing adjustment.

## 05a Output Schema

When `status == "ok"`, `05a_region_pseudobulk.R` writes `formal_de.tsv` with edgeR `topTags` columns such as `gene`, `logFC`, `logCPM`, `F`, `PValue`, and `FDR`.

When `status` is anything else, including `skipped_replicate_gate` or `failed_singular_design`, it writes `exploratory_only.tsv` with `status`, `reason`, and `exploratory_fallback`. Downstream readers should branch on either the file name or the `status` column in `pseudobulk_de_manifest.tsv`. The `exploratory_fallback` field points to the matching `05_spatial_de.R` output for the same comparison, layer, and region.

The `formal_results_tsv` column in `composition_manifest.tsv` refers to the formal composition result path regardless of whether the gate passed. When the gate skipped, that file contains a status/reason stub instead of formal statistics; the manifest `status` column distinguishes these cases.
