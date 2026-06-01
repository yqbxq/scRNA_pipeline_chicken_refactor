# Spatial Deconvolution Compare

`07e_deconvolution_compare.R` reads the method manifests from 07a-07d and writes a comparison manifest plus `recommended_method.txt`.

When two or more methods complete with `status == ok`, 07e loads each method's `spot_celltype_proportions.tsv`, aligns by `deconv_id`, `section`, `spot_id`, and `cell_type`, then records Pearson correlation, RMSE, and Jensen-Shannon distance in `method_comparison_matrix.tsv`. When fewer than two methods finish, the status is `skipped_too_few_methods`, but `recommended_method.txt` is still written so downstream wiring has a stable contract.

Outputs:

| File | Purpose |
| --- | --- |
| `deconv_compare_manifest.tsv` | Single-row compare status. |
| `method_ranking.tsv` | Recommendation ranking based on completion score, consensus Pearson, and primary-method preference. |
| `method_summary.tsv` | Alias of the method ranking contract for downstream consumers. |
| `method_comparison_matrix.tsv` | Pairwise real metric rows for method pairs and cell types. |
| `method_pairwise_metrics.tsv` | Alias of pairwise Pearson/RMSE/JSD metrics. |
| `celltype_method_consistency.tsv` | Cell type-specific method consistency metrics. |
| `spot_method_consistency.tsv` | Spot-level dominant cell type agreement, entropy, and pairwise cosine consistency. |
| `recommended_method_by_celltype.tsv` | Cell type-level recommendation contract. |
| `deconv_evidence_tier.tsv` | Method evidence tier from completion and consensus. |
| `recommended_method.txt` | Method name consumed by post-07 niche derivation. |
