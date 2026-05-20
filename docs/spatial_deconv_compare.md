# Spatial Deconvolution Compare

`07e_deconvolution_compare.R` reads the method manifests from 07a-07d and writes a comparison manifest plus `recommended_method.txt`.

When two or more methods complete with `status == ok`, 07e records method-pair rows in `method_comparison_matrix.tsv`. When fewer than two methods finish, the status is `skipped_too_few_methods`, but `recommended_method.txt` is still written using `SPATIAL_DECONV_PRIMARY` so downstream wiring has a stable contract.

Outputs:

| File | Purpose |
| --- | --- |
| `deconv_compare_manifest.tsv` | Single-row compare status. |
| `method_ranking.tsv` | Completion-based method ranking. |
| `method_comparison_matrix.tsv` | Pairwise comparison placeholder or completed comparison rows. |
| `recommended_method.txt` | Method name consumed by post-07 niche derivation. |
