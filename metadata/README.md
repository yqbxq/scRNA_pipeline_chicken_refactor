# Metadata Directory

This directory is reserved for the project metadata architecture.

## Source Of Truth

`analysis_questions.tsv` is the planned human-maintained source of truth for
analysis intent. Generated module input tables should be derived from it.

## Planned Tables

| File | Role |
|---|---|
| `samples.tsv` | scRNA sample design table |
| `sections.tsv` | Future spatial transcriptomics section table |
| `analysis_questions.tsv` | Tier 1 question schema |
| `comparisons.tsv` | Generated DEG/composition inputs |
| `communication_pairs.tsv` | Generated communication inputs |
| `trajectory_pairs.tsv` | Generated trajectory inputs |
| `scenic_targets.tsv` | Generated regulation inputs |
| `enrichment_targets.tsv` | Generated enrichment inputs |
| `gene_program_targets.tsv` | Generated 05 gene-program availability contract |
| `deconv_pairs.tsv` | Future generated ST deconvolution inputs |
| `spatial_pairs.tsv` | Future generated ST spatial inputs |

M1-M5 should replace the placeholder tables with generated, validated content.

## `trajectory_pairs.tsv` v4 Schema

`trajectory_pairs.tsv` is generated from `analysis_questions.tsv`; do not edit it
directly. The v4 schema has 17 columns:

```text
trajectory_id, source_question_id, layer_scope, root_group, terminal_group,
condition_split_var, condition_split_values, method, tools_to_run,
methods_extra, outlier_qc_policy, regress_cell_cycle, coarse_label_var,
fine_label_var, split_mode, enabled, notes
```

Default `tools_to_run` values are `slingshot,monocle3,paga_dpt,tradeseq,palantir`
for `method=trajectory` and `scvelo_dynamical,scvelo_stochastic,velocyto,cellrank`
for `method=velocity`. Use `methods_extra` for row-level overrides, for example
`+monocle2`, `-palantir`, or `-cellrank`.
