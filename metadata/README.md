# Metadata Directory

This directory is reserved for the project metadata architecture.

## Source Of Truth

`analysis_questions.tsv` is the planned human-maintained source of truth for
analysis intent. Generated module input tables should be derived from it.

## Planned Tables

| File | Role |
|---|---|
| `samples.tsv` | scRNA + ST sample design table; `modality=scrna` is the backward-compatible default |
| `sections.tsv` | ST physical section registry keyed by `section_id` |
| `spatial_reference_inventory.tsv` | scRNA reference selection/freezing contract for ST deconvolution and joint stages |
| `h5ad_export_contract.tsv` | H5AD mirror export contract for obs/var/obsm/layers/uns/spatial fields |
| `analysis_questions.tsv` | Tier 1 question schema |
| `comparisons.tsv` | Generated DEG/composition inputs |
| `communication_pairs.tsv` | Generated communication inputs |
| `trajectory_pairs.tsv` | Generated trajectory inputs |
| `scenic_targets.tsv` | Generated regulation inputs |
| `enrichment_targets.tsv` | Generated enrichment inputs |
| `gene_program_targets.tsv` | Generated 05 gene-program availability contract |
| `deconv_pairs.tsv` | Generated ST deconvolution inputs |
| `spatial_pairs.tsv` | Generated ST spatial inputs |

M1-M5 should replace the placeholder tables with generated, validated content.

## ST Intake Schema

`samples.tsv` supports these ST columns in addition to the scRNA columns:

```text
modality, section_id, chip_id, bundle_layout, image_path,
run_spatial, run_deconv, run_joint
```

Use `metadata/templates/samples_spatial.tsv.template` and
`metadata/templates/sections.tsv.template` as the starting point for spatial
projects. `workflow/03stages/validate_metadata.sh` writes the expanded canonical
schema and checks that every spatial `section_id` is registered in
`sections.tsv`.

The ST intake audit writes `reports/intake/spatial_input_inventory.tsv` from
`config/spatial_intake_contract.tsv`. SAW rows are presence-checked in M2 and
marked `true_pending_intake_python`; actual GEF conversion is intentionally
deferred to the later SAW bridge.

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
