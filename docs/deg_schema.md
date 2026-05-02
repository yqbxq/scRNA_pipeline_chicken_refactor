# DEG / Gene Program Schema

`metadata/comparisons.tsv` is generated from `metadata/analysis_questions.tsv`.
05 uses `analysis_mode` to decide the execution path:

| `analysis_mode` | Main use | Formal status |
|---|---|---|
| `annotation_cluster_marker` | 03/04 raw cluster annotation evidence only. | Annotation-only; not executed by 05a. |
| `subtype_marker` | Receiver subtype or broad cell-type marker, usually one-vs-rest. | Cell-level exploratory. |
| `subtype_pairwise` | Subtype identity comparison. | Cell-level exploratory. |
| `condition_within_type` | syf-vs-f5 within a cell type/subtype. | Formal only when pseudobulk replicate gate passes. |
| `composition` | Sample-level proportion test. | Not a gene program. |
| `qc_composition` | QC-only capture/composition check. | Not a gene program. |
| `global_context` | Whole-layer/global syf-vs-f5 context signature. | Contextual only; not NicheNet receiver DEG. |

D05 is the canonical `global_context` row. It may produce enrichment results,
but those results are contextual background only and must not be used as
cell-type-specific mechanism evidence. D05 可以产生富集结果，但这些结果只能作为
全局背景参考，不能作为具体细胞类型机制证据。

E03 is the canonical QC composition row. It keeps `question_id=E03_layer_compo`
but uses `output_alias=E03_scRNA_GC_TC_capture_balance` for the QC mirror table
and report. E03 is a scRNA captured-cell balance QC result, not a spatial or
tissue abundance estimate. Do not use E03 to claim GC/TC tissue proportion
changes; use ST region annotation, ST deconvolution, or histology/image
quantification for tissue abundance.

`__rest__` is only valid for `annotation_cluster_marker` and `subtype_marker`.
Rows using `group_var=cell_type` are checked against `cell_type`; rows using
`group_var=cell_subtype` are checked against `cell_subtype`.

05d writes `results/tables/deg/gene_program_registry.tsv`. Downstream modules
must use this registry instead of guessing paths from `pair_id` or
`comparison_id` naming. Key fields:

| Column | Meaning |
|---|---|
| `comparison_id` | Stable ID from `comparisons.tsv` / `gene_program_targets.tsv`. |
| `result_level` | `pseudobulk_formal`, `cell_level_exploratory`, or `unavailable`. |
| `formal_status` | Raw formal inference status from 05b when available. |
| `nichenet_eligible` / `nichenet_usage` | Whether and how 07b may use the row. |
| `enrichment_eligible` / `enrichment_usage` | Whether and where 06 may enrich the row. |
| `annotation_only`, `qc_only`, `global_context_only` | Hard separation flags for non-mechanism outputs. |
| `deg_tsv` | Formal pseudobulk DEG table path when available. |
| `marker_tsv` | Cell-level marker / exploratory table path. |
| `top_gene_tsv` | Preferred table for downstream use. |
| `warning` | Explains exploratory fallback, especially when formal was preferred. |

For E03, registry rows must remain explicitly ineligible for downstream gene
program use: `result_status=qc_only`,
`skip_reason=qc_composition_does_not_produce_gene_program`,
`preferred_for_downstream=no`, and `deg_tsv` / `marker_tsv` / `top_gene_tsv`
must be empty. `composition_tsv` may point to the QC mirror table.

06, 07, and future 08 must propagate `formal`, `exploratory_only`, or
`exploratory_forced` into manifests and plot/report labels.

06c must place `global_context` rows such as D05 in the Global Context
Enrichment section, not Core Mechanism Enrichment. 07 must reject any
`receiver_deg_comparison_id` that points to `global_context`.

A5 hardens the registry-only contract: 06 and 07 must not run inline
`FindMarkers`, scan result directories, infer paths from IDs, or fall back across
registry path columns. They may consume only the explicit `top_gene_tsv` selected
by `gene_program_registry.tsv`; missing registry rows or missing files are
reported as `missing_gene_program`.
