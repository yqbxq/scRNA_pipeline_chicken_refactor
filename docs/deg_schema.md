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

06, 07, and future 08 must propagate `formal`, `exploratory_only`, or
`exploratory_forced` into manifests and plot/report labels.
