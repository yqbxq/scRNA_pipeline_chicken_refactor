# Analysis Questions Schema

`metadata/analysis_questions.tsv` is the Tier 1 source of truth for analysis
intent. Humans edit this file; M3 will fan out active rows into generated Tier 2
module input tables.

The conceptual design covers 87 questions. The committed TSV contains 83
materialized rows because F23-F26 differential communication rows are derived
from split communication questions and must be generated automatically rather
than hand-maintained.

## Columns

| Column | Type | Example | Rule |
|---|---|---|---|
| `question_id` | stable ID | `F02_TC_GC_bidir_split` | Unique; never rename after publication use. |
| `question_zh` | text | `TC和GC大类双向通讯syf_f5_split` | Human-readable summary. |
| `scope` | enum | `panorama` | One of `panorama`, `GC_subcluster`, `TC_subcluster`, `ST_section`. |
| `sender_groups` | DSL | `[TC]` | Primary group expression. |
| `receiver_groups` | DSL | `[pGC]\|[eGC]` | Secondary group expression; use `-` when absent. |
| `condition_split` | split DSL | `group_id:syf,f5` | `-`, `syf_only`, `f5_only`, or `column:value1,value2`. |
| `contrast_axis` | enum | `bidirectional` | Selects fan-out rule. |
| `tools_to_run` | tool list | `cellchat+nichenet` | `+` separated tool names. |
| `priority` | enum | `P0` | One of `P0`, `P1`, `P2`, `P3`. |
| `status` | enum | `active` | One of `active`, `planned`. Planned rows are skipped by semantic checks. |
| `depends_on` | ID list | `H01_GC_trajectory_baseline` | Optional `;` separated upstream question IDs. |
| `notes` | text | `wait_TC_subcluster` | Free text, no tab characters; blank is allowed. |
| `activation_policy` | enum | `auto_if_min_cells` | Tier1 default policy; Tier2 resolved rows are authoritative at runtime. |
| `min_sender_cells` | positive integer | `20` | Minimum sender cells for `auto_if_min_cells`. |
| `min_receiver_cells` | positive integer | `20` | Minimum receiver cells for `auto_if_min_cells`. |
| `min_cells_per_condition` | positive integer | `40` | Minimum `sender_n + receiver_n` after pair subset and condition split. |
| `fallback_pair_id` | ID | `F07_GC_dev_seq_baseline` | Tier1 fallback parent; M3 resolves concrete Tier2 fallback task IDs. |
| `derived_from_pair_id` | ID list | `F08...__pGC_to_eGC` | Tier2 derived rows list concrete split input pair IDs. |
| `run_baseline_if_split_fails` | boolean | `yes` | Report baseline fallback when an auto-gated split is skipped. |
| `display_question_id` | optional ID | `E03_scRNA_GC_TC_capture_balance` | Optional display ID; blank defaults to `question_id`. E03 must set this alias. |
| `output_alias` | optional ID | `E03_scRNA_GC_TC_capture_balance` | Optional output file alias; blank defaults to `question_id`. E03 must set this alias. |
| `report_title` | optional text | `scRNA GC/TC captured-cell balance QC` | Optional report title; blank defaults to an automatic title. E03 must set this title. |

## Scope Semantics

| Scope | Meaning |
|---|---|
| `panorama` | Whole integrated scRNA object; broad GC/TC labels and backfilled subtype labels. |
| `GC_subcluster` | GC-focused subcluster object. |
| `TC_subcluster` | TC-focused subcluster object; most rows remain `planned` until TC annotation exists. |
| `ST_section` | Spatial transcriptomics sections; all current ST rows are `planned`. |

## Group DSL

The validator checks syntax, while M3 will implement full expansion.

| Pattern | Meaning |
|---|---|
| `-` | No sender/receiver for this axis. |
| `*` | All valid groups in the current scope. |
| `[A]` | One named group. |
| `[A,B,C]` | Explicit group set. |
| `[A]\|[B]\|[C]` | Independent fan-out tasks. |
| `[A]->[B]->[C]` | Ordered sequence. Unicode `→` is also accepted. |
| `prefix_*` | Prefix selector, such as `TC_*`. |
| reserved tokens | `all_cells`, `GC_subtypes`, `TC_subtypes`, `regions`, `all_spots`, `auto`, `region`, `synthetic`, `panorama_ref`, `GC_sub_ref`, `TC_sub_ref`, `H01_results`, `F02_results`, `I13_results`, `I19_results`, `I19+I20`, `[from_I05,I06]`. |

`[A,B]` and `[A]|[B]` are intentionally different. `[pGC,eGC]` is one
set-valued expression and is used by fan-out rules such as pairwise comparison.
`[pGC]|[eGC]` means two independent fan-out tasks. For example, `B02` uses
`[pGC,eGC,rgGC,lGC]` so M3 creates all pairwise comparisons from one set, while
`F03` uses `[pGC]|[eGC]|[rgGC]|[lGC]` so M3 creates one TC-to-target
communication task per receiver.

## Condition Split DSL

| Pattern | Meaning |
|---|---|
| `-` | No split. |
| `syf_only` | Subset to syf. |
| `f5_only` | Subset to f5. |
| `group_id:syf,f5` | Split by metadata column `group_id` with values `syf` and `f5`. |
| `section:syf,f5` | Future ST split by section metadata. |

When an annotated object exists, M2 checks that the referenced split column is
present in object metadata. Planned rows are skipped.

## contrast_axis Values

| Value | Generated target | Fan-out intent |
|---|---|---|
| `cluster_marker` | `annotation_marker_targets.tsv` | Raw cluster marker evidence for 03/04 annotation only. |
| `identity_marker` | `comparisons.tsv` | Post-annotation cell type/subtype identity marker, such as A01/A02. |
| `directional_DEG` | `comparisons.tsv` | Directed sender vs receiver DEG. |
| `pairwise` | `comparisons.tsv` | All pairwise combinations in a group set. |
| `contrast_only` | `comparisons.tsv` | Within-group condition contrast. |
| `global_stage_context` | `comparisons.tsv` | Global syf/f5 context signature; not receiver DEG. |
| `composition` | `comparisons.tsv` | Composition/proportion test rows. |
| `qc_composition` | `comparisons.tsv` | QC-only composition/capture-balance rows. |
| `bidirectional` | `communication_pairs.tsv` | Sender to receiver and receiver to sender. |
| `directional` | `communication_pairs.tsv` | Directed sender to receiver tasks. |
| `sequential` | `communication_pairs.tsv` | Adjacent transitions in an ordered chain. |
| `symmetric` | `communication_pairs.tsv` | Full matrix or all-vs-all task. |
| `pairwise_comm` | `communication_pairs.tsv` | Pairwise communication among a set. |
| `regulation_per` | `scenic_targets.tsv` | Per group/layer regulation activity. |
| `regulation_pair` | `scenic_targets.tsv` | Pairwise regulation contrast. |
| `regulation_stage` | `scenic_targets.tsv` | Per group stage regulation contrast. |
| `lineage` | `trajectory_pairs.tsv` | Root-to-terminal trajectory. |
| `velocity` | `trajectory_pairs.tsv` | RNA velocity task. |
| `spatial_clustering` | `spatial_pairs.tsv` | Future ST clustering. |
| `spatial_neighbor` | `spatial_pairs.tsv` | Future spatial neighborhood task. |
| `spatial_overlay` | `spatial_pairs.tsv` | Future ST/scRNA overlay task. |
| `deconv_pair` | `deconv_pairs.tsv` | Future deconvolution task. |
| `SVG` | `spatial_pairs.tsv` | Spatial variable gene task. |
| `SVG_diff` | `spatial_pairs.tsv` | Stage comparison for SVGs. |
| `deconv_validation` | `deconv_pairs.tsv` | Future deconvolution validation. |
| `enrichment_target` | `enrichment_targets.tsv` | Explicit planned enrichment target. |

## Generated Tables

M3 owns these files under `metadata/`; do not hand-edit them after the generator
exists.

| File | ID column | Purpose |
|---|---|---|
| `comparisons.tsv` | `comparison_id` | DEG, composition, and comparison-like tasks. |
| `annotation_marker_targets.tsv` | `target_id` | 03/04 annotation-only raw cluster marker targets. |
| `communication_pairs.tsv` | `pair_id` | CellChat/NicheNet tasks. |
| `trajectory_pairs.tsv` | `trajectory_id` | Trajectory and velocity tasks. |
| `scenic_targets.tsv` | `target_id` | SCENIC/decoupleR regulation tasks. |
| `enrichment_targets.tsv` | `target_id` | GO/KEGG/gProfiler enrichment tasks. |
| `deconv_pairs.tsv` | `deconv_id` | Future ST deconvolution tasks. |
| `spatial_pairs.tsv` | `spatial_pair_id` | Future ST neighborhood/overlay tasks. |

Every generated file starts with:

```text
# AUTO-GENERATED. Edit metadata/analysis_questions.tsv instead.
```

Current generated columns:

| File | Columns |
|---|---|
| `comparisons.tsv` | `comparison_id`, `source_question_id`, `display_question_id`, `output_alias`, `report_title`, `layer_scope`, `contrast_axis`, `analysis_mode`, `analysis_unit`, `stat_level`, `group_var`, `ident_1`, `ident_2`, `subset_column`, `subset_value`, `aggregation_group_var`, `composition_group_var`, `batch_var`, `enabled`, `min_biological_replicates`, `force_exploratory`, `min_cells_per_group`, `logfc_threshold`, `produces_gene_program`, `gene_program_role`, `notes` |
| `annotation_marker_targets.tsv` | `target_id`, `source_question_id`, `layer_scope`, `object_layer`, `cluster_column`, `annotation_label_column`, `group_var`, `ident_1`, `ident_2`, `analysis_mode`, `gene_program_role`, `output_dir`, `annotation_only`, `enabled`, `notes` |
| `communication_pairs.tsv` | `pair_id`, `source_question_id`, `layer_scope`, `sender`, `receiver`, `condition_split_var`, `condition_split_values`, `tool`, `communication_mode`, `activation_policy`, `min_sender_cells`, `min_receiver_cells`, `min_cells_per_condition`, `fallback_pair_id`, `derived_from_pair_id`, `run_baseline_if_split_fails`, `requires_all_derived_inputs_pass`, `receiver_gene_program_source`, `baseline_marker_comparison_id`, `receiver_deg_comparison_id`, `direction_filter`, `requires_cell_subtype`, `notes`, `enabled` |
| `trajectory_pairs.tsv` | `trajectory_id`, `source_question_id`, `layer_scope`, `root_group`, `terminal_group`, `condition_split_var`, `condition_split_values`, `method`, `tools_to_run`, `methods_extra`, `outlier_qc_policy`, `regress_cell_cycle`, `coarse_label_var`, `fine_label_var`, `split_mode`, `enabled`, `notes` |
| `scenic_targets.tsv` | `target_id`, `source_question_id`, `layer_scope`, `cell_subset`, `contrast_axis`, `condition_split_var`, `condition_split_values`, `method`, `enabled`, `notes` |
| `enrichment_targets.tsv` | `target_id`, `source_question_id`, `comparison_id`, `layer_scope`, `analysis_mode`, `gene_program_role`, `organism`, `database`, `enrichment_eligible`, `enrichment_usage`, `min_genes`, `enabled`, `notes` |
| `gene_program_targets.tsv` | `comparison_id`, `source_question_id`, `layer_scope`, `analysis_mode`, `gene_program_role`, `produces_gene_program`, `annotation_only`, `qc_only`, `global_context_only`, `nichenet_eligible`, `nichenet_usage`, `enrichment_eligible`, `enrichment_usage`, `preferred_for_downstream`, `expected_result_level`, `formal_preferred`, `formal_status`, `result_status`, `skip_reason`, `eligible_reason`, `ineligible_reason`, `notes` |
| `deconv_pairs.tsv` | `deconv_id`, `source_question_id`, `st_scope`, `reference_scope`, `section_filter`, `condition_split_var`, `condition_split_values`, `tool`, `enabled`, `notes` |
| `spatial_pairs.tsv` | `spatial_pair_id`, `source_question_id`, `st_scope`, `sender`, `receiver`, `contrast_axis`, `condition_split_var`, `condition_split_values`, `tool`, `enabled`, `notes` |

## Marker / DEG Layering Contract

`analysis_mode` and `gene_program_role` are the authoritative runtime contract.
Modules must not infer downstream use from the presence of a TSV/RDS path.

| `analysis_mode` | `gene_program_role` | Producer | Downstream use |
|---|---|---|---|
| `annotation_cluster_marker` | `annotation_marker` | 03/04 annotation | Annotation evidence only; no 06 mechanism enrichment or 07 NicheNet. |
| `subtype_marker` | `receiver_marker` | 05a | Receiver identity/baseline marker only. |
| `subtype_pairwise` | `subtype_pairwise_deg` | 05a | Subtype pairwise enrichment; not default NicheNet receiver DEG. |
| `condition_within_type` | `condition_deg` | 05b | Only valid `receiver_deg_comparison_id` source for split NicheNet. |
| `composition` | `none` | 05c | Composition result only; no gene program. |
| `qc_composition` | `qc_only` | 05c | QC report only; no gene program or biological abundance conclusion. |
| `global_context` | `global_context` | 05b/05d | Global context enrichment only; not cell-type-specific receiver DEG. |

A01/A02 are post-annotation identity markers. They no longer use raw
`cluster_marker` semantics. E03 is `qc_composition`; it must not be NicheNet or
enrichment eligible. D05 is `global_context`; it must not be used as
`receiver_deg_comparison_id`.

D05 may produce enrichment results, but those results are contextual background
only and must not be used as cell-type-specific mechanism evidence. D05 可以产生
富集结果，但这些结果只能作为全局背景参考，不能作为具体细胞类型机制证据。

E03 keeps canonical `question_id=E03_layer_compo` for compatibility, but its
display/output alias is `E03_scRNA_GC_TC_capture_balance`. 05c keeps the normal
composition manifest output and additionally writes:

- `results/tables/composition/E03_layer_compo__composition_GC_TC.tsv`
- `results/tables/qc/composition/E03_scRNA_GC_TC_capture_balance.tsv`
- `reports/qc/composition/E03_scRNA_GC_TC_capture_balance.md`

The E03 QC report must state that the result only describes captured scRNA
GC/TC cell composition. It cannot be interpreted as true tissue GC/TC
proportion change. Tissue abundance must be evaluated by ST region annotation,
ST deconvolution, or histology/image quantification.

07b NicheNet source rules:

- Baseline communication uses `baseline_marker_comparison_id`, which must point
  to `subtype_marker + receiver_marker` with `nichenet_usage=baseline_receiver_marker`.
- Condition-split communication uses `receiver_deg_comparison_id`, which must
  point to `condition_within_type + condition_deg` with
  `nichenet_usage=receiver_condition_deg`.
- `annotation_marker`, `subtype_pairwise_deg`, `qc_only`, `global_context`, and
  `none` are forbidden as receiver DEG sources.

06 enrichment report sections are assigned from resolved metadata rather than
from file presence: `condition_deg`/`mechanism_enrichment` rows enter Core
Mechanism Enrichment, marker/pairwise rows enter Subtype / Identity Enrichment,
`global_context`/`global_context_enrichment` rows enter Global Context
Enrichment, and QC/non-gene-program rows enter QC / Non-gene-program.

## M3 Fan-Out Contract

M3 reads all 83 materialized rows and processes active rows for executable
scRNA modules. A small set of planned communication rows is emitted with
`enabled=no` so the final communication schema stays stable. The 24 ST rows are
all planned, so M3 must write empty `deconv_pairs.tsv` and `spatial_pairs.tsv`
with valid headers until ST-H activates and implements ST fan-out.

`contrast_axis` is the dispatch key:

- `annotation_marker_targets.tsv`: `cluster_marker`; raw cluster marker evidence is produced by 03/04 annotation modules
- `comparisons.tsv`: `identity_marker`, `directional_DEG`, `pairwise`, `contrast_only`, `global_stage_context`, `composition`, `qc_composition`
- `communication_pairs.tsv`: `bidirectional`, `directional`, `sequential`, `symmetric`, `pairwise_comm`
- `scenic_targets.tsv`: `regulation_per`, `regulation_pair`, `regulation_stage`
- `trajectory_pairs.tsv`: `lineage`, `velocity`
- `enrichment_targets.tsv`: derived from generated comparisons whose source question has `tools_to_run` containing `enrichment`
- `gene_program_targets.tsv`: derived from generated comparisons and annotation targets; includes positive and explicitly ineligible rows so downstream modules do not infer usage from file paths
- `deconv_pairs.tsv` / `spatial_pairs.tsv`: ST-H only; current implementation is schema-only

Reserved tokens such as `all_cells`, `GC_subtypes`, `TC_subtypes`, `regions`,
`all_spots`, `auto`, `panorama_ref`, `GC_sub_ref`, `TC_sub_ref`,
`H01_results`, `F02_results`, `I13_results`, `I19_results`, `I19+I20`, and
`[from_I05,I06]` are not ordinary cell subtype names. M3 either expands them by
scope-specific rule or passes them to the future ST fan-out path.

M2 semantic checks look for concrete subtype tokens in `cell_subtype` when an
annotated object exposes that column. M5 provides `cell_subtype` by defaulting
panorama cells to broad `cell_type` labels and overwriting subcluster cells from
04b annotations before subtype-level sender/receiver rows such as
pGC/eGC/rgGC/lGC are treated as fully semantic-checked runtime inputs.

For 05 and 07, generated schemas are runtime contracts:

- `__rest__` is valid only for `annotation_cluster_marker` and
  `subtype_marker`.
- `group_var=cell_type` is checked against object `cell_type`; `group_var=cell_subtype`
  is checked against object `cell_subtype`.
- NicheNet-capable rows must reference explicit gene-program comparison IDs in
  `baseline_marker_comparison_id` and, for condition split tasks,
  `receiver_deg_comparison_id`.
- `direction_filter=no` keeps full CellChat networks; `direction_filter=yes`
  keeps sender-to-receiver LR rows for consensus.
- `requires_cell_subtype=yes` is strict and cannot fall back to broad labels.
- 07 runtime consumes only resolved Tier2 `communication_pairs.tsv` policy
  fields. Tier1 policy fields are defaults for M3 fan-out and are never
  interpreted directly by `07a_cellchat.R` or `07b_nichenet.R`.
- `auto_if_min_cells` gates use pair-specific cells after subset, condition
  split, and sender/receiver resolution: `condition_pair_cell_n = sender_n +
  receiver_n`.
- Split fallback is by reference only. Skipped split rows keep method result
  RDS paths empty, set `success=false`, `result_copied=no`, and report the
  baseline through separate fallback fields.

## Derived Communication Rows

The following conceptual rows are not written to
`metadata/analysis_questions.tsv`:

| Derived ID | Source |
|---|---|
| `F23_TC_GC_diff_overall` | Derived from `F02_TC_GC_bidir_split`. |
| `F24_TC_GC_diff_subtype_*` | Derived from `F04_TC_to_GC_each_split` and `F06_GC_each_to_TC_split`; generated as directional subtype subrows. |
| `F25_GC_internal_diff` | Derived from `F08_GC_dev_seq_split`. |
| `F26_panorama_screen_diff` | Derived from `F15_panorama_screen_split`. |

M3 should create these downstream differential summaries automatically for split
communication questions. `F25_GC_internal_diff` uses
`requires_all_derived_inputs_pass=yes`: partial `F08` split evidence may be
reported as partial, but it cannot produce a GC internal differential
conclusion.

## M2 Validation Contract

`workflow/05single_script/96_validate_metadata_consistency.R` validates:

- required columns and enum values,
- unique `question_id`,
- DSL syntax,
- ST rows remain `planned`,
- every active `contrast_axis` has a known Tier 2 target mapping,
- generated table ID uniqueness when those tables exist,
- generated table schema checks for 05/06/07 Tier 2 contracts,
- 05 `analysis_mode`, `gene_program_role`, `__rest__`, condition, composition,
  QC composition, global context, and group-var semantic checks,
- gene program eligibility matrix checks for annotation-only, QC-only, global
  context, NicheNet usage, and enrichment usage,
- 07 explicit `baseline_marker_comparison_id` / `receiver_deg_comparison_id`
  role checks,
- optional semantic checks against annotated RDS metadata when available,
- hard failure when `requires_cell_subtype=yes` rows lack `cell_subtype`,
- NicheNet receiver/DEG consistency when generated comparisons exist.

The validator writes a machine-readable report to:

```text
results/00_validation/metadata_validation.json
```
