# Analysis Questions Schema

`metadata/analysis_questions.tsv` is the Tier 1 source of truth for analysis
intent. Humans edit this file; M3 will fan out active rows into generated Tier 2
module input tables.

The conceptual design covers 84 questions. The committed TSV contains 80
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
| `cluster_marker` | `comparisons.tsv` | Each group vs all others. |
| `directional_DEG` | `comparisons.tsv` | Directed sender vs receiver DEG. |
| `pairwise` | `comparisons.tsv` | All pairwise combinations in a group set. |
| `contrast_only` | `comparisons.tsv` | Within-group condition contrast. |
| `composition` | `comparisons.tsv` | Composition/proportion test rows. |
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
| `comparisons.tsv` | `comparison_id`, `source_question_id`, `layer_scope`, `contrast_axis`, `ident_1`, `ident_2`, `subset_column`, `subset_value`, `group_var`, `batch_var`, `enabled`, `min_biological_replicates`, `force_exploratory`, `min_cells_per_group`, `logfc_threshold`, `notes` |
| `communication_pairs.tsv` | `pair_id`, `source_question_id`, `layer_scope`, `sender`, `receiver`, `subset_column`, `subset_value`, `condition_split_var`, `condition_split_values`, `tool`, `enabled`, `notes` |
| `trajectory_pairs.tsv` | `trajectory_id`, `source_question_id`, `layer_scope`, `root_group`, `terminal_group`, `condition_split_var`, `condition_split_values`, `method`, `enabled`, `notes` |
| `scenic_targets.tsv` | `target_id`, `source_question_id`, `layer_scope`, `cell_subset`, `contrast_axis`, `condition_split_var`, `condition_split_values`, `method`, `enabled`, `notes` |
| `enrichment_targets.tsv` | `target_id`, `source_question_id`, `comparison_id`, `layer_scope`, `organism`, `database`, `min_genes`, `enabled`, `notes` |
| `deconv_pairs.tsv` | `deconv_id`, `source_question_id`, `st_scope`, `reference_scope`, `section_filter`, `condition_split_var`, `condition_split_values`, `tool`, `enabled`, `notes` |
| `spatial_pairs.tsv` | `spatial_pair_id`, `source_question_id`, `st_scope`, `sender`, `receiver`, `contrast_axis`, `condition_split_var`, `condition_split_values`, `tool`, `enabled`, `notes` |

## M3 Fan-Out Contract

M3 reads all 80 materialized rows but processes only `status=active`. Current
counts are 38 active and 42 planned rows. The 24 ST rows are all planned, so M3
must write empty `deconv_pairs.tsv` and `spatial_pairs.tsv` with valid headers
until ST-H activates and implements ST fan-out.

`contrast_axis` is the dispatch key:

- `comparisons.tsv`: `cluster_marker`, `directional_DEG`, `pairwise`, `contrast_only`, `composition`
- `communication_pairs.tsv`: `bidirectional`, `directional`, `sequential`, `symmetric`, `pairwise_comm`
- `scenic_targets.tsv`: `regulation_per`, `regulation_pair`, `regulation_stage`
- `trajectory_pairs.tsv`: `lineage`, `velocity`
- `enrichment_targets.tsv`: derived from generated comparisons whose source question has `tools_to_run` containing `enrichment`
- `deconv_pairs.tsv` / `spatial_pairs.tsv`: ST-H only; current implementation is schema-only

Reserved tokens such as `all_cells`, `GC_subtypes`, `TC_subtypes`, `regions`,
`all_spots`, `auto`, `panorama_ref`, `GC_sub_ref`, `TC_sub_ref`,
`H01_results`, `F02_results`, `I13_results`, `I19_results`, `I19+I20`, and
`[from_I05,I06]` are not ordinary cell subtype names. M3 either expands them by
scope-specific rule or passes them to the future ST fan-out path.

M2 semantic checks look for concrete subtype tokens in `cell_subtype` when an
annotated object exposes that column. M5 panorama backfill must provide
`cell_subtype` before subtype-level sender/receiver rows such as pGC/eGC/rgGC/lGC
are treated as fully semantic-checked runtime inputs.

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
communication questions.

## M2 Validation Contract

`workflow/05single_script/96_validate_metadata_consistency.R` validates:

- required columns and enum values,
- unique `question_id`,
- DSL syntax,
- ST rows remain `planned`,
- every active `contrast_axis` has a known Tier 2 target mapping,
- generated table ID uniqueness when those tables exist,
- optional semantic checks against annotated RDS metadata when available,
- NicheNet receiver/DEG consistency when generated comparisons exist.

The validator writes a machine-readable report to:

```text
results/00_validation/metadata_validation.json
```
