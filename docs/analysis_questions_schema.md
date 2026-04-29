# Analysis Questions Schema

`metadata/analysis_questions.tsv` is the planned Tier 1 source of truth for
analysis intent. Each row describes one biological question. M3 will fan out
active rows into generated Tier 2 module input tables.

## Columns

| Column | Meaning |
|---|---|
| `question_id` | Stable question ID used in reports and generated row IDs |
| `question_zh` | Human-readable question text |
| `scope` | Analysis scope, such as `panorama`, `GC_subcluster`, `TC_subcluster`, or `ST_section` |
| `sender_groups` | DSL expression for sender or primary groups |
| `receiver_groups` | DSL expression for receiver or secondary groups |
| `condition_split` | Split condition, such as `group:syf,f5`, `syf_only`, `f5_only`, or `-` |
| `contrast_axis` | Fan-out rule selector |
| `tools_to_run` | Tool set, such as `deg+enrichment` or `cellchat+nichenet` |
| `priority` | `P0`, `P1`, `P2`, or `P3` |
| `status` | `active` or `planned` |
| `depends_on` | Optional upstream `question_id` |
| `notes` | Free text |

## DSL Notes

The planned parser should support:

- `*` for all valid groups in the current scope
- `[A]` for one group
- `[A,B,C]` for a set
- `[A]->[B]->[C]` for sequential transitions
- `[A]|[B]|[C]` for independent fan-out tasks
- Prefix selectors such as `TC_*` where supported by the scope metadata

## Generated Tables

M3 should write the generated tables under `metadata/`:

- `comparisons.tsv`
- `communication_pairs.tsv`
- `trajectory_pairs.tsv`
- `scenic_targets.tsv`
- `enrichment_targets.tsv`
- `deconv_pairs.tsv`
- `spatial_pairs.tsv`

Generated files must keep an `AUTO-GENERATED` marker and should not be edited
manually.
