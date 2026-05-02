# Regulation Output Schema

08 writes regulation outputs under `results/tables/regulation`, `results/figures/regulation`, `results/checkpoints/regulation`, and `results/manifests`.

## decoupleR

Per layer:

- `results/tables/regulation/decoupler/<layer>/tf_activity.tsv`
- `results/tables/regulation/decoupler/<layer>/pathway_activity.tsv`
- `results/tables/regulation/decoupler/<layer>/network_mapping_summary.tsv`
- `results/figures/regulation/decoupler/<layer>/tf_activity_heatmap.png`
- `results/figures/regulation/decoupler/<layer>/pathway_activity_heatmap.png`
- `results/checkpoints/regulation/decoupler/<layer>/decoupler_object.rds`

Aggregate:

- `results/tables/regulation/decoupler/decoupler_index.tsv`
- `results/tables/regulation/decoupler/resource_fingerprint.tsv`
- `results/tables/regulation/decoupler/network_mapping_summary.tsv`
- `results/tables/regulation/decoupler/tf_activity.tsv`
- `results/tables/regulation/decoupler/pathway_activity.tsv`
- `results/manifests/08e_decoupler/_manifest.json`

`network_mapping_summary.tsv` contains:

| Column | Meaning |
| --- | --- |
| `resource` | `dorothea` or `progeny`. |
| `source_id` | TF or pathway. |
| `human_target_n` | Original human target count. |
| `mapped_chicken_target_n` | Targets mapped to chicken symbols. |
| `mapping_rate` | `mapped_chicken_target_n / human_target_n`. |
| `kept` | Whether the source passed `DECOUPLER_MIN_TARGETS`. |
| `drop_reason` | Empty when kept; otherwise the reason. |

`resource_fingerprint.tsv` contains:

| Column | Meaning |
| --- | --- |
| `resource` | `dorothea` or `progeny`. |
| `source` | `cache` or `package`. |
| `species_origin` | Human prior-network origin. |
| `mapped_species` | Chicken target gene space. |
| `row_n_raw` | Raw prior rows. |
| `row_n_mapped` | Mapped network rows after source filtering. |
| `package_version` | Runtime package versions. |
| `resource_md5` | Fingerprint of the mapped cached resource. |

## 08f Report

08f writes:

- `reports/eda/regulation/report.md`
- `results/tables/regulation/module_status.tsv`
- `results/tables/regulation/resource_fingerprint.tsv`
- `results/tables/regulation/ortholog_mapping_coverage.tsv`
- `results/tables/regulation/scenic_summary.tsv`
- `results/tables/regulation/decoupler_summary.tsv`
- `results/tables/regulation/scenic_decoupler_tf_overlap.tsv`
- `results/tables/regulation/triage.tsv`
- `results/figures/regulation/scenic_decoupler_tf_overlap.png`
- `results/manifests/08f_regulation_eda/_manifest.json`

`module_status.tsv` contains:

| Column | Meaning |
| --- | --- |
| `source_module` | Producing upstream module. |
| `input_name` | Logical input name. |
| `path` | Expected file path. |
| `available` | `yes` or `no`. |
| `status` | `ok`, `not_available`, or `missing_required`. |
| `reason` | Explanation when unavailable. |

08f reads 05/06/07 summaries only as existing evidence. It must not reclassify DEG roles, enrichment sections, or communication eligibility.
