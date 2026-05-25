# H5AD Export Contract

P-R02-B adds a project-editable contract for H5AD mirror exports. The contract is a TSV so it follows the same metadata style as `samples.tsv`, `sections.tsv`, and module target tables.

`init_project.sh` copies `metadata/h5ad_export_contract.tsv.template` to each project as `metadata/h5ad_export_contract.tsv`. Validators read the project file first and fall back to the repository template.

## Columns

| Column | Meaning |
|---|---|
| `field_path` | Contract path such as `obs.sample_id`, `var.gene_id`, `obsm.X_pca`, `layers.counts`, `uns.module`, or `spatial.tissue_positions`. |
| `field_type` | One of `scalar`, `matrix`, `dict`, or `dataframe`. |
| `required` | `yes`, `no`, or `conditional`. P-R02-B enforces `yes`; later hookup PRs add conditional rules. |
| `dtype` | Expected logical dtype: `str`, `float32`, `int64`, `dict`, `dataframe`, or `-` when not checked. |
| `gate` | `primary` violations fail the primary contract; `exploratory` violations warn only. |
| `modality` | `scrna`, `spatial`, or `both`. |
| `semantics` | Human-readable field purpose. |
| `notes` | Rule notes and future conditional logic. |

## Policy

Primary required fields include sample and condition metadata, primary cell type annotations, gene IDs, counts, export module metadata, and spatial section/image fields for spatial exports. Exploratory fields such as UMAP may warn without blocking export.

`H5AD_CONTRACT_FAIL_ON` controls validator behavior:

- `error`: raise on primary violations.
- `warn`: return structured violations and emit warnings.
- `none`: return structured violations without warning.

Validators return `passed`, `violations`, and `summary`, and can serialize violations into manifest-compatible metadata with `contract_violations_to_manifest()`.
