# H5AD Export Pipeline

P-R02-C adds an isolated `90_export_h5ad` stage. It mirrors existing RDS source-of-truth objects into H5AD files for Python tools without changing any business stage.

## Stage Interface

```bash
bash workflow/03stages/90_export_h5ad.sh \
  --upstream-manifest results/manifests/03d_annotate/_manifest.json \
  --output-key annotated_object \
  --modality scrna \
  --target-dir results/90a_export_h5ad/03d_annotate
```

The stage resolves `--output-key` from the upstream manifest, selects the scRNA or spatial exporter, validates the object against `metadata/h5ad_export_contract.tsv`, writes an export summary TSV, and writes its own `_manifest.json`.

## Mapping

| Source | H5AD target |
|---|---|
| Seurat `meta.data` | `obs` |
| assay feature metadata | `var` |
| reductions | `obsm` using `X_<reduction>` names, with `spatial` kept as `spatial` |
| counts/data/SCT slots or mock layers | `layers` |
| module metadata, timestamps, fingerprints, contract results | `uns` |
| spatial coordinates and scale factors | `uns["spatial"]` and `obsm["spatial"]` |

Spatial exports default to one H5AD per `section_id` because squidpy and COMMOT usually build spatial graphs section by section.

## Backends

`H5AD_EXPORT_BACKEND` defaults to `zellkonverter`. `sceasy` is retained as a fallback option. The smoke test uses `mock` so local validation does not depend on Seurat, zellkonverter, or anndata.

## Business Hookup

P-R02-E connects the `90_export_h5ad` mirror to business stages without changing their primary RDS outputs. Hooks use `export_h5ad_for_gate` from `workflow/02lib/common.sh`.

Default hooks run when `H5AD_EXPORT_GATES` is empty:

| Gate | Stage | Source manifest output | Target directory |
|---|---|---|---|
| `03d_annotation` | `03_panorama.sh` | `03d_annotate:annotated_object` | `results/90a_export_h5ad/03d_panorama` |
| `04b_subcluster` | `04_subcluster.sh` | `04b_subcluster_annotate:panorama_cell_subtype_rds` | `results/90a_export_h5ad/04b_subcluster` |
| `06d_enrichment` | `06_enrichment.sh` | `03d_annotate:annotated_object` | `results/90a_export_h5ad/06d_enrichment` |
| `07c_communication` | `07_communication.sh` | `03d_annotate:annotated_object` | `results/90a_export_h5ad/07c_communication` |
| `spatial_03_region` | `50_run_spatial_pipeline.sh` | `spatial_03_region_annotation:panorama_annotated` | `results/90a_export_h5ad/spatial_03_region` |
| `spatial_04b_subcluster` | `50_run_spatial_pipeline.sh` | `spatial_04b_subcluster_annotate:panorama_subannotated` | `results/90a_export_h5ad/spatial_04b_subcluster` |
| `spatial_05_marker` | `50_run_spatial_pipeline.sh` | `spatial_03_region_annotation:panorama_annotated` | `results/90a_export_h5ad/spatial_05_marker` |
| `spatial_06_extensions` | `60_run_spatial_extensions.sh` | `spatial_06e_niche:panorama_niched` | `results/90a_export_h5ad/spatial_06_extensions` |

Set `H5AD_EXPORT_GATES=03d_annotation,spatial_03_region` to export only selected gates. `H5AD_EXPORT_ON_FAILURE=skip` is the default, so mirror export failures do not block the business stage; use `warn` for louder logs or `error` when a downstream Python workflow requires H5AD.

## Python Consumers

Python-side helpers now live under `workflow/04python/helpers/`:

- `scrna_io.py`: `find_scrna_h5ad`, `resolve_scrna_h5ad_path`, `load_scrna_h5ad`
- `spatial_io.py`: `find_spatial_h5ad`, `resolve_spatial_h5ad_path`, `load_spatial_h5ad`, `load_all_spatial_h5ad`

Consumers should resolve H5AD through these helpers instead of reconstructing AnnData from Seurat RDS. During the transition, `H5AD_PYTHON_FALLBACK_TO_RDS=yes` allows explicit legacy paths to remain usable and emits a deprecation warning.
