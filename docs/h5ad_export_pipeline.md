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

Business stage hookup is intentionally deferred to P-R02-E.
