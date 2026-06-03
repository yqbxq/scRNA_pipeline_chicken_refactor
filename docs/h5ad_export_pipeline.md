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
| `07e_communication` | `07_communication.sh` | `03d_annotate:annotated_object` | `results/90a_export_h5ad/07e_communication` |
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

## H5AD-First Consumers

P-R02-H promotes only H5AD-native or cross-language modules to H5AD-first execution. R-native DE, enrichment, CellChat, NicheNet, MultiNicheNet, Seurat clustering/annotation, RCTD, TransferData, CARD, and TSV-only evidence/report stages continue to use their existing contracts.

H5AD-first modules now use strict resolvers:

- `resolve_scrna_h5ad_path_strict()` for scRNA Python consumers.
- `resolve_spatial_h5ad_path_strict()` for spatial Python consumers.

Strict resolvers accept only real `.h5ad` files from `results/90a_export_h5ad/<module>/` or an explicit `.h5ad` override. They do not fall back to RDS.

Current H5AD-first and mirror targets:

| Area | Module | Contract |
| --- | --- | --- |
| scRNA velocity | `10c_scvelo_dynamical.py` | `VELOCITY_INPUT_MODE=auto|h5ad|loom`; H5AD mode requires `layers.spliced`, `layers.unspliced`, and `layers.counts`. |
| scRNA regulation | `decoupler_scrna.py`, `pyscenic_pipeline.py` | Reads scRNA H5AD directly and writes activity TSV plus optional result H5AD. |
| scDesign3 | `04e_scdesign3_engine.R` | Keeps RDS/TSV outputs and adds `scdesign3_04e` H5AD mirror manifests for synthetic references. |
| ST neighborhood | `06d_spatial_neighborhood.R` | Prefers `spatial_03_region` H5AD and uses legacy transient conversion only when no H5AD is present. |
| ST niche | `06e_niche_derivation.R` | Computes from TSV/RDS inputs but exports `spatial_06e_niche` H5AD mirror. |
| ST cell2location | `07d_deconvolution_cell2location.R` | Reads scRNA and spatial H5AD through Python helpers and writes `cell2location_result.h5ad`. |
| ST validation | `07f_deconvolution_validation.R` | Requires ok synthetic H5AD manifest plus at least two successful synthetic rerun methods for scDesign3-mode I11 PASS. |
| ST regulation/SVG/joint sidecars | `spatial_decoupler.py`, `spatial_pyscenic.py`, `spatialde2_svg.py`, `joint_spatial_sidecar.py` | Resolve spatial H5AD strictly and write schema-valid manifests when runtime dependencies are unavailable. |
