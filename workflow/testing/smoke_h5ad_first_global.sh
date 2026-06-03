#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile \
  workflow/04python/10c_scvelo_dynamical.py \
  workflow/04python/decoupler_scrna.py \
  workflow/04python/pyscenic_pipeline.py \
  workflow/04python/spatial_decoupler.py \
  workflow/04python/spatial_pyscenic.py \
  workflow/04python/spatialde2_svg.py \
  workflow/04python/joint_spatial_sidecar.py \
  workflow/04python/write_matrix_h5ad.py \
  workflow/04python/write_synthetic_h5ad.py

Rscript -e 'files <- c("workflow/05single_script/04e_scdesign3_engine.R","workflow/05single_script/04f_scdesign3_finalize.R","workflow/05single_script/spatial/06d_spatial_neighborhood.R","workflow/05single_script/spatial/06e_niche_derivation.R","workflow/05single_script/spatial/07f_deconvolution_validation.R","workflow/05single_script/spatial/08a_spatial_decoupler.R","workflow/05single_script/spatial/08b_spatial_pyscenic.R","workflow/05single_script/spatial/09a_spatialde2_svg.R","workflow/05single_script/spatial/09b_sparkx_svg.R","workflow/05single_script/spatial/09c_svg_consensus_report.R","workflow/05single_script/spatial/helpers/spatial_svg_utils.R"); for (f in files) parse(f)'

grep -q '^layers.spliced' metadata/h5ad_export_contract.tsv
grep -q '^obsm.X_tf_activity' metadata/h5ad_export_contract.tsv
grep -q '^obs.spatial_niche' metadata/h5ad_export_contract.tsv
grep -q 'resolve_scrna_h5ad_path_strict' workflow/04python/10c_scvelo_dynamical.py
grep -q 'spatial_06e_niche' workflow/05single_script/spatial/06e_niche_derivation.R
grep -q 'synthetic_h5ad_manifest_ok' workflow/05single_script/spatial/07f_deconvolution_validation.R
grep -q '09a_spatialde2_svg.R' workflow/03stages/60_run_spatial_extensions.sh
grep -q '09b_sparkx_svg.R' workflow/03stages/60_run_spatial_extensions.sh
grep -q '09c_svg_consensus_report.R' workflow/03stages/60_run_spatial_extensions.sh

echo "smoke_h5ad_first_global_ok"
