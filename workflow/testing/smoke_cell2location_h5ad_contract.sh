#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/cell2location_pipeline.py
grep -q 'resolve_scrna_h5ad_path' workflow/04python/cell2location_pipeline.py
grep -q 'resolve_spatial_h5ad_path' workflow/04python/cell2location_pipeline.py
grep -q 'cell2location_result.h5ad' workflow/04python/cell2location_pipeline.py
grep -q 'skipped_no_h5ad' workflow/05single_script/spatial/07d_deconvolution_cell2location.R
grep -q 'skipped_no_cell2location' workflow/05single_script/spatial/07d_deconvolution_cell2location.R

echo "smoke_cell2location_h5ad_contract_ok"
