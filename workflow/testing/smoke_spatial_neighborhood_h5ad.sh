#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/spatial_neighborhood.py
Rscript -e 'parse("workflow/05single_script/spatial/06d_spatial_neighborhood.R")'
grep -q 'SPATIAL_NEIGHBORHOOD_H5AD_MODULE' workflow/05single_script/spatial/06d_spatial_neighborhood.R
grep -q 'region_label' workflow/05single_script/spatial/06d_spatial_neighborhood.R
grep -q 'resolve_spatial_h5ad_path' workflow/04python/spatial_neighborhood.py

echo "smoke_spatial_neighborhood_h5ad_ok"
