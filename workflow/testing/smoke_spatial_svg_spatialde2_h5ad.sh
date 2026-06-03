#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/spatialde2_svg.py
Rscript -e 'parse("workflow/05single_script/spatial/09a_spatialde2_svg.R")'
grep -q 'resolve_spatial_h5ad_path_strict' workflow/04python/spatialde2_svg.py
grep -q 'spatialde2_manifest.tsv' workflow/05single_script/spatial/09a_spatialde2_svg.R
grep -q 'skipped_no_spatialde2' workflow/04python/spatialde2_svg.py
grep -q 'SPATIALDE2_ALLOW_PROXY' workflow/04python/spatialde2_svg.py

echo "smoke_spatial_svg_spatialde2_h5ad_ok"
