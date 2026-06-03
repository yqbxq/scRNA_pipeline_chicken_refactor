#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/decoupler_scrna.py workflow/04python/pyscenic_pipeline.py
grep -q 'resolve_scrna_h5ad_path_strict' workflow/04python/decoupler_scrna.py
grep -q 'regulation_result.h5ad' workflow/04python/pyscenic_pipeline.py

echo "smoke_scrna_regulation_h5ad_ok"
