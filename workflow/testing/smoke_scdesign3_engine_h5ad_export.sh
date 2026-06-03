#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/write_matrix_h5ad.py
Rscript -e 'parse("workflow/05single_script/04e_scdesign3_engine.R"); parse("workflow/05single_script/04f_scdesign3_finalize.R")'
grep -q 'scdesign3_h5ad_manifest.tsv' workflow/05single_script/04e_scdesign3_engine.R
grep -q 'write_matrix_h5ad.py' workflow/05single_script/04e_scdesign3_engine.R
grep -q 'engine_h5ad_manifest' workflow/05single_script/04f_scdesign3_finalize.R

echo "smoke_scdesign3_engine_h5ad_export_ok"
