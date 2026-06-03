#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/10c_scvelo_dynamical.py
grep -q 'VELOCITY_INPUT_MODE' workflow/04python/10c_scvelo_dynamical.py
grep -q 'layers: spliced' workflow/04python/10c_scvelo_dynamical.py || grep -q 'spliced' workflow/04python/10c_scvelo_dynamical.py
grep -q 'input_h5ad_path' workflow/04python/10c_scvelo_dynamical.py

echo "smoke_scrna_velocity_h5ad_ok"
