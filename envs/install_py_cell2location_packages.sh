#!/usr/bin/env bash
set -euo pipefail

python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir scvi-tools
python -m pip install --no-cache-dir cell2location

python - <<'PY'
import scanpy
import scvi
import cell2location
print("scanpy", scanpy.__version__)
print("scvi-tools", scvi.__version__)
print("cell2location", getattr(cell2location, "__version__", "ok"))
PY
