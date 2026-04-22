#!/usr/bin/env bash
set -euo pipefail

python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir git+https://github.com/jianhuupenn/SpaGCN.git
python -m pip install --no-cache-dir torch-geometric
python -m pip install --no-cache-dir git+https://github.com/QIFEIDKN/STAGATE_pyG.git

python - <<'PY'
import scanpy
import SpaGCN
import STAGATE_pyG
print("scanpy", scanpy.__version__)
print("SpaGCN", getattr(SpaGCN, "__version__", "ok"))
print("STAGATE_pyG", getattr(STAGATE_pyG, "__version__", "ok"))
PY
