#!/usr/bin/env bash
set -euo pipefail

python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir stlearn

python - <<'PY'
import scanpy
import squidpy
import spatialdata
import stlearn
print("scanpy", scanpy.__version__)
print("squidpy", squidpy.__version__)
print("spatialdata", spatialdata.__version__)
print("stlearn", stlearn.__version__)
PY
