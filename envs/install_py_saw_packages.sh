#!/usr/bin/env bash
set -euo pipefail

python -m pip install --upgrade pip setuptools wheel
python -m pip install --no-cache-dir stereopy

python - <<'PY'
import importlib.metadata as metadata

print("stereopy", metadata.version("stereopy"))
PY
