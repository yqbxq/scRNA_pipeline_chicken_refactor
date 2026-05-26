#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

python3 workflow/04python/07b_liana_consensus.py \
  --mock-run \
  --output "${tmp}/liana_consensus_lr.tsv" \
  --manifest "${tmp}/_manifest.json" \
  --ortholog-lut metadata/ortholog_chicken_human.tsv

test -s "${tmp}/liana_consensus_lr.tsv"
test -s "${tmp}/_manifest.json"
rg -q "lr_axis_id" "${tmp}/liana_consensus_lr.tsv"
rg -q "WNT5A\\|WNT5B\\|FZD3\\|LRP6\\|pGC->rgGC" "${tmp}/liana_consensus_lr.tsv"

python3 - "${tmp}/_manifest.json" <<'PY'
import json
import sys
from pathlib import Path

data = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
assert data["module"] == "07b_liana_consensus"
assert "liana_consensus_tsv" in data["outputs"]
assert data["metrics"]["row_n"] >= 1
print("manifest_ok")
PY

echo "smoke_07b_liana_consensus_ok"
