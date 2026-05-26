#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

fixtures=$'WNT5A\tFZD3\tpGC\trgGC\tWNT5A|FZD3|pGC->rgGC\nWNT5B_WNT5A\tFZD3+LRP6\tpGC\trgGC\tWNT5A|WNT5B|FZD3|LRP6|pGC->rgGC\nTGFB1\tTGFBR2/TGFBR1\tTC\tGC\tTGFB1|TGFBR1|TGFBR2|TC->GC'

while IFS=$'\t' read -r ligand receptor sender receiver expected; do
  [[ -n "${ligand}" ]] || continue
  r_out="$(Rscript -e 'source("workflow/05single_script/helpers/ortholog_lookup_utils.R"); source("workflow/05single_script/helpers/communication_consensus_utils.R"); args <- commandArgs(TRUE); cat(standardize_lr_axis_id(args[1], args[2], args[3], args[4]))' "${ligand}" "${receptor}" "${sender}" "${receiver}")"
  py_out="$(PYTHONPATH="${ROOT_DIR}/workflow/04python" python3 - "${ligand}" "${receptor}" "${sender}" "${receiver}" <<'PY'
import sys
from helpers.lr_axis_id_utils import standardize_lr_axis_id
print(standardize_lr_axis_id(*sys.argv[1:5]), end="")
PY
)"
  [[ "${r_out}" == "${expected}" ]] || { echo "R mismatch: ${r_out} != ${expected}" >&2; exit 1; }
  [[ "${py_out}" == "${expected}" ]] || { echo "Python mismatch: ${py_out} != ${expected}" >&2; exit 1; }
done <<< "${fixtures}"

echo "smoke_standardize_lr_axis_id_ok"
