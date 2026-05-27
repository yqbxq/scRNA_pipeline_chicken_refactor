#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
expr <- parse("workflow/05single_script/07a_cellchat.R")
stopifnot(length(expr) > 0)
expr <- parse("workflow/05single_script/helpers/communication_mapping_utils.R")
stopifnot(length(expr) > 0)
RS

rg -q 'method_evidence_class = "hypothesis_only"' workflow/05single_script/07a_cellchat.R
rg -q 'can_be_primary = "no"' workflow/05single_script/07a_cellchat.R
rg -q 'CellChat is hypothesis_only' workflow/05single_script/07a_cellchat.R
rg -q 'method_evidence_class' workflow/05single_script/helpers/communication_mapping_utils.R
rg -q 'hypothesis_only' docs/communication_method_landscape.md
rg -q 'can_be_primary' docs/communication_consensus_schema.md

echo "smoke_07a_cellchat_hypothesis_only_ok"
