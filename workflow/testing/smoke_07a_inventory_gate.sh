#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/runtime_utils.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/inventory_gate_utils.R", encoding = "UTF-8")

eligibility <- data.frame(
  layer_id = c("panorama", "panorama", "panorama", "panorama"),
  pair_id = "",
  comparison_id = c("cmp1", "cmp1", "cmp2", "cmp2"),
  cluster_id = c("A", "B", "A", "B"),
  evidence_tier = c("primary", "candidate_only", "candidate_only", "skip"),
  stringsAsFactors = FALSE
)

Sys.setenv(INVENTORY_GATE_STRICT_MODE = "no")
gate <- inventory_gate_log_for_labels(c("A", "B"), eligibility, "panorama", "pair1")
stopifnot(gate$passed_inventory_gate[gate$cluster_id == "A"])
stopifnot(!gate$passed_inventory_gate[gate$cluster_id == "B"])

Sys.setenv(INVENTORY_GATE_STRICT_MODE = "yes")
gate_strict <- inventory_gate_log_for_labels(c("A", "B"), eligibility, "panorama", "pair1")
stopifnot(!gate_strict$passed_inventory_gate[gate_strict$cluster_id == "A"])
stopifnot(!gate_strict$passed_inventory_gate[gate_strict$cluster_id == "B"])
RS

rg -q "cellchat_inventory_gate_log_tsv" workflow/05single_script/07a_cellchat.R
rg -q "nichenet_inventory_gate_log_tsv" workflow/05single_script/07c_nichenet.R
rg -q "cluster_eligibility_tsv" workflow/03stages/07_communication.sh
rg -q "INVENTORY_GATE_ALLOWED_TIERS_COMMUNICATION" workflow/02lib/common.sh

echo "smoke_07a_inventory_gate_ok"
