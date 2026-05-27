#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
invisible(parse("workflow/05single_script/07d_communication_consensus.R"))
invisible(parse("workflow/05single_script/helpers/communication_consensus_utils.R"))
source("workflow/05single_script/helpers/communication_consensus_utils.R", encoding = "UTF-8")

cellchat <- data.frame(
  lr_axis_id = c("L1|R1|A->B", "L2|R2|A->C"),
  condition_value = c("syf", "syf"),
  source = c("A", "A"),
  target = c("B", "C"),
  ligand = c("L1", "L2"),
  receptor = c("R1", "R2"),
  method_evidence_class = "hypothesis_only",
  can_be_primary = "no",
  stringsAsFactors = FALSE
)
liana <- data.frame(
  lr_axis_id = c("L1|R1|A->B", "L3|R3|D->E"),
  condition_value = c("syf", "syf"),
  source = c("A", "D"),
  target = c("B", "E"),
  ligand = c("L1", "L3"),
  receptor = c("R1", "R3"),
  liana_consensus_hit = c(TRUE, TRUE),
  n_methods_agreed = c(4, 4),
  liana_consensus_score = c(0.01, 0.02),
  stringsAsFactors = FALSE
)
nichenet <- data.frame(
  lr_axis_id = "L1|R1|A->B",
  condition_value = "syf",
  source = "A",
  target = "B",
  ligand = "L1",
  receptor = "R1",
  nichenet_hit = TRUE,
  method_evidence_class = "downstream_validation",
  stringsAsFactors = FALSE
)

joined <- full_outer_join_lr_tables(cellchat, liana, nichenet)
tiers <- assign_evidence_tier_07d(joined)

primary <- tiers[tiers$lr_axis_id == "L1|R1|A->B", , drop = FALSE]
candidate <- tiers[tiers$lr_axis_id == "L2|R2|A->C", , drop = FALSE]
exploratory <- tiers[tiers$lr_axis_id == "L3|R3|D->E", , drop = FALSE]
stopifnot(primary$evidence_tier[[1]] == "primary")
stopifnot(primary$can_be_primary[[1]] == "yes")
stopifnot(candidate$evidence_tier[[1]] == "candidate")
stopifnot(candidate$can_be_primary[[1]] == "no")
stopifnot(exploratory$evidence_tier[[1]] == "exploratory")
stopifnot(nrow(assign_evidence_tier_07d(full_outer_join_lr_tables(data.frame(), data.frame(), data.frame()))) == 0)
only_liana <- assign_evidence_tier_07d(full_outer_join_lr_tables(data.frame(), liana[1, , drop = FALSE], data.frame()))
stopifnot(only_liana$evidence_tier[[1]] == "exploratory")
RS

rg -q "method_consensus_tsv" workflow/05single_script/07d_communication_consensus.R
rg -q "communication_evidence_tier_tsv" workflow/05single_script/07d_communication_consensus.R
rg -q "CONSENSUS_REQUIRE_NICHENET_FOR_PRIMARY" workflow/02lib/common.sh docs/env_registry.md

echo "smoke_07d_consensus_layer_ok"
