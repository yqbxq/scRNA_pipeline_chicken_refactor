#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

Rscript - "${REPO_ROOT}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
repo <- args[[1]]
source(file.path(repo, "workflow/05single_script/helpers/runtime_utils.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/metadata_io.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/communication_consensus_utils.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/commot_signal_utils.R"), encoding = "UTF-8")

consensus <- data.frame(
  lr_axis_id = c("L1|R1|A->B", "L2|R2|A->B"),
  condition_value = c("syf", "syf"),
  cellchat_hit = c(TRUE, TRUE),
  liana_consensus_hit = c(TRUE, TRUE),
  nichenet_hit = c(TRUE, FALSE),
  commot_spatial_hit = c(FALSE, FALSE),
  stringsAsFactors = FALSE
)
commot <- data.frame(
  lr_axis_id = c("L1|R1|A->B", "L2|R2|A->B"),
  signal_score = c(2.1, 0),
  spatial_support = c("yes", "no"),
  status = c("ok_commot_run", "ok_commot_run"),
  stringsAsFactors = FALSE
)
out <- attach_commot_spatial(consensus, commot)
stopifnot(isTRUE(out$commot_spatial_hit[[1]]))
stopifnot(!isTRUE(out$commot_spatial_hit[[2]]))
tiers <- assign_evidence_tier_07d(out, env_thresholds = list(min_methods_for_primary = 2, require_nichenet_for_primary = FALSE))
stopifnot("primary" %in% tiers$evidence_tier)
cat("smoke_07d_commot_join: ok\n")
RS
