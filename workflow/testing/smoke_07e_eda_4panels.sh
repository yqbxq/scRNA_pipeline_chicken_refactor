#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

Rscript - "${REPO_ROOT}" "${TMP_DIR}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
repo <- args[[1]]
tmp <- args[[2]]
source(file.path(repo, "workflow/05single_script/helpers/runtime_utils.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/metadata_io.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/communication_report_utils.R"), encoding = "UTF-8")
source(file.path(repo, "workflow/05single_script/helpers/communication_report_panels.R"), encoding = "UTF-8")

cfg <- list(
  project_root = tmp,
  table_dir = file.path(tmp, "results/tables"),
  communication_table_dir = file.path(tmp, "results/tables/communication"),
  communication_consensus_figure_dir = file.path(tmp, "results/figures/communication/consensus"),
  communication_eda_panel_dir = file.path(tmp, "results/figures/communication/consensus/07e_panels"),
  communication_report_filter_tier = "auto",
  communication_report_top_n_primary = 30L,
  communication_fallback_summary_tsv = file.path(tmp, "missing_fallback.tsv"),
  nichenet_index_tsv = file.path(tmp, "missing_nichenet.tsv")
)

consensus <- data.frame(
  lr_axis_id = paste0("L", 1:5, "|R", 1:5, "|A->B"),
  condition_value = c("syf", "syf", "f5", "f5", "all"),
  pair_id = c("P1", "P1", "P2", "P2", "P3"),
  layer_id = "panorama",
  source = "TC",
  target = "GC",
  ligand = paste0("L", 1:5),
  receptor = paste0("R", 1:5),
  cellchat_hit = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  liana_consensus_hit = c(TRUE, TRUE, FALSE, TRUE, FALSE),
  nichenet_hit = c(TRUE, FALSE, FALSE, FALSE, FALSE),
  commot_spatial_hit = c(FALSE, FALSE, FALSE, FALSE, FALSE),
  evidence_method_n = c(3, 2, 1, 1, 0),
  evidence_tier = c("primary", "exploratory", "candidate", "exploratory", "blocked"),
  can_be_primary = c("yes", "no", "no", "no", "no"),
  nichenet_n_targets_in_receiver_de = c(12, 0, 0, 0, 0),
  stringsAsFactors = FALSE
)
pairs <- data.frame(
  pair_id = c("P1", "P2", "P3"),
  methods_required = c("cellchat,liana,nichenet", "cellchat,liana", "cellchat"),
  min_methods_agreed = c("2", "2", "1"),
  require_downstream_de = c("auto", "auto", "no"),
  require_spatial_support = c("auto", "auto", "no"),
  evidence_tier_required = c("exploratory", "candidate", "blocked"),
  cellchat_min_samples_consistent = c("", "", ""),
  stringsAsFactors = FALSE
)

filtered <- communication_report_filter_by_tier(consensus, pairs, "auto")
stopifnot(nrow(filtered) == 5)
panels <- list(
  panel_01_tier_distribution = render_panel_tier_distribution(consensus, cfg),
  panel_02_method_agreement = render_panel_method_venn(consensus, cfg),
  panel_03_downstream_chain = render_panel_downstream_chain(consensus, cfg)
)
panel_tsv <- communication_report_panels_to_tsv(panels)
stopifnot(nrow(panel_tsv) == 3)
for (panel in panels) {
  for (path in unlist(panel$plots, use.names = FALSE)) {
    stopifnot(file.exists(path), file.info(path)$size > 0)
  }
}
report <- communication_report_assemble(panels, consensus, filtered, pairs, cfg)
stopifnot(any(grepl("Panel 1", report, fixed = TRUE)))
stopifnot(any(grepl("Panel 2", report, fixed = TRUE)))
stopifnot(any(grepl("Panel 3", report, fixed = TRUE)))
stopifnot(!any(grepl("panel_04", names(panels), fixed = TRUE)))
cat("smoke_07e_eda_panels: ok\n")
RS
