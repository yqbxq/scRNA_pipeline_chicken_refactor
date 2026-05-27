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
  communication_report_filter_tier = "primary",
  communication_report_top_n_primary = 1L,
  communication_fallback_summary_tsv = file.path(tmp, "fallback.tsv"),
  nichenet_index_tsv = file.path(tmp, "nichenet.tsv")
)
ensure_dir(dirname(cfg$communication_fallback_summary_tsv))
write_tsv_local(data.frame(pair_id = "P1", fallback_pair_id = "P0", stringsAsFactors = FALSE), cfg$communication_fallback_summary_tsv)
write_tsv_local(data.frame(pair_id = "P2", status = "missing_gene_program", deg_status = "missing_top_gene_tsv", stringsAsFactors = FALSE), cfg$nichenet_index_tsv)

consensus <- data.frame(
  lr_axis_id = c("L1|R1|A->B", "L2|R2|A->B"),
  pair_id = c("P1", "P2"),
  condition_value = c("syf", "f5"),
  cellchat_hit = c("yes", "yes"),
  liana_consensus_hit = c("yes", "no"),
  nichenet_hit = c("yes", "no"),
  commot_spatial_hit = c("no", "no"),
  evidence_tier = c("primary", "candidate"),
  can_be_primary = c("yes", "no"),
  stringsAsFactors = FALSE
)
pairs <- data.frame(pair_id = c("P1", "P2"), evidence_tier_required = c("primary", "candidate"), stringsAsFactors = FALSE)
filtered <- communication_report_filter_by_tier(consensus, pairs, "primary")
stopifnot(nrow(filtered) == 1, filtered$pair_id[[1]] == "P1")

fallback <- render_panel_legacy_fallback(cfg)
missing_gp <- render_panel_legacy_missing_gene_program(cfg)
stopifnot(identical(fallback$status, "available"))
stopifnot(identical(missing_gp$status, "available"))
cat("smoke_07e_legacy_compat: ok\n")
RS

Rscript -e "parse('${REPO_ROOT}/workflow/05single_script/07e_communication_eda.R'); parse('${REPO_ROOT}/workflow/05single_script/helpers/communication_report_utils.R'); parse('${REPO_ROOT}/workflow/05single_script/helpers/communication_report_panels.R')" >/dev/null
