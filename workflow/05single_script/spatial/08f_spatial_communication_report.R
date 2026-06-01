#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/manifest_utils.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/report_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_communication_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08f_spatial_communication_report"
prepare_dirs_spatial(cfg)

consensus <- st07_read_tsv(file.path(cfg$spatial_communication_table_dir, "spatial_communication_consensus.tsv"))
if (nrow(consensus) == 0) consensus <- st08_empty_consensus()

tier_counts <- if (nrow(consensus) == 0) data.frame(spatial_evidence_tier = "empty", n = 0L, stringsAsFactors = FALSE) else as.data.frame(table(spatial_evidence_tier = consensus$spatial_evidence_tier), stringsAsFactors = FALSE)
level_counts <- if (nrow(consensus) == 0) data.frame(final_interpretation_level = "empty", n = 0L, stringsAsFactors = FALSE) else as.data.frame(table(final_interpretation_level = consensus$final_interpretation_level), stringsAsFactors = FALSE)

has_space <- any(consensus$final_interpretation_level %in% c("core", "exploratory"), na.rm = TRUE)
has_neighbor <- any(tolower(consensus$spatial_evidence_tier) %in% c("spatial_primary", "spatial_supporting"), na.rm = TRUE)
has_final <- any(consensus$final_interpretation_level == "core", na.rm = TRUE) || any(consensus$final_interpretation_level == "exploratory", na.rm = TRUE)

question_gates <- do.call(rbind, list(
  st08_gate_row("I19_comm_in_space", module_name, if (has_space) "WARN" else "FAIL", if (has_space) "yes" else "no", if (has_space) "Spatial communication candidates have ST support; review evidence tier before core interpretation." else "No communication candidate has sufficient ST spatial support.", if (has_space) "spatial_supporting" else "blocked"),
  st08_gate_row("I20_comm_neighbor_check", module_name, if (has_neighbor) "WARN" else "FAIL", if (has_neighbor) "yes" else "no", if (has_neighbor) "At least one candidate has neighborhood/spatial support." else "No candidate has neighborhood support from 06d.", if (has_neighbor) "spatial_supporting" else "blocked"),
  st08_gate_row("I21_comm_final_figure", module_name, if (has_final) "WARN" else "FAIL", if (has_final) "yes" else "no", if (has_final) "Final figure candidates are available but must be tier-filtered." else "No final spatial communication figure candidate is available.", if (has_final) "spatial_supporting" else "blocked")
))

panels <- consensus[consensus$final_interpretation_level %in% c("core", "exploratory"), c("comm_candidate_id", "lr_axis_id", "section_id", "spatial_evidence_tier", "final_interpretation_level"), drop = FALSE]
report_md <- file.path(cfg$spatial_communication_report_dir, "spatial_communication_report.md")
gate_tsv <- file.path(cfg$spatial_communication_table_dir, "spatial_communication_question_gate_status.tsv")
panels_tsv <- file.path(cfg$spatial_communication_table_dir, "final_spatial_communication_panels.tsv")

spatial_write_tsv(question_gates, gate_tsv)
spatial_write_tsv(panels, panels_tsv)
write_markdown_local(c(
  "# Spatial Communication Evidence Report",
  "",
  "COMMOT is treated as spatial support, not a standalone communication mechanism. CellChat-only candidates remain hypothesis-only unless supported by LIANA/MultiNicheNet/receiver DEG plus ST evidence.",
  "",
  "## Evidence Tier Counts",
  render_markdown_table_local(tier_counts),
  "",
  "## Interpretation Level Counts",
  render_markdown_table_local(level_counts),
  "",
  "## Question Gates",
  render_markdown_table_local(question_gates),
  "",
  "## Final Candidate Preview",
  render_markdown_table_local(utils::head(panels, 50))
), report_md)

st07_write_manifest_local(
  manifest_path = cfg$module_08f_spatial_communication_report_manifest_path,
  new_outputs = list(
    spatial_communication_report = build_output_entry(report_md, "md", module_name, "final ST 08 spatial communication evidence report", base_dir = cfg$project_root),
    spatial_communication_question_gate_status = build_output_entry(gate_tsv, "tsv", module_name, "I19/I20/I21 spatial communication question gates", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gates)),
    final_spatial_communication_panels = build_output_entry(panels_tsv, "tsv", module_name, "final figure candidate rows with sufficient spatial communication evidence", base_dir = cfg$project_root, schema = infer_schema_from_df(panels))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_communication_consensus = file.path(cfg$spatial_communication_table_dir, "spatial_communication_consensus.tsv")),
  version = cfg$module_07_version,
  depends_on = list(spatial_08e_spatial_communication_consensus = cfg$module_08e_spatial_communication_consensus_manifest_path)
)

message(sprintf("08f spatial communication report completed. gates=%d panels=%d", nrow(question_gates), nrow(panels)))
