#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_pairs_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_report_panels.R"))
source_utf8(file.path(.script_dir, "helpers", "scTenifoldKnk_hook.R"))

cfg <- get_single_script_config_07()
module_name <- "07e_communication_eda"
prepare_dirs_07(cfg)

consensus_df <- read_tsv_optional(cfg$communication_method_consensus_tsv)
if (nrow(consensus_df) == 0) {
  evidence_df <- read_tsv_optional(cfg$communication_evidence_tier_tsv)
  consensus_df <- evidence_df
}
consensus_df <- communication_report_normalize_consensus(consensus_df)

pairs <- read_communication_pairs(cfg)
filtered_df <- communication_report_filter_by_tier(
  consensus_df,
  pairs,
  filter_tier = cfg$communication_report_filter_tier
)

write_tsv_local(filtered_df, cfg$communication_eda_filtered_tsv)

panels <- list(
  panel_01_tier_distribution = render_panel_tier_distribution(consensus_df, cfg),
  panel_02_method_agreement = render_panel_method_venn(consensus_df, cfg),
  panel_03_downstream_chain = render_panel_downstream_chain(consensus_df, cfg),
  panel_04_scTenifoldKnk_hook = render_panel_scTenifoldKnk_hook(consensus_df, cfg),
  panel_legacy_scdesign3 = render_panel_legacy_scdesign3(cfg),
  panel_legacy_fallback = render_panel_legacy_fallback(cfg),
  panel_legacy_missing_gene_program = render_panel_legacy_missing_gene_program(cfg)
)

panels_tsv <- communication_report_panels_to_tsv(panels)
write_tsv_local(panels_tsv, cfg$communication_eda_panels_tsv)

report_lines <- communication_report_assemble(
  panels = panels,
  consensus_df = consensus_df,
  filtered_df = filtered_df,
  pairs = pairs,
  cfg = cfg
)

writeLines(report_lines, cfg$communication_report_md, useBytes = TRUE)

html_enabled <- tolower(cfg$communication_report_html_enabled) %in% c("yes", "true", "1", "on")
if (isTRUE(html_enabled)) {
  communication_report_write_basic_html(cfg$communication_report_md, cfg$communication_report_html)
} else if (file.exists(cfg$communication_report_html)) {
  unlink(cfg$communication_report_html)
}

outputs <- list(
  report_md = build_output_entry(
    cfg$communication_report_md,
    "md",
    module_name,
    "human-readable 07e communication EDA report",
    base_dir = cfg$project_root
  ),
  communication_eda_panels_tsv = build_output_entry(
    cfg$communication_eda_panels_tsv,
    "tsv",
    module_name,
    "one row per rendered communication EDA panel",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(panels_tsv)
  ),
  communication_eda_filtered_tsv = build_output_entry(
    cfg$communication_eda_filtered_tsv,
    "tsv",
    module_name,
    "communication consensus axes after report-tier filtering",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(filtered_df)
  )
)

if (isTRUE(html_enabled)) {
  outputs$report_html <- build_output_entry(
    cfg$communication_report_html,
    "html",
    module_name,
    "simple HTML companion for the 07e communication EDA report",
    base_dir = cfg$project_root
  )
}

if (file.exists(cfg$module_07e_manifest_path)) {
  unlink(cfg$module_07e_manifest_path)
}
communication_report_write_manifest(
  manifest_path = cfg$module_07e_manifest_path,
  outputs = outputs,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    method_consensus_tsv = cfg$communication_method_consensus_tsv,
    evidence_tier_tsv = cfg$communication_evidence_tier_tsv,
    communication_pairs_tsv = cfg$communication_pairs_sheet,
    communication_report_filter_tier = cfg$communication_report_filter_tier,
    communication_report_top_n_primary = cfg$communication_report_top_n_primary
  ),
  version = cfg$module_07e_eda_version,
  depends_on = list(
    module_07d = cfg$module_07d_manifest_path
  )
)

message(sprintf(
  "07e completed. axes=%d filtered=%d report=%s",
  nrow(consensus_df),
  nrow(filtered_df),
  cfg$communication_report_md
))
