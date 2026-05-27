#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/manifest_utils.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/metadata_io.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/commot_signal_utils.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_08a_spatial_communication_io"
prepare_dirs_spatial(cfg)

enabled_flag <- tolower(cfg$commot_enabled)
sections <- spatial_read_tsv(cfg$section_sheet)
if (nrow(sections) == 0) {
  sections <- data.frame(section_id = "all", commot_enabled = enabled_flag, stringsAsFactors = FALSE)
}
for (col in c("section_id", "commot_enabled")) {
  if (!col %in% colnames(sections)) sections[[col]] <- ifelse(col == "commot_enabled", "auto", "")
}
sections$commot_enabled[!nzchar(sections$commot_enabled)] <- "auto"
sections <- sections[sections$commot_enabled != "no" & enabled_flag != "no", , drop = FALSE]

resolve_h5ad_for_section <- function(section_id) {
  explicit <- Sys.getenv("COMMOT_INPUT_H5AD", unset = "")
  if (nzchar(explicit) && file.exists(explicit)) return(explicit)
  candidates <- c(
    file.path(cfg$spatial_commot_input_dir, sprintf("%s_commot_input.h5ad", section_id)),
    file.path(cfg$spatial_commot_input_dir, sprintf("%s_spatial.h5ad", section_id))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) hit[[1]] else ""
}

input_rows <- if (nrow(sections) == 0) {
  data.frame(section_id = character(), h5ad_path = character(), label_col = character(), status = character(), reason = character(), stringsAsFactors = FALSE)
} else {
  do.call(rbind, lapply(seq_len(nrow(sections)), function(i) {
    section_id <- spatial_cell(sections[i, , drop = FALSE], "section_id", sprintf("section_%s", i))
    h5ad_path <- resolve_h5ad_for_section(section_id)
    status <- if (nzchar(h5ad_path)) "ready" else "skipped_missing_h5ad"
    reason <- if (nzchar(h5ad_path)) "" else "No COMMOT input H5AD was found; run 90_export_h5ad or provide COMMOT_INPUT_H5AD."
    data.frame(
      section_id = section_id,
      h5ad_path = h5ad_path,
      label_col = Sys.getenv("COMMOT_CLUSTER_COL", unset = "cell_type_main"),
      status = status,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }))
}

consensus <- read_tsv_optional(cfg$commot_lr_candidates_tsv)
lr_candidates <- prepare_commot_lr_candidates(consensus)

input_manifest_tsv <- file.path(cfg$spatial_commot_table_dir, "commot_input_manifest.tsv")
lr_candidates_tsv <- file.path(cfg$spatial_commot_table_dir, "commot_lr_candidates.tsv")
spatial_write_tsv(input_rows, input_manifest_tsv)
spatial_write_tsv(lr_candidates, lr_candidates_tsv)

if (file.exists(cfg$module_08a_spatial_communication_io_manifest_path)) {
  unlink(cfg$module_08a_spatial_communication_io_manifest_path)
}
st07_write_manifest_local(
  manifest_path = cfg$module_08a_spatial_communication_io_manifest_path,
  new_outputs = list(
    commot_input_manifest = build_output_entry(input_manifest_tsv, "tsv", module_name, "one row per section-level COMMOT H5AD input", base_dir = cfg$project_root, schema = infer_schema_from_df(input_rows)),
    commot_lr_candidates_tsv = build_output_entry(lr_candidates_tsv, "tsv", module_name, "one row per LR axis candidate passed to COMMOT", base_dir = cfg$project_root, schema = infer_schema_from_df(lr_candidates))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(section_sheet = cfg$section_sheet, commot_lr_candidates_source = cfg$commot_lr_candidates_tsv),
  version = cfg$module_07_version,
  depends_on = list(spatial_07e_deconvolution_compare = cfg$module_07e_deconvolution_compare_manifest_path)
)

message(sprintf("08a COMMOT IO completed. sections=%d lr_candidates=%d", nrow(input_rows), nrow(lr_candidates)))
