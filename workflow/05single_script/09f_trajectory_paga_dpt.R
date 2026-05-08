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
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_09.R"))
source_utf8(file.path(.script_dir, "helpers", "trajectory_utils.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09f_trajectory_paga_dpt"
prepare_dirs_09(cfg)

paga_index_cols <- c(trajectory_method_index_cols_09, "h5ad_path")
units <- trajectory_execution_units_09(cfg)
job_rows <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  input_dir <- trajectory_python_input_dir_09(cfg, "paga_dpt", pair_id, split_value)
  out_csv <- trajectory_method_output_path_09(cfg, "paga_dpt", pair_id, split_value, "paga_pseudotime")
  conn_tsv <- trajectory_method_output_path_09(cfg, "paga_dpt", pair_id, split_value, "paga_connectivity", "tsv")
  h5ad_path <- trajectory_method_output_path_09(cfg, "paga_dpt", pair_id, split_value, "paga_dpt", "h5ad")
  fig_png <- trajectory_method_figure_path_09(cfg, "paga_dpt", pair_id, split_value, "Figure_PAGA")
  status_path <- trajectory_method_output_path_09(cfg, "paga_dpt", pair_id, split_value, "paga_status", "txt")
  ensure_dir(dirname(out_csv))
  ensure_dir(dirname(fig_png))

  tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])
    trajectory_export_python_input_09(seu, input_dir, unit$coarse_label_var[[1]], root, unit$terminal_group[[1]])
    job_rows[[length(job_rows) + 1L]] <<- data.frame(
      pair_id = pair_id,
      split_value = split_value,
      input_rds = unit$input_rds[[1]],
      input_dir = input_dir,
      label_var = unit$coarse_label_var[[1]],
      n_cells = ncol(seu),
      pseudotime_csv = out_csv,
      connectivity_tsv = conn_tsv,
      h5ad_path = h5ad_path,
      figure_png = fig_png,
      status = status_path,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0), label = character(0)), out_csv, row.names = FALSE)
    write_tsv_local(data.frame(source = character(0), target = character(0), connectivity = numeric(0)), conn_tsv)
    writeLines(paste("failed", conditionMessage(e), sep = "\t"), status_path, useBytes = TRUE)
    job_rows[[length(job_rows) + 1L]] <<- data.frame(
      pair_id = pair_id,
      split_value = split_value,
      input_rds = unit$input_rds[[1]],
      input_dir = input_dir,
      label_var = unit$coarse_label_var[[1]],
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      pseudotime_csv = out_csv,
      connectivity_tsv = conn_tsv,
      h5ad_path = h5ad_path,
      figure_png = fig_png,
      status = status_path,
      stringsAsFactors = FALSE
    )
  })
}

jobs <- if (length(job_rows) > 0) dplyr::bind_rows(job_rows) else trajectory_empty_df_09(c(
  "pair_id", "split_value", "input_rds", "input_dir", "label_var", "n_cells",
  "pseudotime_csv", "connectivity_tsv", "h5ad_path", "figure_png", "status"
))
write_tsv_local(jobs, cfg$trajectory_paga_jobs_tsv)

python_script <- normalizePath(file.path(.script_dir, "..", "04python", "paga_dpt.py"), winslash = "/", mustWork = FALSE)
if (nrow(jobs) > 0) {
  ok <- trajectory_run_python_bridge_09(python_script, cfg$trajectory_paga_jobs_tsv, cfg$trajectory_paga_index_tsv)
  if (!ok) {
    index_fallback <- data.frame(
      pair_id = jobs$pair_id,
      split_value = jobs$split_value,
      method = "paga_dpt",
      methods_enabled = "yes",
      input_rds = jobs$input_rds,
      output_path = jobs$pseudotime_csv,
      extra_path = jobs$connectivity_tsv,
      figure_path = jobs$figure_png,
      n_cells = jobs$n_cells,
      status = "skipped_python_unavailable",
      reason = "velocity python bridge could not be executed",
      runtime_s = "",
      h5ad_path = jobs$h5ad_path,
      stringsAsFactors = FALSE
    )
    write_tsv_local(index_fallback, cfg$trajectory_paga_index_tsv)
  }
} else {
  write_tsv_local(trajectory_empty_df_09(paga_index_cols), cfg$trajectory_paga_index_tsv)
}

index_df <- read_tsv_optional(cfg$trajectory_paga_index_tsv)
if (nrow(index_df) == 0) {
  index_df <- trajectory_empty_df_09(paga_index_cols)
}

dynamic_outputs <- list()
if (nrow(jobs) > 0) {
  for (i in seq_len(nrow(jobs))) {
    dynamic_outputs[[sprintf("paga_dpt_%s_pseudotime", trajectory_unit_file_id_09(jobs$pair_id[[i]], jobs$split_value[[i]]))]] <- build_output_entry(
      jobs$pseudotime_csv[[i]],
      "csv",
      module_name,
      "PAGA-DPT pseudotime",
      base_dir = cfg$project_root
    )
  }
}
outputs <- c(
  list(
    paga_dpt_index_tsv = build_output_entry(cfg$trajectory_paga_index_tsv, "tsv", module_name, "PAGA-DPT method status by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
    paga_dpt_jobs_tsv = build_output_entry(cfg$trajectory_paga_jobs_tsv, "tsv", module_name, "PAGA-DPT Python bridge job table", base_dir = cfg$project_root, schema = infer_schema_from_df(jobs))
  ),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09f_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path, python_script = python_script),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path)
)

message("09f completed. PAGA-DPT index: ", cfg$trajectory_paga_index_tsv)
