#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09g_trajectory_palantir"
prepare_dirs_09(cfg)

palantir_index_cols <- c(trajectory_method_index_cols_09, "fate_path")
units_all <- trajectory_execution_units_09(cfg)
enabled_mask <- if (nrow(units_all) > 0) {
  vapply(seq_len(nrow(units_all)), function(i) {
    trajectory_method_enabled_09(units_all[i, , drop = FALSE], "palantir")
  }, logical(1))
} else {
  rep(FALSE, nrow(units_all))
}
units <- units_all[enabled_mask, , drop = FALSE]
job_rows <- list()
index_rows <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  input_dir <- trajectory_python_input_dir_09(cfg, "palantir", pair_id, split_value)
  pseudo_csv <- trajectory_method_output_path_09(cfg, "palantir", pair_id, split_value, "palantir_pseudotime")
  fate_csv <- trajectory_method_output_path_09(cfg, "palantir", pair_id, split_value, "palantir_fate")
  entropy_csv <- trajectory_method_output_path_09(cfg, "palantir", pair_id, split_value, "palantir_entropy")
  fig_png <- trajectory_method_figure_path_09(cfg, "palantir", pair_id, split_value, "Figure_Palantir")
  status_path <- trajectory_method_output_path_09(cfg, "palantir", pair_id, split_value, "palantir_status", "txt")
  ensure_dir(dirname(pseudo_csv))
  ensure_dir(dirname(fig_png))

  result <- tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])
    trajectory_export_python_input_09(seu, input_dir, unit$coarse_label_var[[1]], root, unit$terminal_group[[1]])
    list(
      job = data.frame(
        pair_id = pair_id,
        split_value = split_value,
        input_rds = unit$input_rds[[1]],
        input_dir = input_dir,
        label_var = unit$coarse_label_var[[1]],
        n_cells = ncol(seu),
        pseudotime_csv = pseudo_csv,
        fate_csv = fate_csv,
        entropy_csv = entropy_csv,
        figure_png = fig_png,
        status = status_path,
        stringsAsFactors = FALSE
      ),
      index = NULL
    )
  }, error = function(e) {
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), pseudo_csv, row.names = FALSE)
    write.csv(data.frame(cell_id = character(0), terminal_state = character(0), probability = numeric(0)), fate_csv, row.names = FALSE)
    write.csv(data.frame(cell_id = character(0), entropy = numeric(0)), entropy_csv, row.names = FALSE)
    writeLines(paste("failed", conditionMessage(e), sep = "\t"), status_path, useBytes = TRUE)
    list(
      job = NULL,
      index = data.frame(
        pair_id = pair_id,
        split_value = split_value,
        method = "palantir",
        methods_enabled = "yes",
        input_rds = unit$input_rds[[1]],
        output_path = pseudo_csv,
        extra_path = entropy_csv,
        figure_path = fig_png,
        n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
        status = "failed",
        reason = conditionMessage(e),
        runtime_s = "",
        fate_path = fate_csv,
        stringsAsFactors = FALSE
      )
    )
  })
  if (!is.null(result$job)) {
    job_rows[[length(job_rows) + 1L]] <- result$job
  }
  if (!is.null(result$index)) {
    index_rows[[length(index_rows) + 1L]] <- result$index
  }
}

jobs <- if (length(job_rows) > 0) dplyr::bind_rows(job_rows) else trajectory_empty_df_09(c(
  "pair_id", "split_value", "input_rds", "input_dir", "label_var", "n_cells",
  "pseudotime_csv", "fate_csv", "entropy_csv", "figure_png", "status"
))
write_tsv_local(jobs, cfg$trajectory_palantir_jobs_tsv)

python_script <- normalizePath(file.path(.script_dir, "..", "04python", "palantir.py"), winslash = "/", mustWork = FALSE)
if (nrow(jobs) > 0) {
  ok <- trajectory_run_python_bridge_09(python_script, cfg$trajectory_palantir_jobs_tsv, cfg$trajectory_palantir_index_tsv)
  if (!ok) {
    index_fallback <- data.frame(
      pair_id = jobs$pair_id,
      split_value = jobs$split_value,
      method = "palantir",
      methods_enabled = "yes",
      input_rds = jobs$input_rds,
      output_path = jobs$pseudotime_csv,
      extra_path = jobs$entropy_csv,
      figure_path = jobs$figure_png,
      n_cells = jobs$n_cells,
      status = "skipped_python_unavailable",
      reason = "velocity python bridge could not be executed",
      runtime_s = "",
      fate_path = jobs$fate_csv,
      stringsAsFactors = FALSE
    )
    write_tsv_local(index_fallback, cfg$trajectory_palantir_index_tsv)
  }
} else {
  write_tsv_local(trajectory_empty_df_09(palantir_index_cols), cfg$trajectory_palantir_index_tsv)
}

if (nrow(units_all) > 0 && any(!enabled_mask)) {
  disabled <- units_all[!enabled_mask, , drop = FALSE]
  for (i in seq_len(nrow(disabled))) {
    unit <- disabled[i, , drop = FALSE]
    reason <- trajectory_method_disabled_reason_09(unit, "palantir")
    index_rows[[length(index_rows) + 1L]] <- data.frame(
      pair_id = unit$pair_id[[1]],
      split_value = display_scalar_value(unit$split_value[[1]], "pooled"),
      method = "palantir",
      methods_enabled = "no",
      input_rds = unit$input_rds[[1]],
      output_path = "",
      extra_path = "",
      figure_path = "",
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = "skipped_disabled",
      reason = reason,
      runtime_s = 0,
      fate_path = "",
      stringsAsFactors = FALSE
    )
  }
}

bridge_index <- read_tsv_optional(cfg$trajectory_palantir_index_tsv)
if (length(index_rows) > 0) {
  bridge_index <- dplyr::bind_rows(bridge_index, dplyr::bind_rows(index_rows))
}
if (nrow(bridge_index) == 0) {
  bridge_index <- trajectory_empty_df_09(palantir_index_cols)
}
write_tsv_local(bridge_index, cfg$trajectory_palantir_index_tsv)

dynamic_outputs <- list()
if (nrow(jobs) > 0) {
  for (i in seq_len(nrow(jobs))) {
    dynamic_outputs[[sprintf("palantir_%s_pseudotime", trajectory_unit_file_id_09(jobs$pair_id[[i]], jobs$split_value[[i]]))]] <- build_output_entry(jobs$pseudotime_csv[[i]], "csv", module_name, "Palantir pseudotime", base_dir = cfg$project_root)
    dynamic_outputs[[sprintf("palantir_%s_entropy", trajectory_unit_file_id_09(jobs$pair_id[[i]], jobs$split_value[[i]]))]] <- build_output_entry(jobs$entropy_csv[[i]], "csv", module_name, "Palantir entropy", base_dir = cfg$project_root)
    dynamic_outputs[[sprintf("palantir_%s_fate", trajectory_unit_file_id_09(jobs$pair_id[[i]], jobs$split_value[[i]]))]] <- build_output_entry(jobs$fate_csv[[i]], "csv", module_name, "Palantir fate probability", base_dir = cfg$project_root)
  }
}
outputs <- c(
  list(
    palantir_index_tsv = build_output_entry(cfg$trajectory_palantir_index_tsv, "tsv", module_name, "Palantir method status by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(bridge_index)),
    palantir_jobs_tsv = build_output_entry(cfg$trajectory_palantir_jobs_tsv, "tsv", module_name, "Palantir Python bridge job table", base_dir = cfg$project_root, schema = infer_schema_from_df(jobs))
  ),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09g_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path, trajectory_pairs = cfg$trajectory_pairs_sheet, python_script = python_script),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path)
)

message("09g completed. Palantir index: ", cfg$trajectory_palantir_index_tsv)
