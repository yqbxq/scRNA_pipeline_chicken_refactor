#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "scdesign3_engine_helpers.R"))

load_required_packages(c("dplyr", "jsonlite"))

cfg <- get_single_script_config_04()
module_name <- "04e_scdesign3_engine"
prepare_dirs_04(cfg)
set.seed(cfg$random_seed)

table_dir <- file.path(cfg$table_dir, module_name)
figure_root <- file.path(cfg$figure_dir, module_name)
report_dir <- file.path(cfg$results_dir, "reports")
checkpoint_dir <- file.path(cfg$checkpoint_dir, module_name)
manifest_path <- file.path(cfg$manifest_dir, module_name, "_manifest.json")
ensure_dir(table_dir)
ensure_dir(figure_root)
ensure_dir(report_dir)
ensure_dir(checkpoint_dir)
ensure_dir(dirname(manifest_path))

targets_path <- env_or_default_03("SCDESIGN3_TARGETS_SHEET", file.path(cfg$metadata_dir, "scdesign3_targets.tsv"))
thresholds_path <- env_or_default_03("SCDESIGN3_THRESHOLDS_SHEET", file.path(cfg$metadata_dir, "scdesign3_thresholds.tsv"))
simulation_designs_path <- env_or_default_03("SCDESIGN3_SIMULATION_DESIGNS_SHEET", file.path(cfg$metadata_dir, "scdesign3_simulation_designs.tsv"))
preflight_target_gate_path <- file.path(cfg$table_dir, "04d_cluster_robustness", "target_gate_status.tsv")

targets <- read_tsv_optional(targets_path)
required_target_cols <- c(
  "target_id", "target_type", "layer_id", "input_object", "truth_col",
  "questions_covered", "n_simulations", "resolution_grid",
  "mixture_design", "primary_metric", "pass_threshold", "warn_threshold",
  "fail_threshold", "output_dir", "status"
)
missing_target_cols <- setdiff(required_target_cols, colnames(targets))
if (length(missing_target_cols) > 0) {
  stop(sprintf("scdesign3_targets.tsv missing columns: %s", paste(missing_target_cols, collapse = ", ")), call. = FALSE)
}

preflight <- read_tsv_optional(preflight_target_gate_path)
if (nrow(preflight) > 0) {
  for (col in c("target_id", "resolved_input_path", "gate_status", "gate_level", "reason")) {
    if (!col %in% colnames(preflight)) preflight[[col]] <- ""
  }
  preflight <- preflight[, c("target_id", "resolved_input_path", "gate_status", "gate_level", "reason"), drop = FALSE]
  colnames(preflight) <- c("target_id", "resolved_input_path", "preflight_gate_status", "preflight_gate_level", "preflight_reason")
  targets <- dplyr::left_join(targets, preflight, by = "target_id")
} else {
  targets$resolved_input_path <- ""
  targets$preflight_gate_status <- ""
  targets$preflight_gate_level <- ""
  targets$preflight_reason <- ""
}
for (col in c("resolved_input_path", "preflight_gate_status", "preflight_gate_level", "preflight_reason")) {
  if (!col %in% colnames(targets)) targets[[col]] <- ""
  targets[[col]][is.na(targets[[col]])] <- ""
}

targets_to_process <- targets[
  targets$status == "active" | targets$target_type == "cluster_robustness",
  ,
  drop = FALSE
]

engine_cfg <- list(
  seed = cfg$random_seed,
  n_simulations_default = env_integer_03("SCDESIGN3_N_SIM", 5L),
  n_cores = env_integer_03("SCDESIGN3_N_CORES", max(1L, env_integer_03("MAIN_THREADS", 4L))),
  family_use = env_or_default_03("SCDESIGN3_FAMILY", "nb"),
  max_cells_per_label = env_integer_03("SCDESIGN3_MAX_CELLS_PER_LABEL", 2000L),
  n_hvg = env_integer_03("SCDESIGN3_N_HVG", 2000L),
  n_pcs = env_integer_03("SCDESIGN3_N_PCS", 30L),
  resolution_default = scd_parse_resolution_grid(env_or_default_03("SCDESIGN3_RESOLUTION_GRID", "0.6")),
  checkpoint_dir = checkpoint_dir,
  figure_root = figure_root
)
force_rerun <- scd_bool(Sys.getenv("SCDESIGN3_ENGINE_FORCE", unset = "no"))

target_results <- list()
if (nrow(targets_to_process) > 0) {
  for (idx in seq_len(nrow(targets_to_process))) {
    target <- targets_to_process[idx, , drop = FALSE]
    target_id <- scd_scalar(target$target_id)
    checkpoint_path <- file.path(checkpoint_dir, paste0(target_id, ".rds"))
    if (!force_rerun && file.exists(checkpoint_path)) {
      result <- readRDS(checkpoint_path)
    } else {
      result <- run_cluster_target(target, engine_cfg)
      saveRDS(result, checkpoint_path)
    }
    target_results[[length(target_results) + 1L]] <- result
  }
}

target_metrics <- scd_bind_rows(lapply(target_results, `[[`, "target_metrics"), scd_target_metrics_cols)
per_simulation <- scd_bind_rows(lapply(target_results, `[[`, "per_simulation"), scd_per_simulation_cols)
per_label <- scd_bind_rows(lapply(target_results, `[[`, "per_label"), scd_per_label_cols)
engine_status <- scd_bind_rows(lapply(target_results, `[[`, "engine_status"), scd_engine_status_cols)

write_scdesign3_h5ad_mirror <- function(result, cfg, module_name) {
  synthetic <- result$synthetic_reference
  target_id <- if (nrow(result$engine_status) > 0) result$engine_status$target_id[[1]] else ""
  target_id <- scd_scalar(target_id, "unknown_target")
  mirror_dir <- file.path(cfg$results_dir, "90a_export_h5ad", "scdesign3_04e")
  ensure_dir(mirror_dir)
  target_dir <- file.path(table_dir, target_id)
  ensure_dir(target_dir)
  manifest_tsv <- file.path(target_dir, "synthetic_h5ad_manifest.tsv")
  h5ad_path <- file.path(mirror_dir, paste0(target_id, ".h5ad"))
  if (is.null(synthetic) || is.null(synthetic$counts) || is.null(synthetic$meta)) {
    row <- data.frame(
      artifact_id = target_id,
      artifact_role = "scdesign3_synthetic_reference",
      h5ad_path = h5ad_path,
      source_object = "scDesign3_04e",
      target_id = target_id,
      question_ids = "",
      n_obs = 0L,
      n_vars = 0L,
      fingerprint = "",
      status = "not_available",
      reason = "04e did not produce a synthetic reference count matrix for this target.",
      stringsAsFactors = FALSE
    )
    write_tsv_local(row, manifest_tsv)
    return(row)
  }
  counts <- synthetic$counts
  meta <- synthetic$meta
  counts_rds <- file.path(target_dir, "synthetic_counts.rds")
  meta_tsv <- file.path(target_dir, "synthetic_cell_metadata.tsv")
  genes_tsv <- file.path(target_dir, "synthetic_genes.tsv")
  mtx_path <- file.path(target_dir, "synthetic_counts.mtx")
  saveRDS(counts, counts_rds)
  write_tsv_local(meta, meta_tsv)
  genes <- data.frame(gene_id = rownames(counts), gene_symbol = rownames(counts), stringsAsFactors = FALSE)
  write_tsv_local(genes, genes_tsv)
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    row <- data.frame(
      artifact_id = target_id,
      artifact_role = "scdesign3_synthetic_reference",
      h5ad_path = h5ad_path,
      source_object = "scDesign3_04e",
      target_id = target_id,
      question_ids = "",
      n_obs = 0L,
      n_vars = 0L,
      fingerprint = "",
      status = "skipped_no_matrix_package",
      reason = "Matrix package is required to write Matrix Market inputs for H5AD mirror.",
      stringsAsFactors = FALSE
    )
    write_tsv_local(row, manifest_tsv)
    return(row)
  }
  Matrix::writeMM(Matrix::Matrix(counts, sparse = TRUE), mtx_path)
  writer <- file.path(cfg$pipeline_root, "workflow", "04python", "write_matrix_h5ad.py")
  run <- tryCatch(system2(Sys.getenv("PY_SPATIAL_BIN", unset = "python3"), args = c(
    writer,
    "--mtx", mtx_path,
    "--genes", genes_tsv,
    "--obs", meta_tsv,
    "--obs-id-col", "synthetic_cell_id",
    "--out-h5ad", h5ad_path,
    "--manifest", manifest_tsv,
    "--artifact-id", target_id,
    "--artifact-role", "scdesign3_synthetic_reference",
    "--source-object", "scDesign3_04e"
  ), stdout = TRUE, stderr = TRUE), error = function(e) structure(conditionMessage(e), status = 127))
  if (!file.exists(manifest_tsv)) {
    row <- data.frame(
      artifact_id = target_id,
      artifact_role = "scdesign3_synthetic_reference",
      h5ad_path = h5ad_path,
      source_object = "scDesign3_04e",
      target_id = target_id,
      question_ids = "",
      n_obs = 0L,
      n_vars = 0L,
      fingerprint = "",
      status = "failed_write_h5ad",
      reason = paste(as.character(run), collapse = " "),
      stringsAsFactors = FALSE
    )
    write_tsv_local(row, manifest_tsv)
  }
  row <- read_tsv_optional(manifest_tsv)
  for (col in c("target_id", "question_ids")) {
    if (!col %in% colnames(row)) row[[col]] <- if (identical(col, "target_id")) target_id else ""
  }
  row
}

scdesign3_h5ad_manifest <- scd_bind_rows(
  lapply(target_results, write_scdesign3_h5ad_mirror, cfg = cfg, module_name = module_name),
  c("artifact_id", "artifact_role", "h5ad_path", "source_object", "target_id", "question_ids", "n_obs", "n_vars", "fingerprint", "status", "reason")
)

paths <- list(
  target_metrics_tsv = file.path(table_dir, "target_metrics.tsv"),
  per_simulation_metrics_tsv = file.path(table_dir, "per_simulation_metrics.tsv"),
  per_label_metrics_tsv = file.path(table_dir, "per_label_metrics.tsv"),
  engine_status_tsv = file.path(table_dir, "engine_status.tsv"),
  scdesign3_h5ad_manifest_tsv = file.path(table_dir, "scdesign3_h5ad_manifest.tsv"),
  report_md = file.path(report_dir, "04e_scdesign3_engine.md"),
  figure_dir = figure_root
)

write_tsv_local(target_metrics, paths$target_metrics_tsv)
write_tsv_local(per_simulation, paths$per_simulation_metrics_tsv)
write_tsv_local(per_label, paths$per_label_metrics_tsv)
write_tsv_local(engine_status, paths$engine_status_tsv)
write_tsv_local(scdesign3_h5ad_manifest, paths$scdesign3_h5ad_manifest_tsv)

status_summary <- if (nrow(target_metrics) > 0) {
  target_metrics %>%
    dplyr::count(target_type, status, gate_status, name = "target_n") %>%
    dplyr::arrange(target_type, status, gate_status)
} else {
  scd_empty_df(c("target_type", "status", "gate_status", "target_n"))
}

scored_summary <- if (nrow(target_metrics) > 0) {
  target_metrics %>%
    dplyr::select(target_id, target_type, n_simulations_done, best_resolution, median_ARI, median_NMI, median_max_jaccard, primary_metric_value, gate_status, status, reason)
} else {
  scd_empty_df(c("target_id", "target_type", "n_simulations_done", "best_resolution", "median_ARI", "median_NMI", "median_max_jaccard", "primary_metric_value", "gate_status", "status", "reason"))
}

missing_packages <- scd_missing_engine_packages()
report_lines <- build_report_lines_v04(
  title = "04e scDesign3 Engine",
  header_bullets = c(
    "M1 executes cluster_robustness targets only; other target types remain registered for later milestones.",
    sprintf("targets_seen: `%s`; targets_processed: `%s`", nrow(targets), nrow(targets_to_process)),
    sprintf("n_simulations_default: `%s`; effective override: `%s`", engine_cfg$n_simulations_default, Sys.getenv("SCDESIGN3_ENGINE_N_SIM", unset = "")),
    sprintf("missing_engine_packages: `%s`", ifelse(length(missing_packages) == 0, "none", paste(missing_packages, collapse = ",")))
  ),
  key_files = list(
    scdesign3_targets = targets_path,
    scdesign3_thresholds = thresholds_path,
    scdesign3_simulation_designs = simulation_designs_path,
    target_metrics = paths$target_metrics_tsv,
    per_simulation_metrics = paths$per_simulation_metrics_tsv,
    per_label_metrics = paths$per_label_metrics_tsv,
    engine_status = paths$engine_status_tsv,
    scdesign3_h5ad_manifest = paths$scdesign3_h5ad_manifest_tsv,
    figures = paths$figure_dir
  ),
  review_focus = c(
    "Review target rows with status `ok`; these are the only rows eligible for 04f PASS/WARN/FAIL upgrade.",
    "Rows with `waiting_input` or `missing_engine_dependency` preserve the 04d preflight gate rather than becoming biological failures.",
    "Approve the `scdesign3_validated` gate only after ARI/NMI distributions and recovery plots are acceptable."
  ),
  extra_sections = list(
    "Target Status Summary" = render_markdown_table_local(status_summary),
    "Target Metrics" = render_markdown_table_local(scored_summary)
  )
)
write_markdown_local(report_lines, paths$report_md)

if (file.exists(manifest_path)) {
  unlink(manifest_path)
}
write_manifest_local(
  manifest_path = manifest_path,
  new_outputs = list(
    target_metrics_tsv = build_output_entry(paths$target_metrics_tsv, "tsv", module_name, "one row per scDesign3 target engine result", base_dir = cfg$project_root, schema = infer_schema_from_df(target_metrics)),
    per_simulation_metrics_tsv = build_output_entry(paths$per_simulation_metrics_tsv, "tsv", module_name, "one row per target, simulation, and resolution", base_dir = cfg$project_root, schema = infer_schema_from_df(per_simulation)),
    per_label_metrics_tsv = build_output_entry(paths$per_label_metrics_tsv, "tsv", module_name, "truth-label recovery metrics per simulation", base_dir = cfg$project_root, schema = infer_schema_from_df(per_label)),
    engine_status_tsv = build_output_entry(paths$engine_status_tsv, "tsv", module_name, "fit/simulate/score runtime status per target", base_dir = cfg$project_root, schema = infer_schema_from_df(engine_status)),
    scdesign3_h5ad_manifest_tsv = build_output_entry(paths$scdesign3_h5ad_manifest_tsv, "tsv", module_name, "H5AD mirror manifest for 04e synthetic references", base_dir = cfg$project_root, schema = infer_schema_from_df(scdesign3_h5ad_manifest)),
    figure_dir = build_output_entry(paths$figure_dir, "directory", module_name, "per-target scDesign3 diagnostic figures", base_dir = cfg$project_root),
    report = build_output_entry(paths$report_md, "md", module_name, "04e scDesign3 engine report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    scdesign3_targets = targets_path,
    scdesign3_thresholds = thresholds_path,
    scdesign3_simulation_designs = simulation_designs_path,
    preflight_target_gate_status = preflight_target_gate_path
  ),
  version = cfg$module_version,
  depends_on = list(
    metadata = list(targets = targets_path, thresholds = thresholds_path, simulation_designs = simulation_designs_path),
    module_04d = cfg$module_04d_manifest_path
  )
)

message("04e scDesign3 engine completed. targets processed: ", nrow(targets_to_process), "; status rows: ", nrow(engine_status))
