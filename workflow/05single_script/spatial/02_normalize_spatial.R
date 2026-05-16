#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")

load_required_packages(c("Seurat", "Matrix", "ggplot2", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_02_normalize"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

method_catalog <- c(
  "m0_no_normalization",
  "m1_lognormalize",
  "m2_sct_v1",
  "m3_sct_v2",
  "m4_pearson_residuals"
)
method_alias <- c(
  m0 = "m0_no_normalization",
  m1 = "m1_lognormalize",
  m2 = "m2_sct_v1",
  m3 = "m3_sct_v2",
  m4 = "m4_pearson_residuals",
  none = "m0_no_normalization",
  lognormalize = "m1_lognormalize",
  lognorm = "m1_lognormalize",
  sct = "m3_sct_v2",
  sct_v1 = "m2_sct_v1",
  sct_v2 = "m3_sct_v2",
  pearson = "m4_pearson_residuals",
  pearson_residuals = "m4_pearson_residuals"
)

parse_cli <- function(args) {
  out <- list(layer_file = cfg$spatial_object_layer_file, methods = "all", override_file = cfg$selected_method_override_file)
  i <- 1L
  positional <- character()
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--layer-file") && i < length(args)) {
      out$layer_file <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--methods") && i < length(args)) {
      out$methods <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--override-file") && i < length(args)) {
      out$override_file <- args[[i + 1L]]
      i <- i + 2L
    } else {
      positional <- c(positional, arg)
      i <- i + 1L
    }
  }
  if (length(positional) > 0 && file.exists(positional[[1]])) {
    out$layer_file <- positional[[1]]
  }
  out
}

parse_methods_arg <- function(value) {
  value <- trimws(as.character(value %||% "all"))
  if (!nzchar(value) || identical(tolower(value), "all")) {
    return(method_catalog)
  }
  tokens <- trimws(unlist(strsplit(value, "[,;[:space:]]+", perl = TRUE), use.names = FALSE))
  tokens <- tokens[nzchar(tokens)]
  resolved <- vapply(tokens, function(token) {
    key <- tolower(token)
    if (key %in% names(method_alias)) {
      unname(method_alias[[key]])
    } else {
      token
    }
  }, character(1))
  unknown <- setdiff(resolved, method_catalog)
  if (length(unknown) > 0) {
    stop(sprintf("unknown normalization method(s): %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  unique(resolved)
}

layer_defaults <- function(layer_file) {
  layers <- spatial_read_tsv(layer_file)
  defaults <- list(hvg_n = cfg$default_hvg_n, normalization_methods = "")
  if (nrow(layers) == 0) {
    return(defaults)
  }
  enabled_col <- if ("enabled" %in% colnames(layers)) tolower(layers$enabled) %in% c("yes", "true", "1", "on", "auto") else rep(TRUE, nrow(layers))
  layer <- layers[enabled_col, , drop = FALSE]
  if (nrow(layer) == 0) {
    layer <- layers[1, , drop = FALSE]
  } else {
    layer <- layer[1, , drop = FALSE]
  }
  hvg_col <- intersect(c("hvg_n", "hvg_nfeatures", "hvg_nfeatures_default"), colnames(layer))[1]
  if (!is.na(hvg_col)) {
    defaults$hvg_n <- as.integer(spatial_numeric_or(layer[[hvg_col]], cfg$default_hvg_n))
  }
  method_col <- intersect(c("normalization_method", "normalization_methods"), colnames(layer))[1]
  if (!is.na(method_col)) {
    defaults$normalization_methods <- spatial_cell(layer, method_col, "")
  }
  defaults
}

run_one_method <- function(method, obj, hvg_n) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    obj_m <- switch(
      method,
      m0_no_normalization = run_normalize_m0(obj, hvg_n = hvg_n),
      m1_lognormalize = run_normalize_m1(obj, hvg_n = hvg_n),
      m2_sct_v1 = run_normalize_m2(obj, hvg_n = hvg_n),
      m3_sct_v2 = run_normalize_m3(obj, hvg_n = hvg_n),
      m4_pearson_residuals = run_normalize_m4_py_bridge(obj, py_bin = cfg$py_spatial_bin, script_path = cfg$pearson_residuals_py, hvg_n = hvg_n)
    )
    obj_m@misc$normalization <- list(
      method = method,
      hvg = Seurat::VariableFeatures(obj_m),
      hvg_n = length(Seurat::VariableFeatures(obj_m)),
      timestamp = as.character(Sys.time())
    )
    list(status = "ok", obj = obj_m, message = "")
  }, error = function(e) {
    status <- if (identical(method, "m4_pearson_residuals")) "failed_py_bridge" else "failed"
    list(status = status, obj = NULL, message = conditionMessage(e))
  })
  result$timing_sec <- round(proc.time()[["elapsed"]] - started, 3)
  result
}

write_section_tables <- function(section_id, hvg_iou, pca_corr, coherence, timing) {
  if (nrow(hvg_iou) > 0) hvg_iou$section_id <- section_id
  if (nrow(pca_corr) > 0) pca_corr$section_id <- section_id
  if (nrow(coherence) > 0) coherence$section_id <- section_id
  if (nrow(timing) > 0) timing$section_id <- section_id
  list(hvg_iou = hvg_iou, pca_corr = pca_corr, coherence = coherence, timing = timing)
}

bind_tables <- function(tables) {
  tables <- tables[vapply(tables, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(tables) == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  all_cols <- unique(unlist(lapply(tables, colnames), use.names = FALSE))
  tables <- lapply(tables, function(df) {
    missing <- setdiff(all_cols, colnames(df))
    for (col in missing) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, tables)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
layer_info <- layer_defaults(args$layer_file)
methods_to_run <- parse_methods_arg(args$methods)
if (!identical(tolower(args$methods), "all") && length(methods_to_run) == 0 && nzchar(layer_info$normalization_methods)) {
  methods_to_run <- parse_methods_arg(layer_info$normalization_methods)
}
hvg_n <- layer_info$hvg_n

post_qc_paths <- list.files(cfg$spatial_checkpoint_dir, pattern = "_post_qc\\.rds$", full.names = TRUE)
if (length(post_qc_paths) == 0) {
  stop(sprintf("no post-QC spatial objects found in %s", cfg$spatial_checkpoint_dir), call. = FALSE)
}

all_hvg <- list()
all_pca <- list()
all_coherence <- list()
all_timing <- list()
selected_rows <- list()
canonical_paths <- character()

for (post_qc_path in post_qc_paths) {
  obj <- readRDS(post_qc_path)
  ids <- spatial_object_ids(obj, fallback = sub("_post_qc\\.rds$", "", basename(post_qc_path)))
  stem <- sub("_post_qc\\.rds$", "", basename(post_qc_path))
  section_dir <- file.path(cfg$normalization_variants_dir, stem)
  ensure_dir(section_dir)

  results <- list()
  for (method in method_catalog) {
    if (!method %in% methods_to_run) {
      results[[method]] <- list(status = "skipped", obj = NULL, message = "not requested", timing_sec = NA_real_)
      next
    }
    result <- run_one_method(method, obj, hvg_n = hvg_n)
    if (!is.null(result$obj)) {
      out_path <- file.path(section_dir, sprintf("%s_post_norm.rds", method))
      saveRDS(result$obj, out_path)
      result$output_rds <- out_path
    }
    results[[method]] <- result
  }

  successful <- names(results)[vapply(results, function(x) identical(x$status, "ok") && !is.null(x$obj), logical(1))]
  if (length(successful) == 0) {
    stop(sprintf("all normalization methods failed for section_id=%s sample_id=%s", ids$section_id, ids$sample_id), call. = FALSE)
  }

  override <- read_normalization_override(args$override_file, section_id = ids$section_id)
  final_choice <- select_normalization_method(cfg$default_normalization_method, override, available_methods = successful)
  selected_path <- results[[final_choice]]$output_rds
  canonical_path <- file.path(cfg$spatial_checkpoint_dir, sprintf("%s_post_norm.rds", stem))
  if (!file.copy(selected_path, canonical_path, overwrite = TRUE)) {
    stop(sprintf("failed to copy selected normalization object for section_id=%s", ids$section_id), call. = FALSE)
  }
  canonical_paths <- c(canonical_paths, canonical_path)

  hvg_iou <- compute_hvg_iou_matrix(results)
  pca_corr <- compute_pca_confounder_correlation(results, c("nCount_Spatial", "percent.mito"))
  coherence <- compute_marker_spatial_coherence(results, marker_panel = NULL)
  timing <- summarize_timing(results)
  table_set <- write_section_tables(ids$section_id, hvg_iou, pca_corr, coherence, timing)
  all_hvg[[length(all_hvg) + 1L]] <- table_set$hvg_iou
  all_pca[[length(all_pca) + 1L]] <- table_set$pca_corr
  all_coherence[[length(all_coherence) + 1L]] <- table_set$coherence
  all_timing[[length(all_timing) + 1L]] <- table_set$timing

  selected_rows[[length(selected_rows) + 1L]] <- data.frame(
    sample_id = ids$sample_id,
    section_id = ids$section_id,
    default_choice = cfg$default_normalization_method,
    override_choice = override,
    final_choice = final_choice,
    decision_timestamp = as.character(Sys.time()),
    decided_by = if (nzchar(override)) "override_file" else "pipeline_default",
    canonical_rds = relative_path_local(canonical_path, cfg$project_root),
    stringsAsFactors = FALSE
  )
}

hvg_df <- bind_tables(all_hvg)
pca_df <- bind_tables(all_pca)
coherence_df <- bind_tables(all_coherence)
timing_df <- bind_tables(all_timing)
selected_df <- do.call(rbind, selected_rows)

hvg_path <- file.path(cfg$normalization_compare_dir, "hvg_iou_matrix.tsv")
pca_path <- file.path(cfg$normalization_compare_dir, "pca_confounder_corr.tsv")
coherence_path <- file.path(cfg$normalization_compare_dir, "marker_spatial_coherence.tsv")
timing_path <- file.path(cfg$normalization_compare_dir, "method_timing.tsv")
spatial_write_tsv(hvg_df, hvg_path)
spatial_write_tsv(pca_df, pca_path)
spatial_write_tsv(coherence_df, coherence_path)
spatial_write_tsv(timing_df, timing_path)
spatial_write_tsv(selected_df, cfg$selected_method_tsv)

report_lines <- c(
  "# Spatial Normalization Compare Report",
  "",
  sprintf("- Methods requested: `%s`", paste(methods_to_run, collapse = ", ")),
  sprintf("- Default method: `%s`", cfg$default_normalization_method),
  sprintf("- Override file: `%s`", relative_path_local(args$override_file, cfg$project_root)),
  "",
  "## Selected Method",
  render_markdown_table_local(selected_df),
  "",
  "## Method Timing",
  render_markdown_table_local(timing_df),
  "",
  "## HVG IoU Matrix",
  render_markdown_table_local(hvg_df),
  "",
  "## PCA Confounder Correlation",
  render_markdown_table_local(pca_df),
  "",
  "## Spatial Coherence",
  render_markdown_table_local(coherence_df),
  "",
  "## Review Note",
  "Normalization uses an auditable default. During the downstream spatial integration gate, edit the override TSV and rerun this stage if another method is preferred."
)
write_markdown_local(report_lines, file.path(cfg$normalization_compare_dir, "report.md"))

write_manifest_local(
  manifest_path = cfg$module_02_norm_manifest_path,
  new_outputs = list(
    normalization_variant_dir = build_output_entry(cfg$normalization_variants_dir, "directory", module_name, "per-section normalization variants", base_dir = cfg$project_root),
    normalized_object_dir = build_output_entry(cfg$spatial_checkpoint_dir, "directory", module_name, "selected per-section normalized Seurat objects", base_dir = cfg$project_root),
    report = build_output_entry(file.path(cfg$normalization_compare_dir, "report.md"), "md", module_name, "normalization comparison report", base_dir = cfg$project_root),
    selected_method = build_output_entry(cfg$selected_method_tsv, "tsv", module_name, "selected normalization method per section", base_dir = cfg$project_root, schema = infer_schema_from_df(selected_df)),
    hvg_iou_matrix = build_output_entry(hvg_path, "tsv", module_name, "HVG IoU matrices by section", base_dir = cfg$project_root, schema = infer_schema_from_df(hvg_df)),
    pca_confounder_corr = build_output_entry(pca_path, "tsv", module_name, "PCA confounder correlations", base_dir = cfg$project_root, schema = infer_schema_from_df(pca_df)),
    marker_spatial_coherence = build_output_entry(coherence_path, "tsv", module_name, "spatial coherence fallback scores", base_dir = cfg$project_root, schema = infer_schema_from_df(coherence_df)),
    method_timing = build_output_entry(timing_path, "tsv", module_name, "normalization method runtime and status", base_dir = cfg$project_root, schema = infer_schema_from_df(timing_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    post_qc_object_dir = cfg$spatial_checkpoint_dir,
    spatial_object_layer_file = args$layer_file,
    selected_method_override_file = args$override_file
  ),
  version = cfg$module_02_version,
  depends_on = list(
    spatial_01b_qc_filter = cfg$module_02_filter_manifest_path,
    spatial_01c_post_qc_eda = cfg$module_02_eda_manifest_path
  )
)

message(sprintf("spatial normalization complete: %s sections", nrow(selected_df)))
