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
module_name <- "spatial_02b_finalize_clustering"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

backend_catalog <- c("b1_seurat_snn", "b2_bayesspace", "b3_spagcn", "b4_stagate")

parse_cli <- function(args) {
  out <- list(
    layer_file = cfg$spatial_object_layer_file,
    backends = "",
    override_file = cfg$selected_clustering_backend_override_file,
    marker_panel_dir = file.path(cfg$config_dir, "marker_panels")
  )
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--layer-file") && i < length(args)) {
      out$layer_file <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--backends") && i < length(args)) {
      out$backends <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--override-file") && i < length(args)) {
      out$override_file <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--marker-panel-dir") && i < length(args)) {
      out$marker_panel_dir <- args[[i + 1L]]
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

resolve_backends <- function(value, layer_settings) {
  raw <- if (nzchar(value)) value else spatial_cell(as.data.frame(layer_settings, stringsAsFactors = FALSE), "clustering_backends", cfg$clustering_backends)
  tokens <- spatial_tokenize(raw)
  if (length(tokens) == 0 || any(tolower(tokens) == "all")) {
    tokens <- backend_catalog
  }
  tokens <- unique(tokens)
  unknown <- setdiff(tokens, backend_catalog)
  if (length(unknown) > 0) {
    stop(sprintf("unknown spatial clustering backend(s): %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  unique(c("b1_seurat_snn", tokens))
}

write_variant <- function(obj, backend) {
  out_dir <- file.path(cfg$clustering_variants_dir, backend)
  ensure_dir(out_dir)
  out_path <- file.path(out_dir, sprintf("spatial_panorama_clustered_%s.rds", backend))
  saveRDS(obj, out_path)
  out_path
}

run_backend <- function(backend, panorama, reduction, dims, target_k, res_range, res_fine_step) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    if (identical(backend, "b1_seurat_snn")) {
      obj <- cluster_b1_seurat(panorama, reduction = reduction, dims = dims, target_k = target_k, res_range = res_range, res_fine_step = res_fine_step)
      list(status = "ok", obj = obj, cluster_col = "cluster_b1_seurat_snn", message = "")
    } else if (identical(backend, "b2_bayesspace")) {
      obj <- cluster_b2_bayesspace(panorama, target_k = target_k, dims = dims)
      list(status = "ok", obj = obj, cluster_col = "cluster_b2_bayesspace", message = "")
    } else if (identical(backend, "b3_spagcn")) {
      bridge <- cluster_b3_spagcn_bridge(panorama, cfg$spagcn_py_bin, cfg$spagcn_py_script, target_k = target_k, work_dir = file.path(cfg$clustering_variants_dir, backend, "bridge"), reduction = reduction)
      obj <- panorama
      obj$cluster_b3_spagcn <- bridge$clusters
      obj@misc$b3_spagcn_bridge <- bridge
      list(status = "ok", obj = obj, cluster_col = "cluster_b3_spagcn", message = "")
    } else if (identical(backend, "b4_stagate")) {
      bridge <- cluster_b4_stagate_bridge(panorama, cfg$stagate_py_bin, cfg$stagate_py_script, target_k = target_k, work_dir = file.path(cfg$clustering_variants_dir, backend, "bridge"), reduction = reduction)
      obj <- panorama
      obj$cluster_b4_stagate <- bridge$clusters
      obj@misc$b4_stagate_bridge <- bridge
      list(status = "ok", obj = obj, cluster_col = "cluster_b4_stagate", message = "")
    } else {
      stop(sprintf("unsupported backend: %s", backend), call. = FALSE)
    }
  }, error = function(e) {
    bayesspace_dep_missing <- identical(backend, "b2_bayesspace") && (
      !requireNamespace("BayesSpace", quietly = TRUE) ||
        !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE)
    )
    status <- if (backend %in% c("b3_spagcn", "b4_stagate")) "failed_py_bridge" else if (bayesspace_dep_missing) "failed_dependency" else "failed"
    list(status = status, obj = NULL, cluster_col = "", message = conditionMessage(e))
  })
  result$timing_sec <- round(proc.time()[["elapsed"]] - started, 3)
  result$reduction <- reduction
  result
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!file.exists(cfg$spatial_panorama_integrated_rds)) {
  stop(sprintf("missing integrated panorama: %s", cfg$spatial_panorama_integrated_rds), call. = FALSE)
}
panorama <- readRDS(cfg$spatial_panorama_integrated_rds)
layer_settings <- spatial_layer_settings(args$layer_file, "panorama_st")
backends <- resolve_backends(args$backends, layer_settings)
target_k <- as.integer(spatial_numeric_or(layer_settings$target_clusters, cfg$clustering_target))
dims <- spatial_parse_dims(layer_settings$pca_dims %||% "1:30", fallback = 1:30)
res_range <- spatial_parse_numeric_vector(layer_settings$res_range %||% "", fallback = c(0.4, 0.6, 0.8))
res_fine_step <- spatial_numeric_or(layer_settings$res_fine_step, 0.05)
integration_mode <- select_integration_mode(args$layer_file, default = cfg$default_integration_mode)
reduction <- switch(integration_mode, harmony = "harmony", cca = "cca_integrated", none = "pca_none", integration_mode)
if (!reduction %in% names(panorama@reductions)) {
  reduction <- "pca_none"
}
default_backend <- spatial_cell(as.data.frame(layer_settings, stringsAsFactors = FALSE), "selected_clustering_backend", cfg$default_clustering_backend)
if (!nzchar(default_backend)) {
  default_backend <- cfg$default_clustering_backend
}

panel <- tryCatch(read_spatial_region_panel(args$marker_panel_dir, layer_id = "panorama_st"), error = function(e) data.frame(layer_id = character(), celltype = character(), gene = character(), evidence_source = character()))

results <- list()
variant_paths <- list()
for (backend in backend_catalog) {
  if (!backend %in% backends) {
    results[[backend]] <- list(status = "skipped", obj = NULL, cluster_col = "", message = "not requested", timing_sec = NA_real_, reduction = reduction)
    next
  }
  result <- run_backend(backend, panorama, reduction, dims, target_k, res_range, res_fine_step)
  if (identical(backend, "b1_seurat_snn") && !identical(result$status, "ok")) {
    stop(sprintf("required b1_seurat_snn backend failed: %s", result$message), call. = FALSE)
  }
  if (identical(result$status, "ok")) {
    result$output_rds <- write_variant(result$obj, backend)
    variant_paths[[backend]] <- result$output_rds
  }
  results[[backend]] <- result
}

selected_backend <- select_clustering_backend(args$override_file, results, default = default_backend)
if (!nzchar(selected_backend)) {
  stop("no successful spatial clustering backend is available", call. = FALSE)
}
selected_obj <- results[[selected_backend]]$obj
selected_col <- results[[selected_backend]]$cluster_col
selected_obj$cluster_default <- factor(selected_obj@meta.data[[selected_col]])
selected_obj@misc$clustering <- list(
  selected_backend = selected_backend,
  selected_cluster_col = selected_col,
  integration_mode = integration_mode,
  reduction = reduction,
  target_clusters = target_k,
  decided_by = if (nzchar(read_clustering_override(args$override_file))) "override_file" else "pipeline_default"
)
selected_variant_path <- write_variant(selected_obj, selected_backend)
results[[selected_backend]]$obj <- selected_obj
results[[selected_backend]]$output_rds <- selected_variant_path
if (!file.copy(selected_variant_path, cfg$spatial_panorama_clustered_rds, overwrite = TRUE)) {
  stop(sprintf("failed to copy selected clustered panorama to %s", cfg$spatial_panorama_clustered_rds), call. = FALSE)
}

cluster_metrics <- summarize_clustering(results, panel)
selected_df <- data.frame(
  default_choice = default_backend,
  override_choice = read_clustering_override(args$override_file),
  final_choice = selected_backend,
  cluster_col = selected_col,
  integration_mode = integration_mode,
  reduction = reduction,
  target_clusters = target_k,
  canonical_rds = relative_path_local(cfg$spatial_panorama_clustered_rds, cfg$project_root),
  selected_variant_rds = relative_path_local(selected_variant_path, cfg$project_root),
  decision_timestamp = as.character(Sys.time()),
  decided_by = if (nzchar(read_clustering_override(args$override_file))) "override_file" else "pipeline_default",
  stringsAsFactors = FALSE
)

spatial_write_tsv(cluster_metrics, file.path(cfg$clustering_compare_dir, "cluster_metrics.tsv"))
spatial_write_tsv(selected_df, cfg$selected_backend_tsv)

report_lines <- c(
  "# Spatial Clustering Backend Compare",
  "",
  sprintf("- Integration mode: `%s`", integration_mode),
  sprintf("- Reduction: `%s`", reduction),
  sprintf("- Target clusters: `%d`", target_k),
  sprintf("- Selected backend: `%s`", selected_backend),
  "",
  "## Backend Metrics",
  render_markdown_table_local(cluster_metrics),
  "",
  "## Selected Backend",
  render_markdown_table_local(selected_df)
)
write_markdown_local(report_lines, file.path(cfg$clustering_compare_dir, "report.md"))

write_manifest_local(
  manifest_path = cfg$module_02b_clustering_manifest_path,
  new_outputs = list(
    panorama_clustered = build_output_entry(cfg$spatial_panorama_clustered_rds, "rds", module_name, "canonical clustered spatial panorama", base_dir = cfg$project_root),
    clustering_variants_dir = build_output_entry(cfg$clustering_variants_dir, "directory", module_name, "per-backend clustering variant RDS files", base_dir = cfg$project_root),
    cluster_metrics = build_output_entry(file.path(cfg$clustering_compare_dir, "cluster_metrics.tsv"), "tsv", module_name, "one row per clustering backend metrics and status", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_metrics)),
    selected_backend = build_output_entry(cfg$selected_backend_tsv, "tsv", module_name, "selected clustering backend decision", base_dir = cfg$project_root, schema = infer_schema_from_df(selected_df)),
    report = build_output_entry(file.path(cfg$clustering_compare_dir, "report.md"), "md", module_name, "spatial clustering backend comparison report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(panorama_integrated = cfg$spatial_panorama_integrated_rds, spatial_object_layer_file = args$layer_file, clustering_override_file = args$override_file),
  version = cfg$module_03_version,
  depends_on = list(spatial_02a_integration_eda = cfg$module_02a_integration_manifest_path)
)

message(sprintf("spatial clustering complete: selected %s", selected_backend))
