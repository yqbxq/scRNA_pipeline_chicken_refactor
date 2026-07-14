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
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_validation_helpers.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03c_cluster"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

read_selected_value <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少 selected_integration.txt: %s", path), call. = FALSE)
  }
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) {
    stop(sprintf("selected_integration.txt 没有有效选择: %s", path), call. = FALSE)
  }
  lines[[1]]
}

panorama_spec <- panorama_layer_spec(cfg)
selected_value <- read_selected_value(cfg$selected_integration_file)
parts <- strsplit(selected_value, "__", fixed = TRUE)[[1]]
if (length(parts) != 2 || !all(nzchar(parts))) {
  stop(sprintf("selected integration 格式必须是 <normalization>__<integration>: %s", selected_value), call. = FALSE)
}
selected_norm <- parts[[1]]
selected_mode <- parts[[2]]
candidate_key <- paste0("candidate_", selected_norm, "_", selected_mode)

manifest_03a2 <- read_manifest_local(cfg$module_03a2_manifest_path)
candidate_rds <- resolve_output_local(manifest_03a2, candidate_key)
candidate_index <- read_tsv_optional(resolve_output_local(manifest_03a2, "candidate_index_tsv"))
candidate_row <- candidate_index[candidate_index$candidate_key == candidate_key, , drop = FALSE]
if (nrow(candidate_row) != 1) {
  stop(sprintf("候选索引中找不到选择: %s", candidate_key), call. = FALSE)
}

seu <- readRDS(candidate_rds)
seu <- maybe_join_layers(seu)
cluster_col <- paste0(panorama_spec$layer_id, "_cluster")
resolution_search_tsv <- file.path(cfg$integration_table_dir_layer, "resolution_search.tsv")
cluster_result <- run_resolution_search(
  seu,
  reduction_name = candidate_row$reduction_name[[1]],
  layer_spec = panorama_spec,
  seed = cfg$random_seed,
  log_path = resolution_search_tsv,
  decision_path = cfg$cluster_resolution_decision_file
)
selected_umap <- if ("umap_name" %in% colnames(candidate_row)) normalize_scalar_value(candidate_row$umap_name[[1]], "") else ""
seu <- finalize_layer_object(
  cluster_result$seu,
  reduction_name = candidate_row$reduction_name[[1]],
  umap_name = selected_umap,
  cluster_col = cluster_col
)
seu$cell_type <- ""
seu$cell_type_confidence <- "未定"
seu$annotation_relation <- "无关"
seu@misc$panorama_layer_spec <- panorama_spec
seu@misc$selected_integration <- selected_value
seu@misc$selected_resolution <- cluster_result$selected_resolution

ensure_dir(dirname(cfg$panorama_clustered_rds))
saveRDS(seu, cfg$panorama_clustered_rds)
ensure_dir(dirname(cfg$compat_clustered_rds))
if (!identical(normalizePath(cfg$panorama_clustered_rds, winslash = "/", mustWork = FALSE), normalizePath(cfg$compat_clustered_rds, winslash = "/", mustWork = FALSE))) {
  invisible(file.copy(cfg$panorama_clustered_rds, cfg$compat_clustered_rds, overwrite = TRUE))
}

cluster_summary <- build_cluster_summary_local(seu, cluster_col)
cluster_summary_csv <- file.path(cfg$integration_table_dir_layer, "cluster_summary.csv")
write_csv_local(cluster_summary, cluster_summary_csv)

selected_resolution_txt <- file.path(cfg$integration_table_dir_layer, "selected_resolution.txt")
writeLines(
  c(
    sprintf("selected_resolution=%s", cluster_result$selected_resolution),
    sprintf("selected_cluster_count=%s", cluster_result$selected_cluster_count),
    sprintf("target_clusters=%s", panorama_spec$target_clusters),
    sprintf("cluster_selection_mode=%s", cluster_result$selection_mode %||% panorama_spec$cluster_selection_mode),
    sprintf("fallback_used=%s", ifelse(cluster_result$fallback_used, "true", "false")),
    sprintf("selected_integration=%s", selected_value)
  ),
  selected_resolution_txt,
  useBytes = TRUE
)

layer_status_row <- data.frame(
  layer_id = panorama_spec$layer_id,
  status = "built",
  parent_layer = "",
  cluster_column = cluster_col,
  selected_normalization = selected_norm,
  selected_integration = selected_mode,
  selected_reduction = candidate_row$reduction_name[[1]],
  selected_umap = selected_umap,
  selected_resolution = as.character(cluster_result$selected_resolution),
  selected_cluster_count = as.character(cluster_result$selected_cluster_count),
  fallback_used = ifelse(cluster_result$fallback_used, "true", "false"),
  clustered_rds = normalizePath(cfg$panorama_clustered_rds, winslash = "/", mustWork = FALSE),
  annotated_rds = "",
  stringsAsFactors = FALSE
)
invisible(upsert_layer_status(cfg$layer_status_file, layer_status_row))

triage_df <- if (isTRUE(cluster_result$fallback_used)) {
  make_triage_row(
    sample_id = panorama_spec$layer_id,
    severity = "medium",
    signal_id = "resolution_fallback_used",
    evidence = sprintf("target_clusters=%s; selected_cluster_count=%s; selected_resolution=%s", panorama_spec$target_clusters, cluster_result$selected_cluster_count, cluster_result$selected_resolution)
  )
} else {
  empty_triage_df(include_sample = TRUE)
}
triage_tsv <- file.path(cfg$integration_table_dir_layer, "resolution_triage.tsv")
write_tsv_local(triage_df, triage_tsv)

cluster_selection_outputs <- list()
if (!is.null(cluster_result$output_paths)) {
  path_descriptions <- c(
    resolution_metrics_tsv = "per-resolution clustering validation metrics",
    adjacent_ari_tsv = "adjacent-resolution ARI stability",
    seed_stability_tsv = "same-resolution seed ARI stability",
    seed_ari_pairs_tsv = "pairwise seed ARI values",
    cluster_size_metrics_tsv = "cluster size and small-cluster metrics",
    silhouette_metrics_tsv = "silhouette metrics by resolution",
    ch_metrics_tsv = "Calinski-Harabasz metrics by resolution",
    stability_plateaus_tsv = "detected stable resolution plateaus",
    resolution_ranking_tsv = "ranked cluster resolution candidates",
    cluster_resolution_recommendation_tsv = "recommended cluster resolution",
    cluster_resolution_decision_tsv = "manual cluster resolution decision table"
  )
  for (entry_name in names(path_descriptions)) {
    out_path <- cluster_result$output_paths[[entry_name]] %||% ""
    if (nzchar(out_path) && file.exists(out_path)) {
      cluster_selection_outputs[[entry_name]] <- build_output_entry(out_path, "tsv", module_name, path_descriptions[[entry_name]], base_dir = cfg$project_root)
    }
  }
}

write_manifest_local(
  manifest_path = cfg$module_03c_manifest_path,
  new_outputs = c(list(
    clustered_object = build_output_entry(cfg$panorama_clustered_rds, "rds", module_name, "panorama clustered Seurat object", base_dir = cfg$project_root),
    compatibility_clustered_object = build_output_entry(cfg$compat_clustered_rds, "rds", module_name, "compatibility checkpoint for downstream legacy modules", base_dir = cfg$project_root),
    resolution_search_tsv = build_output_entry(resolution_search_tsv, "tsv", module_name, "resolution search records", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_result$search_table)),
    cluster_summary_csv = build_output_entry(cluster_summary_csv, "csv", module_name, "cluster composition summary by sample/group", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_summary)),
    selected_resolution_txt = build_output_entry(selected_resolution_txt, "txt", module_name, "selected resolution summary", base_dir = cfg$project_root),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root),
    resolution_triage_tsv = build_output_entry(triage_tsv, "tsv", module_name, "resolution triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
  ), cluster_selection_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03a2_manifest = cfg$module_03a2_manifest_path,
    selected_integration_file = cfg$selected_integration_file,
    selected_candidate = candidate_rds
  ),
  version = cfg$module_version,
  depends_on = list(module_03a2 = cfg$module_03a2_manifest_path, module_03b = cfg$module_03b_manifest_path)
)

message("03c 完成。clustered_object: ", cfg$panorama_clustered_rds)
