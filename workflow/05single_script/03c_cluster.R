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
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03c_cluster"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

read_selected_value <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Missing selected_integration.txt: %s", path), call. = FALSE)
  }
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) {
    stop(sprintf("selected_integration.txt has no valid choice: %s", path), call. = FALSE)
  }
  lines[[1]]
}

panorama_spec <- panorama_layer_spec(cfg)
selected_value <- read_selected_value(cfg$selected_integration_file)
parts <- strsplit(selected_value, "__", fixed = TRUE)[[1]]
if (length(parts) != 2 || !all(nzchar(parts))) {
  stop(sprintf("selected integration must be <normalization>__<integration>: %s", selected_value), call. = FALSE)
}
selected_norm <- parts[[1]]
selected_mode <- parts[[2]]
candidate_key <- paste0("candidate_", selected_norm, "_", selected_mode)

manifest_03a2 <- read_manifest_local(cfg$module_03a2_manifest_path)
candidate_rds <- resolve_output_local(manifest_03a2, candidate_key)
candidate_index <- read_tsv_optional(resolve_output_local(manifest_03a2, "candidate_index_tsv"))
candidate_row <- candidate_index[candidate_index$candidate_key == candidate_key, , drop = FALSE]
if (nrow(candidate_row) != 1) {
  stop(sprintf("Selected candidate missing in index: %s", candidate_key), call. = FALSE)
}

seu <- readRDS(candidate_rds)
seu <- maybe_join_layers(seu)
cluster_result <- run_quality_resolution_search(
  seu,
  reduction_name = candidate_row$reduction_name[[1]],
  layer_spec = panorama_spec,
  seed = cfg$random_seed,
  max_metric_cells = cfg$cluster_metric_max_cells,
  seed_repeats = cfg$cluster_seed_repeats,
  subsample_repeats = cfg$cluster_subsample_repeats,
  subsample_fraction = cfg$cluster_subsample_fraction,
  candidate_top_n = cfg$cluster_candidate_top_n
)

seu <- cluster_result$seu
seu$cell_type <- ""
seu$cell_type_confidence <- "undetermined"
seu$annotation_relation <- "none"
seu@misc$panorama_layer_spec <- panorama_spec
seu@misc$selected_integration <- selected_value
seu@misc$provisional_selected_resolution <- cluster_result$selected_resolution
seu@misc$provisional_selected_cluster_col <- cluster_result$selected_cluster_col

ensure_dir(dirname(cfg$panorama_candidate_clustered_rds))
saveRDS(seu, cfg$panorama_candidate_clustered_rds)

resolution_metrics_tsv <- file.path(cfg$integration_table_dir_layer, "resolution_metrics.tsv")
resolution_ranking_tsv <- file.path(cfg$integration_table_dir_layer, "resolution_ranking.tsv")
adjacent_ari_tsv <- file.path(cfg$integration_table_dir_layer, "adjacent_ari.tsv")
cluster_assignments_tsv <- file.path(cfg$integration_table_dir_layer, "cluster_assignments.tsv")
selected_resolution_tsv <- file.path(cfg$integration_table_dir_layer, "selected_resolution_provisional.tsv")
cluster_summary_csv <- file.path(cfg$integration_table_dir_layer, "cluster_summary_provisional.csv")

write_tsv_local(cluster_result$resolution_metrics, resolution_metrics_tsv)
write_tsv_local(cluster_result$resolution_ranking, resolution_ranking_tsv)
write_tsv_local(cluster_result$adjacent_ari, adjacent_ari_tsv)
write_tsv_local(cluster_result$cluster_assignments, cluster_assignments_tsv)

selected_row <- cluster_result$resolution_ranking[cluster_result$resolution_ranking$rank == 1L, , drop = FALSE]
selected_row$selected_integration <- selected_value
selected_row$selection_stage <- "provisional_clustering_metrics_only"
write_tsv_local(selected_row, selected_resolution_tsv)

cluster_summary <- build_cluster_summary_local(seu, cluster_result$selected_cluster_col)
write_csv_local(cluster_summary, cluster_summary_csv)

layer_status_row <- data.frame(
  layer_id = panorama_spec$layer_id,
  status = "cluster_candidates",
  parent_layer = "",
  cluster_column = "",
  candidate_cluster_column = cluster_result$selected_cluster_col,
  selected_normalization = selected_norm,
  selected_integration = selected_mode,
  selected_reduction = candidate_row$reduction_name[[1]],
  selected_umap = candidate_row$umap_name[[1]],
  selected_resolution = as.character(cluster_result$selected_resolution),
  selected_cluster_count = as.character(cluster_result$selected_cluster_count),
  fallback_used = "false",
  clustered_rds = "",
  candidate_clustered_rds = normalizePath(cfg$panorama_candidate_clustered_rds, winslash = "/", mustWork = FALSE),
  annotated_rds = "",
  stringsAsFactors = FALSE
)
invisible(upsert_layer_status(cfg$layer_status_file, layer_status_row))

triage_rows <- list()
if (nrow(cluster_result$resolution_ranking) > 0) {
  weak <- cluster_result$resolution_ranking[cluster_result$resolution_ranking$rank == 1L, , drop = FALSE]
  if (nrow(weak) == 1 && is.finite(weak$clustering_score[[1]]) && weak$clustering_score[[1]] < 0.50) {
    triage_rows[[length(triage_rows) + 1]] <- make_triage_row(
      sample_id = panorama_spec$layer_id,
      severity = "medium",
      signal_id = "weak_resolution_selection_score",
      evidence = sprintf("resolution=%s; clustering_score=%.3f", weak$resolution[[1]], weak$clustering_score[[1]])
    )
  }
}
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
triage_tsv <- file.path(cfg$integration_table_dir_layer, "resolution_triage.tsv")
write_tsv_local(triage_df, triage_tsv)

write_manifest_local(
  manifest_path = cfg$module_03c_manifest_path,
  new_outputs = list(
    candidate_clustered_object = build_output_entry(cfg$panorama_candidate_clustered_rds, "rds", module_name, "Seurat object carrying candidate resolution cluster columns", base_dir = cfg$project_root),
    resolution_metrics_tsv = build_output_entry(resolution_metrics_tsv, "tsv", module_name, "one row per resolution with stability/separation/fragmentation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_result$resolution_metrics)),
    resolution_ranking_tsv = build_output_entry(resolution_ranking_tsv, "tsv", module_name, "resolution metrics sorted by clustering score", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_result$resolution_ranking)),
    adjacent_ari_tsv = build_output_entry(adjacent_ari_tsv, "tsv", module_name, "ARI between neighboring resolutions", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_result$adjacent_ari)),
    cluster_assignments_tsv = build_output_entry(cluster_assignments_tsv, "tsv", module_name, "cell-level candidate cluster assignments", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_result$cluster_assignments)),
    selected_resolution_tsv = build_output_entry(selected_resolution_tsv, "tsv", module_name, "provisional clustering-metric selected resolution", base_dir = cfg$project_root, schema = infer_schema_from_df(selected_row)),
    cluster_summary_csv = build_output_entry(cluster_summary_csv, "csv", module_name, "provisional cluster composition summary by sample/group", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_summary)),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root),
    resolution_triage_tsv = build_output_entry(triage_tsv, "tsv", module_name, "resolution triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
  ),
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

message("03c complete. Candidate clustered object: ", cfg$panorama_candidate_clustered_rds)
