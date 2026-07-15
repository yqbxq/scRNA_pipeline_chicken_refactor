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
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "cluster_marker_audit_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03c2_cluster_marker_audit"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03c <- read_manifest_local(cfg$module_03c_manifest_path)
candidate_rds <- resolve_output_local(manifest_03c, "candidate_clustered_object")
ranking_path <- resolve_output_local(manifest_03c, "resolution_ranking_tsv")

seu <- readRDS(candidate_rds)
seu <- maybe_join_layers(seu)
ranking <- read_tsv_optional(ranking_path)
if (nrow(ranking) == 0) {
  stop(sprintf("Missing resolution ranking rows: %s", ranking_path), call. = FALSE)
}
ranking <- ranking[order(as.integer(ranking$rank)), , drop = FALSE]
candidate_rows <- ranking[ranking$marker_audit_candidate %in% c(TRUE, "TRUE", "true", "yes", "1"), , drop = FALSE]
if (nrow(candidate_rows) == 0) {
  candidate_rows <- ranking[seq_len(min(nrow(ranking), cfg$cluster_marker_candidate_top_n)), , drop = FALSE]
}
candidate_rows <- candidate_rows[seq_len(min(nrow(candidate_rows), cfg$cluster_marker_candidate_top_n)), , drop = FALSE]

candidate_quality_rows <- list()
candidate_index_rows <- list()
for (i in seq_len(nrow(candidate_rows))) {
  row <- candidate_rows[i, , drop = FALSE]
  cluster_col <- row$candidate_cluster_col[[1]]
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(sprintf("Candidate cluster column missing: %s", cluster_col), call. = FALSE)
  }
  message("Marker audit candidate: ", row$candidate_id[[1]], " (", cluster_col, ")")
  audit <- run_cluster_marker_audit(
    seu,
    cluster_col = cluster_col,
    cfg = cfg,
    assay = cfg$cluster_marker_assay
  )
  out_dir <- file.path(cfg$cluster_marker_table_dir_layer, "candidates", row$candidate_id[[1]])
  paths <- write_cluster_marker_audit_outputs(audit, out_dir)
  quality <- candidate_marker_quality_row(row, audit, cfg, out_dir)
  candidate_quality_rows[[length(candidate_quality_rows) + 1]] <- quality
  candidate_index_rows[[length(candidate_index_rows) + 1]] <- data.frame(
    candidate_id = row$candidate_id[[1]],
    resolution = row$resolution[[1]],
    candidate_cluster_col = cluster_col,
    marker_dir = out_dir,
    markers_raw = paths$markers_raw,
    markers_scored = paths$markers_scored,
    markers_anno = paths$markers_anno,
    markers_strict = paths$markers_strict,
    cluster_marker_qc = paths$cluster_marker_qc,
    top20_markers = paths$top20_markers,
    top50_markers = paths$top50_markers,
    stringsAsFactors = FALSE
  )
}

candidate_quality <- dplyr::bind_rows(candidate_quality_rows)
candidate_quality$marker_quality_score <- marker_quality_score(candidate_quality)
candidate_quality <- candidate_quality %>%
  dplyr::arrange(dplyr::desc(marker_quality_score), dplyr::desc(marker_supported_cluster_fraction), resolution)
candidate_index <- dplyr::bind_rows(candidate_index_rows)

candidate_quality_tsv <- file.path(cfg$cluster_marker_table_dir_layer, "candidate_marker_quality.tsv")
candidate_index_tsv <- file.path(cfg$cluster_marker_table_dir_layer, "candidate_marker_audit_index.tsv")
report_path <- file.path(cfg$cluster_marker_report_dir_layer, "report.md")

write_tsv_local(candidate_quality, candidate_quality_tsv)
write_tsv_local(candidate_index, candidate_index_tsv)

report_lines <- c(
  "# 03c2 Cluster Marker Audit",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- assay: `%s`", cfg$cluster_marker_assay),
  sprintf("- candidate_count: `%s`", nrow(candidate_quality)),
  sprintf("- input_object: `%s`", candidate_rds),
  "",
  "## Candidate Marker Quality",
  render_markdown_table_local(candidate_quality[, intersect(c(
    "candidate_id", "resolution", "n_clusters", "marker_supported_cluster_fraction",
    "median_strict_marker_n", "minimum_strict_marker_n", "median_top20_pct_diff",
    "fraction_clusters_with_few_strict_markers", "marker_quality_score"
  ), colnames(candidate_quality)), drop = FALSE]),
  "",
  "## Outputs",
  sprintf("- candidate_marker_quality: `%s`", candidate_quality_tsv),
  sprintf("- candidate_marker_audit_index: `%s`", candidate_index_tsv),
  "",
  "This stage is panel-free. It does not assign cell types and does not write annotation metadata."
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03c2_manifest_path,
  new_outputs = list(
    candidate_marker_quality_tsv = build_output_entry(candidate_quality_tsv, "tsv", module_name, "one row per audited resolution candidate with marker-quality scores", base_dir = cfg$project_root, schema = infer_schema_from_df(candidate_quality)),
    candidate_marker_audit_index_tsv = build_output_entry(candidate_index_tsv, "tsv", module_name, "paths to per-candidate marker audit outputs", base_dir = cfg$project_root, schema = infer_schema_from_df(candidate_index)),
    report = build_output_entry(report_path, "md", module_name, "human-readable candidate marker audit report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03c_manifest = cfg$module_03c_manifest_path,
    candidate_clustered_object = candidate_rds,
    resolution_ranking_tsv = ranking_path
  ),
  version = cfg$module_version,
  depends_on = list(module_03c = cfg$module_03c_manifest_path)
)

message("03c2 complete. Candidate marker quality: ", candidate_quality_tsv)
