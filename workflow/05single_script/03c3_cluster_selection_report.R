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
source_utf8(file.path(.script_dir, "helpers", "cluster_marker_audit_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03c3_cluster_selection_report"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03c <- read_manifest_local(cfg$module_03c_manifest_path)
manifest_03c2 <- read_manifest_local(cfg$module_03c2_manifest_path)
candidate_rds <- resolve_output_local(manifest_03c, "candidate_clustered_object")
ranking_path <- resolve_output_local(manifest_03c, "resolution_ranking_tsv")
quality_path <- resolve_output_local(manifest_03c2, "candidate_marker_quality_tsv")

seu <- readRDS(candidate_rds)
seu <- maybe_join_layers(seu)
ranking <- read_tsv_optional(ranking_path)
quality <- read_tsv_optional(quality_path)
if (nrow(ranking) == 0 || nrow(quality) == 0) {
  stop("Missing clustering ranking or marker quality rows.", call. = FALSE)
}

selection_df <- dplyr::inner_join(
  ranking,
  quality[, setdiff(colnames(quality), intersect(colnames(quality), c("resolution", "candidate_cluster_col", "n_clusters"))), drop = FALSE],
  by = "candidate_id"
)
if (nrow(selection_df) == 0) {
  stop("No overlap between clustering candidates and marker-audited candidates.", call. = FALSE)
}
for (score_col in c("stability_score", "separation_score", "fragmentation_score", "marker_quality_score")) {
  selection_df[[score_col]] <- suppressWarnings(as.numeric(selection_df[[score_col]]))
  selection_df[[score_col]][!is.finite(selection_df[[score_col]])] <- 0.5
}
if (!"eligible" %in% colnames(selection_df)) {
  selection_df$eligible <- TRUE
}
if (!"ineligible_reason" %in% colnames(selection_df)) {
  selection_df$ineligible_reason <- ""
}
marker_supported <- suppressWarnings(as.numeric(selection_df$marker_supported_cluster_fraction))
low_marker_support <- is.finite(marker_supported) & marker_supported < 0.50
selection_df$ineligible_reason[low_marker_support] <- ifelse(
  nzchar(selection_df$ineligible_reason[low_marker_support]),
  paste(selection_df$ineligible_reason[low_marker_support], "marker_supported_clusters_lt_50pct", sep = ";"),
  "marker_supported_clusters_lt_50pct"
)
selection_df$eligible <- as.logical(selection_df$eligible) & !low_marker_support
selection_df$final_selection_score <- 0.50 * selection_df$stability_score +
  0.20 * selection_df$separation_score +
  0.20 * selection_df$marker_quality_score +
  0.10 * selection_df$fragmentation_score
if (any(selection_df$eligible)) {
  selection_df <- selection_df %>%
    dplyr::arrange(dplyr::desc(eligible), dplyr::desc(final_selection_score), dplyr::desc(stability_score), resolution)
} else {
  selection_df <- selection_df %>%
    dplyr::arrange(dplyr::desc(final_selection_score), dplyr::desc(stability_score), resolution)
}
selection_df$final_rank <- seq_len(nrow(selection_df))
selection_df$recommended_auto <- selection_df$final_rank == 1L

recommended_candidate_id <- selection_df$candidate_id[[which(selection_df$recommended_auto)[[1]]]]
decision_cols <- c("layer_id", "recommended_candidate_id", "approved_candidate_id", "status", "note")
selection_control <- read_tsv_optional(cfg$cluster_selection_file)
if (nrow(selection_control) == 0) {
  selection_control <- data.frame(
    layer_id = character(0),
    recommended_candidate_id = character(0),
    approved_candidate_id = character(0),
    status = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )
}
for (col in decision_cols) {
  if (!col %in% colnames(selection_control)) {
    selection_control[[col]] <- ""
  }
}
selection_control <- selection_control[, decision_cols, drop = FALSE]
row_idx <- which(selection_control$layer_id == panorama_spec$layer_id)
if (length(row_idx) == 0L) {
  selection_control <- rbind(
    selection_control,
    data.frame(
      layer_id = panorama_spec$layer_id,
      recommended_candidate_id = recommended_candidate_id,
      approved_candidate_id = "",
      status = "pending",
      note = "",
      stringsAsFactors = FALSE
    )
  )
  row_idx <- nrow(selection_control)
} else {
  row_idx <- row_idx[[1]]
  selection_control$recommended_candidate_id[[row_idx]] <- recommended_candidate_id
  if (!nzchar(normalize_scalar_value(selection_control$status[[row_idx]]))) {
    selection_control$status[[row_idx]] <- "pending"
  }
}

approved_candidate_id <- normalize_scalar_value(selection_control$approved_candidate_id[[row_idx]])
approval_status <- tolower(normalize_scalar_value(selection_control$status[[row_idx]], "pending"))
selected_candidate_id <- recommended_candidate_id
selection_source <- "automatic_provisional"
selection_status <- approval_status
if (identical(approval_status, "approved")) {
  if (nzchar(approved_candidate_id)) {
    if (!approved_candidate_id %in% selection_df$candidate_id) {
      stop(
        sprintf(
          "cluster_selection.tsv approved_candidate_id `%s` is not in current candidates: %s",
          approved_candidate_id,
          paste(selection_df$candidate_id, collapse = ",")
        ),
        call. = FALSE
      )
    }
    selected_candidate_id <- approved_candidate_id
    selection_source <- "manual_approved"
  } else {
    selection_source <- "manual_approved_recommended"
  }
} else if (nzchar(approved_candidate_id)) {
  selection_status <- paste0(approval_status, "_candidate_not_approved")
}
write_tsv_local(selection_control, cfg$cluster_selection_file)

selection_df$selected_final <- selection_df$candidate_id == selected_candidate_id
selection_df$selection_source <- ifelse(selection_df$selected_final, selection_source, "")
selection_df$selection_status <- ifelse(selection_df$selected_final, selection_status, "")

selected <- selection_df[selection_df$selected_final, , drop = FALSE]
source_cluster_col <- selected$candidate_cluster_col[[1]]
if (!source_cluster_col %in% colnames(seu@meta.data)) {
  stop(sprintf("Selected candidate cluster column missing: %s", source_cluster_col), call. = FALSE)
}

cluster_col <- paste0(panorama_spec$layer_id, "_cluster")
selected_reduction <- normalize_scalar_value(seu@misc$selected_reduction, "")
selected_umap <- normalize_scalar_value(seu@misc$selected_umap, "")
if (!nzchar(selected_reduction)) {
  selected_reduction <- "pca"
}
if (!nzchar(selected_umap)) {
  selected_umap <- "umap"
}
manifest_03a2 <- read_manifest_local(cfg$module_03a2_manifest_path)
candidate_index <- read_tsv_optional(resolve_output_local(manifest_03a2, "candidate_index_tsv"))
selected_integration <- normalize_scalar_value(seu@misc$selected_integration, "")
if (nrow(candidate_index) > 0 && nzchar(selected_integration)) {
  key <- paste0("candidate_", gsub("__", "_", selected_integration, fixed = TRUE))
  idx <- candidate_index[candidate_index$candidate_key == key, , drop = FALSE]
  if (nrow(idx) == 1) {
    selected_reduction <- idx$reduction_name[[1]]
    selected_umap <- idx$umap_name[[1]]
  }
}

seu <- finalize_layer_object(
  seu,
  reduction_name = selected_reduction,
  umap_name = selected_umap,
  cluster_col = cluster_col,
  source_cluster_col = source_cluster_col
)
seu$cell_type <- ""
seu$cell_type_confidence <- "undetermined"
seu$annotation_relation <- "none"
seu@misc$selected_resolution <- selected$resolution[[1]]
seu@misc$selected_cluster_col <- cluster_col
seu@misc$selected_candidate_cluster_col <- source_cluster_col
seu@misc$cluster_selection_score <- selected$final_selection_score[[1]]
seu@misc$cluster_selection_source <- selection_source
seu@misc$cluster_selection_status <- selection_status

ensure_dir(dirname(cfg$panorama_clustered_rds))
saveRDS(seu, cfg$panorama_clustered_rds)
ensure_dir(dirname(cfg$compat_clustered_rds))
if (!identical(normalizePath(cfg$panorama_clustered_rds, winslash = "/", mustWork = FALSE), normalizePath(cfg$compat_clustered_rds, winslash = "/", mustWork = FALSE))) {
  invisible(file.copy(cfg$panorama_clustered_rds, cfg$compat_clustered_rds, overwrite = TRUE))
}

final_audit <- run_cluster_marker_audit(
  seu,
  cluster_col = cluster_col,
  cfg = cfg,
  assay = cfg$cluster_marker_assay
)
marker_paths <- write_cluster_marker_audit_outputs(final_audit, cfg$cluster_marker_table_dir_layer)

cluster_summary <- build_cluster_summary_local(seu, cluster_col)
cluster_summary_csv <- file.path(cfg$integration_table_dir_layer, "cluster_summary.csv")
selection_tsv <- file.path(cfg$integration_table_dir_layer, "selected_resolution.tsv")
selection_ranking_tsv <- file.path(cfg$integration_table_dir_layer, "cluster_selection_ranking.tsv")
selected_resolution_txt <- file.path(cfg$integration_table_dir_layer, "selected_resolution.txt")
report_path <- file.path(cfg$cluster_selection_report_dir_layer, "report.md")
write_csv_local(cluster_summary, cluster_summary_csv)
write_tsv_local(selection_df, selection_ranking_tsv)
write_tsv_local(selected, selection_tsv)
writeLines(
  c(
    sprintf("selected_resolution=%s", selected$resolution[[1]]),
    sprintf("selected_cluster_count=%s", selected$n_clusters[[1]]),
    sprintf("selected_candidate_id=%s", selected$candidate_id[[1]]),
    sprintf("recommended_candidate_id=%s", recommended_candidate_id),
    sprintf("approved_candidate_id=%s", approved_candidate_id),
    sprintf("selection_source=%s", selection_source),
    sprintf("selection_status=%s", selection_status),
    sprintf("selected_cluster_column=%s", cluster_col),
    sprintf("source_candidate_cluster_column=%s", source_cluster_col),
    sprintf("final_selection_score=%s", selected$final_selection_score[[1]])
  ),
  selected_resolution_txt,
  useBytes = TRUE
)

status_df <- read_tsv_optional(cfg$layer_status_file)
status_row <- if (nrow(status_df) > 0 && any(status_df$layer_id == panorama_spec$layer_id)) {
  status_df[status_df$layer_id == panorama_spec$layer_id, , drop = FALSE][1, , drop = FALSE]
} else {
  data.frame(layer_id = panorama_spec$layer_id, stringsAsFactors = FALSE)
}
status_row$status <- "clustered"
status_row$cluster_column <- cluster_col
status_row$selected_resolution <- as.character(selected$resolution[[1]])
status_row$selected_cluster_count <- as.character(selected$n_clusters[[1]])
status_row$cluster_selection_status <- selection_status
status_row$cluster_selection_source <- selection_source
status_row$clustered_rds <- normalizePath(cfg$panorama_clustered_rds, winslash = "/", mustWork = FALSE)
status_row$annotated_rds <- ""
invisible(upsert_layer_status(cfg$layer_status_file, status_row))

report_lines <- c(
  "# 03c3 Cluster Selection Report",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- selected_resolution: `%s`", selected$resolution[[1]]),
  sprintf("- selected_cluster_count: `%s`", selected$n_clusters[[1]]),
  sprintf("- recommended_candidate_id: `%s`", recommended_candidate_id),
  sprintf("- approved_candidate_id: `%s`", ifelse(nzchar(approved_candidate_id), approved_candidate_id, "")),
  sprintf("- selection_source: `%s`", selection_source),
  sprintf("- cluster_selection_file: `%s`", cfg$cluster_selection_file),
  sprintf("- final_selection_score: `%s`", fmt_num(selected$final_selection_score[[1]], 3)),
  sprintf("- marker assay: `%s`", cfg$cluster_marker_assay),
  "",
  "## Selection Ranking",
  render_markdown_table_local(selection_df[, intersect(c(
    "candidate_id", "resolution", "n_clusters", "stability_score", "separation_score",
    "marker_quality_score", "fragmentation_score", "eligible", "ineligible_reason",
    "final_selection_score", "recommended_auto", "selected_final", "selection_source", "selection_status"
  ), colnames(selection_df)), drop = FALSE]),
  "",
  "## Final Cluster Marker QC",
  render_markdown_table_local(final_audit$cluster_qc[, intersect(c(
    "cluster", "n_cells", "raw_marker_n", "anno_marker_n", "strict_marker_n",
    "top20_mean_pct_diff", "top20_high_pct2_n", "deg_quality_flag", "risk_level"
  ), colnames(final_audit$cluster_qc)), drop = FALSE]),
  "",
  "## Key Outputs",
  sprintf("- clustered_object: `%s`", cfg$panorama_clustered_rds),
  sprintf("- markers_raw: `%s`", marker_paths$markers_raw),
  sprintf("- markers_strict: `%s`", marker_paths$markers_strict),
  sprintf("- cluster_marker_qc: `%s`", marker_paths$cluster_marker_qc),
  "",
  "Review the final cluster choice and marker audit before approving the `clustering` gate."
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03c3_manifest_path,
  new_outputs = list(
    clustered_object = build_output_entry(cfg$panorama_clustered_rds, "rds", module_name, "final selected panorama clustered Seurat object", base_dir = cfg$project_root),
    compatibility_clustered_object = build_output_entry(cfg$compat_clustered_rds, "rds", module_name, "compatibility checkpoint for downstream legacy modules", base_dir = cfg$project_root),
    cluster_selection_ranking_tsv = build_output_entry(selection_ranking_tsv, "tsv", module_name, "marker-assisted final resolution ranking", base_dir = cfg$project_root, schema = infer_schema_from_df(selection_df)),
    selected_resolution_tsv = build_output_entry(selection_tsv, "tsv", module_name, "selected final resolution row", base_dir = cfg$project_root, schema = infer_schema_from_df(selected)),
    selected_resolution_txt = build_output_entry(selected_resolution_txt, "txt", module_name, "selected final resolution summary", base_dir = cfg$project_root),
    cluster_selection_tsv = build_output_entry(cfg$cluster_selection_file, "tsv", module_name, "manual clustering selection control table", base_dir = cfg$project_root, schema = infer_schema_from_df(selection_control)),
    cluster_summary_csv = build_output_entry(cluster_summary_csv, "csv", module_name, "final cluster composition summary by sample/group", base_dir = cfg$project_root, schema = infer_schema_from_df(cluster_summary)),
    markers_raw_tsv = build_output_entry(marker_paths$markers_raw, "tsv", module_name, "panel-free FindAllMarkers output for final clustering", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$markers_raw)),
    markers_scored_tsv = build_output_entry(marker_paths$markers_scored, "tsv", module_name, "scored final cluster markers", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$markers_scored)),
    markers_anno_tsv = build_output_entry(marker_paths$markers_anno, "tsv", module_name, "annotation-ready final cluster markers", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$markers_anno)),
    markers_strict_tsv = build_output_entry(marker_paths$markers_strict, "tsv", module_name, "strict final cluster markers", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$markers_strict)),
    cluster_marker_qc_tsv = build_output_entry(marker_paths$cluster_marker_qc, "tsv", module_name, "final cluster DEG/marker quality QC", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$cluster_qc)),
    top20_markers_tsv = build_output_entry(marker_paths$top20_markers, "tsv", module_name, "top 20 final markers by specificity rank", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$top20_markers)),
    top50_markers_tsv = build_output_entry(marker_paths$top50_markers, "tsv", module_name, "top 50 final markers by specificity rank", base_dir = cfg$project_root, schema = infer_schema_from_df(final_audit$top50_markers)),
    report = build_output_entry(report_path, "md", module_name, "human-readable final clustering selection report", base_dir = cfg$project_root),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03c_manifest = cfg$module_03c_manifest_path,
    module_03c2_manifest = cfg$module_03c2_manifest_path,
    cluster_selection_file = cfg$cluster_selection_file,
    candidate_clustered_object = candidate_rds,
    resolution_ranking_tsv = ranking_path,
    candidate_marker_quality_tsv = quality_path
  ),
  version = cfg$module_version,
  depends_on = list(module_03c = cfg$module_03c_manifest_path, module_03c2 = cfg$module_03c2_manifest_path)
)

message("03c3 complete. Final clustered object: ", cfg$panorama_clustered_rds)
