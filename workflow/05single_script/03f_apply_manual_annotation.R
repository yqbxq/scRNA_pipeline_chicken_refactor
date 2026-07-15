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
source_utf8(file.path(.script_dir, "helpers", "clustering_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03f_apply_manual_annotation"
prepare_dirs_03(cfg)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03c3 <- read_manifest_local(cfg$module_03c3_manifest_path)
clustered_rds <- resolve_output_local(manifest_03c3, "clustered_object")
seu <- readRDS(clustered_rds)
seu <- maybe_join_layers(seu)

cluster_col <- paste0(panorama_spec$layer_id, "_cluster")
if (!cluster_col %in% colnames(seu@meta.data)) {
  stop(sprintf("Clustered object missing final cluster column: %s", cluster_col), call. = FALSE)
}
cluster_ids <- sort(unique(as.character(seu@meta.data[[cluster_col]])))
cluster_sizes <- as.data.frame(table(as.character(seu@meta.data[[cluster_col]])), stringsAsFactors = FALSE)
colnames(cluster_sizes) <- c("cluster", "n_cells")

template <- data.frame(
  layer_id = panorama_spec$layer_id,
  cluster = cluster_ids,
  manual_annotation = "",
  manual_annotation_zh = "",
  confidence = "",
  reason = "",
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(cluster_sizes, by = "cluster") %>%
  dplyr::select(layer_id, cluster, n_cells, manual_annotation, manual_annotation_zh, confidence, reason)

ensure_dir(dirname(cfg$manual_annotation_file))
manual <- read_tsv_optional(cfg$manual_annotation_file)
if (nrow(manual) == 0) {
  write_tsv_local(template, cfg$manual_annotation_file)
  stop(sprintf("Manual annotation template created: %s. Fill manual_annotation for every cluster and rerun 03_panorama.", cfg$manual_annotation_file), call. = FALSE)
}

required_cols <- c("layer_id", "cluster", "manual_annotation", "manual_annotation_zh", "confidence", "reason")
for (col in required_cols) {
  if (!col %in% colnames(manual)) {
    manual[[col]] <- ""
  }
}
if (!"n_cells" %in% colnames(manual)) {
  manual$n_cells <- NA_integer_
}
manual <- manual[, unique(c(required_cols, "n_cells")), drop = FALSE]
for (col in required_cols) {
  manual[[col]] <- vapply(manual[[col]], normalize_scalar_value, character(1))
}
manual <- manual[manual$layer_id %in% c(panorama_spec$layer_id, "*", ""), , drop = FALSE]
manual$layer_id[!nzchar(manual$layer_id) | manual$layer_id == "*"] <- panorama_spec$layer_id
manual <- manual[!duplicated(manual[, c("layer_id", "cluster")]), , drop = FALSE]

missing_clusters <- setdiff(cluster_ids, manual$cluster)
blank_clusters <- manual$cluster[manual$cluster %in% cluster_ids & !nzchar(manual$manual_annotation)]
extra_clusters <- setdiff(manual$cluster, cluster_ids)
if (length(missing_clusters) > 0 || length(blank_clusters) > 0 || length(extra_clusters) > 0) {
  merged_template <- template
  if (nrow(manual) > 0) {
    keep_cols <- intersect(c("cluster", "manual_annotation", "manual_annotation_zh", "confidence", "reason"), colnames(manual))
    merged_template <- merged_template %>%
      dplyr::select(-manual_annotation, -manual_annotation_zh, -confidence, -reason) %>%
      dplyr::left_join(manual[, keep_cols, drop = FALSE], by = "cluster")
    for (col in c("manual_annotation", "manual_annotation_zh", "confidence", "reason")) {
      if (!col %in% colnames(merged_template)) {
        merged_template[[col]] <- ""
      }
      merged_template[[col]][is.na(merged_template[[col]])] <- ""
    }
    merged_template <- merged_template[, colnames(template), drop = FALSE]
  }
  write_tsv_local(merged_template, cfg$manual_annotation_file)
  stop(sprintf(
    "Manual annotation table is incomplete: missing=%s blank=%s extra=%s. Updated template: %s",
    paste(missing_clusters, collapse = ","),
    paste(blank_clusters, collapse = ","),
    paste(extra_clusters, collapse = ","),
    cfg$manual_annotation_file
  ), call. = FALSE)
}

manual <- manual[match(cluster_ids, manual$cluster), , drop = FALSE]
manual$n_cells <- cluster_sizes$n_cells[match(manual$cluster, cluster_sizes$cluster)]
manual$confidence[!nzchar(manual$confidence)] <- "medium"

cluster_key <- as.character(seu@meta.data[[cluster_col]])
map_annotation <- setNames(manual$manual_annotation, manual$cluster)
map_annotation_zh <- setNames(manual$manual_annotation_zh, manual$cluster)
map_confidence <- setNames(manual$confidence, manual$cluster)
map_reason <- setNames(manual$reason, manual$cluster)

seu$annotation_label <- unname(map_annotation[cluster_key])
seu$cell_type <- seu$annotation_label
seu$cell_type_confidence <- unname(map_confidence[cluster_key])
seu$annotation_confidence <- seu$cell_type_confidence
seu$annotation_relation <- "manual"
seu$annotation_source <- "manual"
seu$annotation_reason <- unname(map_reason[cluster_key])
seu$manual_annotation_zh <- unname(map_annotation_zh[cluster_key])
seu$cell_type_zh <- seu$manual_annotation_zh

ensure_dir(dirname(cfg$panorama_annotated_rds))
saveRDS(seu, cfg$panorama_annotated_rds)
ensure_dir(dirname(cfg$compat_annotated_rds))
if (!identical(normalizePath(cfg$panorama_annotated_rds, winslash = "/", mustWork = FALSE), normalizePath(cfg$compat_annotated_rds, winslash = "/", mustWork = FALSE))) {
  invisible(file.copy(cfg$panorama_annotated_rds, cfg$compat_annotated_rds, overwrite = TRUE))
}

annotation_table <- manual %>%
  dplyr::mutate(
    cluster_id = cluster,
    cluster_column = cluster_col,
    annotation_source = "manual",
    final_annotation = manual_annotation,
    candidate_celltype = manual_annotation,
    relation = "manual",
    evidence_relation = "manual"
  )
annotation_summary <- annotation_table %>%
  dplyr::count(manual_annotation, confidence, name = "cluster_n") %>%
  dplyr::arrange(manual_annotation, confidence)

annotation_table_tsv <- file.path(cfg$annotation_table_dir_layer, "manual_annotation_table.tsv")
annotation_summary_tsv <- file.path(cfg$annotation_table_dir_layer, "manual_annotation_summary.tsv")
report_path <- file.path(cfg$annotation_report_dir_layer, "manual_annotation_apply.md")
write_tsv_local(annotation_table, annotation_table_tsv)
write_tsv_local(annotation_summary, annotation_summary_tsv)

status_df <- read_tsv_optional(cfg$layer_status_file)
status_row <- if (nrow(status_df) > 0 && any(status_df$layer_id == panorama_spec$layer_id)) {
  status_df[status_df$layer_id == panorama_spec$layer_id, , drop = FALSE][1, , drop = FALSE]
} else {
  data.frame(layer_id = panorama_spec$layer_id, stringsAsFactors = FALSE)
}
status_row$status <- "annotated"
status_row$annotated_rds <- normalizePath(cfg$panorama_annotated_rds, winslash = "/", mustWork = FALSE)
status_row$cluster_column <- cluster_col
invisible(upsert_layer_status(cfg$layer_status_file, status_row))

report_lines <- c(
  "# 03f Apply Manual Annotation",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- manual_annotation_file: `%s`", cfg$manual_annotation_file),
  sprintf("- annotated_object: `%s`", cfg$panorama_annotated_rds),
  sprintf("- cluster_count: `%s`", nrow(annotation_table)),
  "",
  "## Manual Annotation Table",
  render_markdown_table_local(annotation_table[, intersect(c(
    "cluster", "n_cells", "manual_annotation", "manual_annotation_zh", "confidence", "reason"
  ), colnames(annotation_table)), drop = FALSE])
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03f_manifest_path,
  new_outputs = list(
    annotated_object = build_output_entry(cfg$panorama_annotated_rds, "rds", module_name, "manually annotated panorama Seurat object", base_dir = cfg$project_root),
    compatibility_annotated_object = build_output_entry(cfg$compat_annotated_rds, "rds", module_name, "compatibility checkpoint for downstream legacy modules", base_dir = cfg$project_root),
    manual_annotation_file = build_output_entry(cfg$manual_annotation_file, "tsv", module_name, "human-edited manual annotation input", base_dir = cfg$project_root),
    annotation_table_tsv = build_output_entry(annotation_table_tsv, "tsv", module_name, "one row per manually annotated cluster", base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_table)),
    annotation_summary_tsv = build_output_entry(annotation_summary_tsv, "tsv", module_name, "manual annotation cluster counts", base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_summary)),
    report = build_output_entry(report_path, "md", module_name, "manual annotation apply report", base_dir = cfg$project_root),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", module_name, "one row per built/annotated object layer", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03c3_manifest = cfg$module_03c3_manifest_path,
    module_03d_marker_risk_manifest = cfg$module_03d_marker_risk_manifest_path,
    module_03e_panel_evidence_manifest = cfg$module_03e_panel_evidence_manifest_path,
    manual_annotation_file = cfg$manual_annotation_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03c3 = cfg$module_03c3_manifest_path, module_03d_marker_risk = cfg$module_03d_marker_risk_manifest_path, module_03e_panel_evidence = cfg$module_03e_panel_evidence_manifest_path)
)

write_manifest_local(
  manifest_path = cfg$module_03d_manifest_path,
  new_outputs = list(
    annotated_object = build_output_entry(cfg$panorama_annotated_rds, "rds", "03f_apply_manual_annotation", "manually annotated panorama Seurat object", base_dir = cfg$project_root),
    compatibility_annotated_object = build_output_entry(cfg$compat_annotated_rds, "rds", "03f_apply_manual_annotation", "compatibility checkpoint for downstream legacy modules", base_dir = cfg$project_root),
    annotation_table_tsv = build_output_entry(annotation_table_tsv, "tsv", "03f_apply_manual_annotation", "one row per manually annotated cluster", base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_table)),
    annotation_summary_tsv = build_output_entry(annotation_summary_tsv, "tsv", "03f_apply_manual_annotation", "manual annotation cluster counts", base_dir = cfg$project_root, schema = infer_schema_from_df(annotation_summary)),
    report = build_output_entry(report_path, "md", "03f_apply_manual_annotation", "manual annotation apply report", base_dir = cfg$project_root),
    layer_status_tsv = build_output_entry(cfg$layer_status_file, "tsv", "03f_apply_manual_annotation", "one row per built/annotated object layer", base_dir = cfg$project_root)
  ),
  module_name = "03d_annotate",
  base_dir = cfg$project_root,
  inputs = list(
    module_03f_manifest = cfg$module_03f_manifest_path,
    manual_annotation_file = cfg$manual_annotation_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03f = cfg$module_03f_manifest_path)
)

message("03f complete. Annotated object: ", cfg$panorama_annotated_rds)
