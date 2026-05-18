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
module_name <- "spatial_04c_subcluster_eda"
prepare_dirs_spatial(cfg)

assignment_schema <- function() {
  data.frame(
    layer_id = character(),
    parent_region = character(),
    cluster = character(),
    sub_region = character(),
    n_spots = integer(),
    module_score_max = numeric(),
    panel_hit_count_winner = integer(),
    evidence_agreement = logical(),
    confidence = character(),
    stringsAsFactors = FALSE
  )
}

plot_spatial_subregion <- function(obj, path) {
  ensure_dir(dirname(path))
  coords <- spatial_filter_coords(obj)
  coords <- coords[intersect(rownames(coords), rownames(obj@meta.data)), , drop = FALSE]
  coords$sub_region <- as.character(obj@meta.data[rownames(coords), "sub_region", drop = TRUE])
  section_col <- if ("section_id" %in% colnames(obj@meta.data)) "section_id" else if ("spatial_section_id" %in% colnames(obj@meta.data)) "spatial_section_id" else ""
  if (nzchar(section_col)) {
    coords$section_id <- as.character(obj@meta.data[rownames(coords), section_col, drop = TRUE])
  }
  p <- ggplot2::ggplot(coords, ggplot2::aes(x = col, y = row, color = sub_region)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.9, na.rm = TRUE) +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_fixed() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "sub_region")
  if ("section_id" %in% colnames(coords)) {
    p <- p + ggplot2::facet_wrap(~section_id)
  }
  ggplot2::ggsave(path, p, width = 6.4, height = 5.2, dpi = 180, bg = "white")
  path
}

plot_umap_subregion <- function(obj, path) {
  ensure_dir(dirname(path))
  emb <- Seurat::Embeddings(obj, reduction = "umap_subcluster")
  plot_df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], sub_region = as.character(obj@meta.data[rownames(emb), "sub_region", drop = TRUE]))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = sub_region)) +
    ggplot2::geom_point(size = 0.9, alpha = 0.9, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(color = "sub_region")
  ggplot2::ggsave(path, p, width = 5.8, height = 4.8, dpi = 180, bg = "white")
  path
}

plot_module_score_boxplot <- function(obj, path) {
  ensure_dir(dirname(path))
  score_cols <- grep("^ms_", colnames(obj@meta.data), value = TRUE)
  if (length(score_cols) == 0) {
    plot_df <- data.frame(sub_region = as.character(obj$sub_region), module = "none", score = 0, stringsAsFactors = FALSE)
  } else {
    plot_df <- do.call(rbind, lapply(score_cols, function(col) {
      data.frame(sub_region = as.character(obj$sub_region), module = col, score = as.numeric(obj@meta.data[[col]]), stringsAsFactors = FALSE)
    }))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = sub_region, y = score, fill = sub_region)) +
    ggplot2::geom_boxplot(outlier.size = 0.2, na.rm = TRUE) +
    ggplot2::facet_wrap(~module, scales = "free_y") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = "module score")
  ggplot2::ggsave(path, p, width = 7.0, height = 4.8, dpi = 180, bg = "white")
  path
}

plot_marker_heatmap <- function(module_score_matrix, path) {
  ensure_dir(dirname(path))
  value_cols <- setdiff(colnames(module_score_matrix), "cluster")
  if (length(value_cols) == 0) {
    plot_df <- data.frame(cluster = "none", module = "none", score = 0, stringsAsFactors = FALSE)
  } else {
    plot_df <- do.call(rbind, lapply(value_cols, function(col) {
      data.frame(cluster = module_score_matrix$cluster, module = col, score = as.numeric(module_score_matrix[[col]]), stringsAsFactors = FALSE)
    }))
  }
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = module, y = cluster, fill = score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(x = "module", y = "cluster", fill = "score")
  ggplot2::ggsave(path, p, width = 6.2, height = 4.6, dpi = 180, bg = "white")
  path
}

section_coherence_table <- function(obj) {
  section_col <- if ("section_id" %in% colnames(obj@meta.data)) "section_id" else if ("spatial_section_id" %in% colnames(obj@meta.data)) "spatial_section_id" else ""
  groups <- if (nzchar(section_col)) split(rownames(obj@meta.data), as.character(obj@meta.data[[section_col]])) else list(panorama = rownames(obj@meta.data))
  rows <- lapply(names(groups), function(section_id) {
    cells <- groups[[section_id]]
    score <- if (length(cells) >= 3L) {
      sub <- tryCatch(Seurat::subset(obj, cells = cells), error = function(e) NULL)
      if (is.null(sub)) NA_real_ else compute_spatial_coherence_score(sub, "sub_region", k = 6L)
    } else {
      NA_real_
    }
    data.frame(section_id = section_id, spatial_coherence_score = score, n_spots = length(cells), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

plot_section_coherence <- function(coherence_df, path) {
  ensure_dir(dirname(path))
  p <- ggplot2::ggplot(coherence_df, ggplot2::aes(x = section_id, y = spatial_coherence_score, fill = section_id)) +
    ggplot2::geom_col(width = 0.7, na.rm = TRUE) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = "spatial coherence")
  ggplot2::ggsave(path, p, width = 5.8, height = 4.2, dpi = 180, bg = "white")
  path
}

index_path <- file.path(cfg$spatial_subcluster_build_table_dir, "region_subset_index.tsv")
assignment_path <- file.path(cfg$spatial_subcluster_annotate_table_dir, "region_subset_assignment.tsv")
evidence_path <- file.path(cfg$spatial_subcluster_annotate_table_dir, "evidence_summary.tsv")
index_df <- spatial_read_tsv(index_path)
assignment_df <- spatial_read_tsv(assignment_path)
evidence_summary <- spatial_read_tsv(evidence_path)
ok_layers <- if (nrow(evidence_summary) > 0 && "status" %in% colnames(evidence_summary)) evidence_summary$layer_id[evidence_summary$status == "ok"] else character(0)
ok_layers <- ok_layers[nzchar(ok_layers)]

layer_results <- list()
combined_lines <- c("# Spatial Subregion Annotation EDA", "")
if (length(ok_layers) == 0) {
  if (nrow(evidence_summary) > 0 && "status" %in% colnames(evidence_summary)) {
    triage_df <- data.frame(
      layer_id = as.character(evidence_summary$layer_id),
      signal = as.character(evidence_summary$status),
      status = ifelse(evidence_summary$status %in% c("no_subsets_built", "skipped_no_enabled_layers", "skipped"), "skipped", "warn"),
      value = as.character(evidence_summary$n_assigned %||% ""),
      threshold = "",
      message = as.character(evidence_summary$reason %||% "no ok subcluster annotation layers are available"),
      stringsAsFactors = FALSE
    )
  } else {
    triage_df <- data.frame(layer_id = "", signal = "no_enabled_region_subsets", status = "skipped", value = "", threshold = "", message = "no ok subcluster annotation layers are available", stringsAsFactors = FALSE)
  }
  spatial_write_tsv(triage_df, file.path(cfg$spatial_subcluster_eda_table_dir, "triage.tsv"))
  report_path <- file.path(cfg$spatial_subregion_annotation_dir, "report.md")
  write_markdown_local(c(combined_lines, "No enabled region subsets were available for sub-region EDA.", "", "## Triage", render_markdown_table_local(triage_df)), report_path)
  write_manifest_local(
    manifest_path = cfg$module_04c_subcluster_eda_manifest_path,
    new_outputs = list(
      triage = build_output_entry(file.path(cfg$spatial_subcluster_eda_table_dir, "triage.tsv"), "tsv", module_name, "subcluster triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
      report = build_output_entry(report_path, "md", module_name, "combined subregion EDA report", base_dir = cfg$project_root)
    ),
    module_name = module_name,
    base_dir = cfg$project_root,
    inputs = list(region_subset_index = index_path, region_subset_assignment = assignment_path),
    version = cfg$module_04_version,
    depends_on = list(spatial_04b_subcluster_annotate = cfg$module_04b_subcluster_annotate_manifest_path)
  )
  message("spatial subcluster EDA skipped: no ok layers")
  quit(status = 0)
}

for (layer_id in ok_layers) {
  layer_report_dir <- file.path(cfg$spatial_subregion_annotation_dir, layer_id)
  layer_figure_dir <- file.path(cfg$spatial_subcluster_eda_figure_dir, layer_id)
  layer_table_dir <- file.path(cfg$spatial_subcluster_eda_table_dir, layer_id)
  ensure_dir(layer_report_dir)
  ensure_dir(layer_figure_dir)
  ensure_dir(layer_table_dir)
  subset_rds <- file.path(cfg$spatial_region_annotated_dir, sprintf("%s.rds", spatial_safe_id(layer_id)))
  if (!file.exists(subset_rds)) {
    layer_results[[layer_id]] <- list(status = "failed", message = "annotated subset RDS is missing", assignment = assignment_schema())
    next
  }
  obj <- readRDS(subset_rds)
  layer_assignment <- assignment_df[assignment_df$layer_id == layer_id, , drop = FALSE]
  module_score_matrix <- spatial_read_tsv(file.path(cfg$spatial_subcluster_annotate_table_dir, layer_id, "module_score_matrix.tsv"))
  if (nrow(module_score_matrix) == 0) {
    module_score_matrix <- data.frame(cluster = character(), stringsAsFactors = FALSE)
  }
  coherence_df <- section_coherence_table(obj)
  overlap_df <- data.frame(
    layer_id = layer_id,
    cluster = as.character(layer_assignment$cluster),
    sub_region = as.character(layer_assignment$sub_region),
    panel_hit_count_winner = as.integer(layer_assignment$panel_hit_count_winner),
    has_panel_overlap = as.integer(layer_assignment$panel_hit_count_winner) > 0,
    stringsAsFactors = FALSE
  )
  spatial_write_tsv(coherence_df, file.path(layer_table_dir, "coherence_per_section.tsv"))
  spatial_write_tsv(overlap_df, file.path(layer_table_dir, "panel_overlap.tsv"))
  figures <- c(
    umap_subregion = plot_umap_subregion(obj, file.path(layer_figure_dir, "umap_subregion.png")),
    spatial_subregion = plot_spatial_subregion(obj, file.path(layer_figure_dir, "spatial_subregion.png")),
    marker_heatmap = plot_marker_heatmap(module_score_matrix, file.path(layer_figure_dir, "marker_heatmap.png")),
    coherence_per_section = plot_section_coherence(coherence_df, file.path(layer_figure_dir, "coherence_per_section.png")),
    module_score_boxplot = plot_module_score_boxplot(obj, file.path(layer_figure_dir, "module_score_boxplot.png"))
  )
  layer_results[[layer_id]] <- list(status = "ok", obj = obj, assignment = layer_assignment, figures = figures)
  layer_report <- c(
    sprintf("# Spatial Subregion EDA: %s", layer_id),
    "",
    sprintf("- Annotated subset: `%s`", relative_path_local(subset_rds, cfg$project_root)),
    sprintf("- Figures: `%s`", relative_path_local(layer_figure_dir, cfg$project_root)),
    "",
    "## Assignment",
    render_markdown_table_local(layer_assignment),
    "",
    "## Section Coherence",
    render_markdown_table_local(coherence_df),
    "",
    "## Panel Overlap",
    render_markdown_table_local(overlap_df)
  )
  write_markdown_local(layer_report, file.path(layer_report_dir, "report.md"))
}

triage_df <- summarize_subcluster_triage(layer_results, cfg)
triage_path <- file.path(cfg$spatial_subcluster_eda_table_dir, "triage.tsv")
spatial_write_tsv(triage_df, triage_path)
failed_layers <- names(layer_results)[vapply(layer_results, function(x) !identical(x$status, "ok"), logical(1))]
combined_lines <- c(
  combined_lines,
  sprintf("- Region subset index: `%s`", relative_path_local(index_path, cfg$project_root)),
  sprintf("- Assignment table: `%s`", relative_path_local(assignment_path, cfg$project_root)),
  sprintf("- Triage table: `%s`", relative_path_local(triage_path, cfg$project_root)),
  if (length(failed_layers) > 0) sprintf("- Failed layers: `%s`", paste(failed_layers, collapse = ", ")) else "- Failed layers: none",
  "",
  "## Evidence Summary",
  render_markdown_table_local(evidence_summary),
  "",
  "## Triage",
  render_markdown_table_local(triage_df),
  "",
  "## Assignment",
  render_markdown_table_local(assignment_df)
)
report_path <- file.path(cfg$spatial_subregion_annotation_dir, "report.md")
write_markdown_local(combined_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_04c_subcluster_eda_manifest_path,
  new_outputs = list(
    triage = build_output_entry(triage_path, "tsv", module_name, "subcluster triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    report = build_output_entry(report_path, "md", module_name, "combined subregion EDA report", base_dir = cfg$project_root),
    figures_dir = build_output_entry(cfg$spatial_subcluster_eda_figure_dir, "directory", module_name, "subregion EDA figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(region_subset_index = index_path, region_subset_assignment = assignment_path),
  version = cfg$module_04_version,
  depends_on = list(spatial_04b_subcluster_annotate = cfg$module_04b_subcluster_annotate_manifest_path)
)

message(sprintf("spatial subcluster EDA complete: %d ok layer(s)", length(ok_layers)))
