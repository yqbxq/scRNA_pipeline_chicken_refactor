empty_spatial_marker_table <- function() {
  data.frame(gene = character(), cluster = character(), group_value = character(), stringsAsFactors = FALSE)
}

empty_spatial_de_table <- function() {
  data.frame(gene = character(), avg_log2FC = numeric(), p_val = numeric(), p_val_adj = numeric(), stringsAsFactors = FALSE)
}

empty_spatial_formal_table <- function() {
  data.frame(gene = character(), logFC = numeric(), PValue = numeric(), FDR = numeric(), stringsAsFactors = FALSE)
}

empty_spatial_manifest_05 <- function() {
  data.frame(
    comparison_id = character(),
    layer_id = character(),
    group_by = character(),
    parent_region = character(),
    region_value = character(),
    analysis_mode = character(),
    analysis_unit = character(),
    gene_program_role = character(),
    group_var = character(),
    subset_column = character(),
    subset_value = character(),
    force_exploratory = character(),
    min_cells_per_group = integer(),
    min_biological_replicates = integer(),
    logfc_threshold = numeric(),
    inference_status = character(),
    replicate_gate_status = character(),
    status = character(),
    reason = character(),
    exploratory_results_tsv = character(),
    exploratory_summary_tsv = character(),
    formal_results_tsv = character(),
    exploratory_fallback = character(),
    stringsAsFactors = FALSE
  )
}

spatial_manifest_row_05 <- function(vars, layer_id, group_by, region_value = "", parent_region = "", inference_status = "exploratory", replicate_gate_status = "", status = "ok", reason = "", exploratory_results_tsv = "", exploratory_summary_tsv = "", formal_results_tsv = "", exploratory_fallback = "", cfg = NULL) {
  base_dir <- if (is.null(cfg)) "" else cfg$project_root
  rel <- function(path) {
    if (!nzchar(path) || is.null(cfg)) "" else relative_path_local(path, base_dir)
  }
  data.frame(
    comparison_id = vars$comparison_id,
    layer_id = layer_id,
    group_by = group_by,
    parent_region = parent_region,
    region_value = region_value,
    analysis_mode = vars$analysis_mode,
    analysis_unit = vars$analysis_unit,
    gene_program_role = vars$gene_program_role,
    group_var = vars$group_var,
    subset_column = vars$subset_column,
    subset_value = vars$subset_value,
    force_exploratory = ifelse(vars$force_exploratory, "yes", "no"),
    min_cells_per_group = vars$min_cells_per_group,
    min_biological_replicates = vars$min_biological_replicates,
    logfc_threshold = vars$logfc_threshold,
    inference_status = inference_status,
    replicate_gate_status = replicate_gate_status,
    status = status,
    reason = reason,
    exploratory_results_tsv = rel(exploratory_results_tsv),
    exploratory_summary_tsv = rel(exploratory_summary_tsv),
    formal_results_tsv = rel(formal_results_tsv),
    exploratory_fallback = rel(exploratory_fallback),
    stringsAsFactors = FALSE
  )
}

plot_spatial_marker_volcano <- function(markers, group_value, out_path) {
  ensure_dir(dirname(out_path))
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(markers) == 0) {
    return("")
  }
  df <- markers
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(df))[1]
  p_col <- intersect(c("p_val_adj", "p_val"), colnames(df))[1]
  if (is.na(fc_col) || is.na(p_col)) {
    return("")
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[fc_col]], y = -log10(pmax(.data[[p_col]], 1e-300)))) +
    ggplot2::geom_point(size = 0.8, alpha = 0.75) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = fc_col, y = paste0("-log10(", p_col, ")"), title = group_value)
  ggplot2::ggsave(out_path, p, width = 5.2, height = 4.2, dpi = 180, bg = "white")
  out_path
}

plot_spatial_top_markers_heatmap <- function(markers, out_path, top_n = 20L) {
  ensure_dir(dirname(out_path))
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(markers) == 0 || !"group_value" %in% colnames(markers)) {
    return("")
  }
  fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))[1]
  if (is.na(fc_col)) {
    return("")
  }
  df <- markers[order(markers$group_value, -as.numeric(markers[[fc_col]])), , drop = FALSE]
  df <- do.call(rbind, lapply(split(df, df$group_value), function(x) utils::head(x, top_n)))
  if (is.null(df) || nrow(df) == 0) {
    return("")
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group_value, y = gene, fill = .data[[fc_col]])) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, fill = fc_col)
  ggplot2::ggsave(out_path, p, width = 7.2, height = 6.0, dpi = 180, bg = "white")
  out_path
}

plot_spatial_composition_stackbar <- function(prop_df, out_path) {
  ensure_dir(dirname(out_path))
  if (!requireNamespace("ggplot2", quietly = TRUE) || nrow(prop_df) == 0) {
    return("")
  }
  p <- ggplot2::ggplot(prop_df, ggplot2::aes(x = section_id, y = proportion, fill = group_value)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~condition, scales = "free_x") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "spot proportion", fill = "group")
  ggplot2::ggsave(out_path, p, width = 6.8, height = 4.8, dpi = 180, bg = "white")
  out_path
}

write_spatial_marker_outputs <- function(out_dir, fig_dir, markers, summary, vars) {
  ensure_dir(out_dir)
  ensure_dir(fig_dir)
  if (nrow(markers) == 0) {
    markers <- empty_spatial_marker_table()
  }
  markers_path <- file.path(out_dir, "markers.tsv")
  summary_path <- file.path(out_dir, "summary.tsv")
  spatial_write_tsv(markers, markers_path)
  spatial_write_tsv(summary, summary_path)
  plot_spatial_top_markers_heatmap(markers, file.path(fig_dir, "top_markers_heatmap.png"), vars$top_n)
  if (nrow(markers) > 0 && "group_value" %in% colnames(markers)) {
    for (group_value in unique(as.character(markers$group_value))) {
      group_df <- markers[as.character(markers$group_value) == group_value, , drop = FALSE]
      plot_spatial_marker_volcano(group_df, group_value, file.path(fig_dir, sprintf("marker_volcano_%s.png", spatial_safe_id(group_value))))
    }
  }
  list(markers = markers_path, summary = summary_path)
}

write_spatial_pseudobulk_outputs <- function(out_dir, formal_res, status, reason, exploratory_fallback_path = "") {
  ensure_dir(out_dir)
  if (identical(status, "ok") && nrow(formal_res) > 0) {
    path <- file.path(out_dir, "formal_de.tsv")
    spatial_write_tsv(formal_res, path)
  } else {
    path <- file.path(out_dir, "exploratory_only.tsv")
    spatial_write_tsv(data.frame(status = status, reason = reason, exploratory_fallback = exploratory_fallback_path, stringsAsFactors = FALSE), path)
  }
  path
}

write_spatial_composition_outputs <- function(out_dir, fig_dir, prop_df, formal_res, status, reason) {
  ensure_dir(out_dir)
  ensure_dir(fig_dir)
  prop_path <- file.path(out_dir, "proportion.tsv")
  formal_path <- file.path(out_dir, "formal_propeller.tsv")
  spatial_write_tsv(prop_df, prop_path)
  if (identical(status, "ok") && nrow(formal_res) > 0) {
    spatial_write_tsv(formal_res, formal_path)
  } else {
    spatial_write_tsv(data.frame(status = status, reason = reason, stringsAsFactors = FALSE), formal_path)
  }
  plot_spatial_composition_stackbar(prop_df, file.path(fig_dir, "composition_stackbar.png"))
  list(proportion = prop_path, formal = formal_path)
}

write_spatial_de_outputs <- function(out_dir, fig_dir, de_res, vars, status, reason, region_value) {
  ensure_dir(out_dir)
  ensure_dir(fig_dir)
  path <- file.path(out_dir, "spatial_de_spotlevel.tsv")
  if (identical(status, "ok") && nrow(de_res) > 0) {
    spatial_write_tsv(de_res, path)
    plot_spatial_marker_volcano(de_res, region_value, file.path(fig_dir, "de_volcano.png"))
  } else {
    spatial_write_tsv(data.frame(status = status, reason = reason, stringsAsFactors = FALSE), path)
  }
  path
}

summarize_spatial_de_manifest <- function(manifest_df, cfg, report_path) {
  status_summary <- as.data.frame(table(status = manifest_df$status), stringsAsFactors = FALSE)
  colnames(status_summary) <- c("status", "n")
  lines <- c(
    "# Spatial Marker / DE Summary",
    "",
    "## Status Summary",
    render_markdown_table_local(status_summary),
    "",
    "## Manifest Preview",
    render_markdown_table_local(utils::head(manifest_df, 50))
  )
  write_markdown_local(lines, report_path)
  report_path
}
