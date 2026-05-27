communication_report_plot_path <- function(cfg, filename) {
  dir <- cfg$communication_eda_panel_dir %||% file.path(cfg$communication_consensus_figure_dir, "07e_panels")
  if (exists("ensure_dir", mode = "function")) ensure_dir(dir) else dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, filename)
}

communication_report_write_empty_png <- function(path, title = "No data") {
  if (exists("ensure_dir", mode = "function")) ensure_dir(dirname(path)) else dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = 900, height = 650)
  on.exit(dev.off(), add = TRUE)
  plot.new()
  title(main = title)
  text(0.5, 0.5, "No rows available")
  invisible(path)
}

communication_report_barplot_png <- function(values, path, title, ylab = "Count", colors = NULL) {
  if (exists("ensure_dir", mode = "function")) ensure_dir(dirname(path)) else dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  png(path, width = 1000, height = 700)
  on.exit(dev.off(), add = TRUE)
  if (length(values) == 0 || sum(values, na.rm = TRUE) == 0) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, "No rows available")
  } else {
    barplot(values, las = 2, main = title, ylab = ylab, col = colors %||% "steelblue")
  }
  invisible(path)
}

render_panel_tier_distribution <- function(consensus_df, cfg) {
  df <- communication_report_normalize_consensus(consensus_df)
  tiers <- c("primary", "exploratory", "candidate", "blocked")
  counts <- table(factor(df$evidence_tier, levels = tiers))
  tier_counts <- data.frame(
    evidence_tier = tiers,
    n = as.integer(counts),
    pct = if (sum(counts) > 0) round(100 * as.integer(counts) / sum(counts), 1) else 0,
    stringsAsFactors = FALSE
  )
  pie_path <- communication_report_plot_path(cfg, "panel_01_tier_distribution_pie.png")
  heat_path <- communication_report_plot_path(cfg, "panel_01_tier_distribution_heatmap.png")
  colors <- communication_report_tier_colors()[tiers]

  if (sum(counts) == 0) {
    communication_report_write_empty_png(pie_path, "Evidence tier distribution")
  } else {
    png(pie_path, width = 900, height = 650)
    pie(as.integer(counts), labels = paste0(tiers, " (", as.integer(counts), ")"), col = colors, main = "Evidence tier distribution")
    dev.off()
  }

  if (nrow(df) == 0) {
    communication_report_write_empty_png(heat_path, "Tier count by pair and condition")
  } else {
    pair_condition <- paste(df$pair_id, df$condition_value, sep = " | ")
    mat <- table(pair_condition, factor(df$evidence_tier, levels = tiers))
    png(heat_path, width = 1100, height = 800)
    image(
      t(as.matrix(mat[nrow(mat):1, , drop = FALSE])),
      axes = FALSE,
      col = colorRampPalette(c("white", "#2C7FB8"))(20),
      main = "Tier count by pair and condition"
    )
    axis(1, at = seq(0, 1, length.out = nrow(mat)), labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
    axis(2, at = seq(0, 1, length.out = length(tiers)), labels = tiers, las = 2)
    dev.off()
  }

  list(
    name = "tier_distribution",
    status = "ok",
    md_lines = c(
      "## Panel 1: Evidence Tier Distribution",
      "",
      communication_report_md_table(tier_counts, max_rows = 10L),
      sprintf("![tier pie](%s)", pie_path),
      sprintf("![tier heatmap](%s)", heat_path)
    ),
    plots = list(pie = pie_path, heatmap = heat_path),
    tsv = tier_counts
  )
}

render_panel_method_venn <- function(consensus_df, cfg) {
  df <- communication_report_normalize_consensus(consensus_df)
  if (!"methods_hit_pattern" %in% colnames(df)) {
    df$methods_hit_pattern <- communication_report_methods_pattern(df)
  }
  pattern_counts <- as.data.frame(table(df$methods_hit_pattern), stringsAsFactors = FALSE)
  names(pattern_counts) <- c("methods_hit_pattern", "n")
  pattern_counts <- pattern_counts[order(-pattern_counts$n, pattern_counts$methods_hit_pattern), , drop = FALSE]
  plot_path <- communication_report_plot_path(cfg, "panel_02_method_agreement.png")
  values <- setNames(pattern_counts$n, pattern_counts$methods_hit_pattern)
  communication_report_barplot_png(values, plot_path, "Method agreement pattern", ylab = "Axis count")
  cellchat_only <- sum(df$methods_hit_pattern == "cellchat", na.rm = TRUE)
  list(
    name = "method_agreement",
    status = "ok",
    md_lines = c(
      "## Panel 2: Method Agreement",
      "",
      sprintf("- CellChat-only candidate axes: %d", cellchat_only),
      "",
      communication_report_md_table(pattern_counts, max_rows = 30L),
      sprintf("![method agreement](%s)", plot_path)
    ),
    plots = list(method_agreement = plot_path),
    tsv = pattern_counts
  )
}

render_panel_downstream_chain <- function(consensus_df, cfg) {
  df <- communication_report_normalize_consensus(consensus_df)
  targets <- suppressWarnings(as.numeric(df$nichenet_n_targets_in_receiver_de %||% NA_real_))
  downstream <- ifelse(df$nichenet_hit | (!is.na(targets) & targets > 0), "yes", "no")
  chain <- data.frame(
    pair_id = df$pair_id %||% character(0),
    condition_value = df$condition_value %||% character(0),
    downstream_support = downstream,
    n_axes = 1L,
    stringsAsFactors = FALSE
  )
  if (nrow(chain) > 0) {
    chain <- aggregate(n_axes ~ pair_id + condition_value + downstream_support, data = chain, FUN = sum)
  }
  plot_path <- communication_report_plot_path(cfg, "panel_03_downstream_chain.png")
  values <- if (nrow(chain) == 0) numeric(0) else tapply(chain$n_axes, chain$downstream_support, sum)
  communication_report_barplot_png(values, plot_path, "Downstream target support", ylab = "Axis count")
  list(
    name = "downstream_target_chain",
    status = "ok",
    md_lines = c(
      "## Panel 3: Downstream Target Chain",
      "",
      "MultiNicheNet/NicheNet downstream support is summarized as ligand-target evidence intersecting receiver DE when available.",
      "",
      communication_report_md_table(chain, max_rows = 60L),
      sprintf("![downstream chain](%s)", plot_path)
    ),
    plots = list(downstream_chain = plot_path),
    tsv = chain
  )
}

render_panel_legacy_scdesign3 <- function(cfg) {
  candidates <- c(
    file.path(cfg$table_dir, "04f_scdesign3_finalize", "communication_scdesign3_gate_post_engine.tsv"),
    file.path(cfg$table_dir, "04d_cluster_robustness", "communication_scdesign3_gate.tsv")
  )
  existing <- candidates[file.exists(candidates)]
  status <- if (length(existing) > 0) "available" else "not_run"
  tsv <- data.frame(status = status, path = paste(existing, collapse = ","), stringsAsFactors = FALSE)
  list(
    name = "legacy_scdesign3_gate",
    status = status,
    md_lines = c(
      "## Legacy: scDesign3 Gate",
      "",
      sprintf("- Status: %s", status),
      if (length(existing) > 0) sprintf("- Source: `%s`", existing[[1]]) else "- Source: not available"
    ),
    plots = list(),
    tsv = tsv
  )
}

render_panel_legacy_fallback <- function(cfg) {
  path <- cfg$communication_fallback_summary_tsv %||% ""
  df <- if (exists("read_tsv_optional", mode = "function")) read_tsv_optional(path) else data.frame(stringsAsFactors = FALSE)
  tsv <- if (nrow(df) == 0) data.frame(status = "not_available", path = path, stringsAsFactors = FALSE) else df
  list(
    name = "legacy_fallback_summary",
    status = if (nrow(df) == 0) "not_available" else "available",
    md_lines = c(
      "## Legacy: Fallback Summary",
      "",
      communication_report_md_table(tsv, max_rows = 40L)
    ),
    plots = list(),
    tsv = tsv
  )
}

render_panel_legacy_missing_gene_program <- function(cfg) {
  path <- cfg$nichenet_index_tsv %||% ""
  df <- if (exists("read_tsv_optional", mode = "function")) read_tsv_optional(path) else data.frame(stringsAsFactors = FALSE)
  if (nrow(df) > 0 && "status" %in% colnames(df)) {
    df <- df[df$status == "missing_gene_program", , drop = FALSE]
  } else {
    df <- data.frame(stringsAsFactors = FALSE)
  }
  tsv <- if (nrow(df) == 0) data.frame(status = "none", stringsAsFactors = FALSE) else df
  list(
    name = "legacy_missing_gene_program",
    status = if (nrow(df) == 0) "none" else "available",
    md_lines = c(
      "## Legacy: Missing Gene Programs",
      "",
      communication_report_md_table(tsv, max_rows = 40L)
    ),
    plots = list(),
    tsv = tsv
  )
}
