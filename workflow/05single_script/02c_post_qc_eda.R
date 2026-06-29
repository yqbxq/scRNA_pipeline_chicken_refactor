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
source_utf8(file.path(.script_dir, "helpers", "gtf_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "plotting_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tidyr", "tibble", "scales", "jsonlite"))

cfg <- get_single_script_config_02()
module_name <- "02c_post_qc_eda"
prepare_dirs_02(cfg)
set.seed(cfg$random_seed)

manifest_01a <- read_manifest_local(cfg$module_01a_manifest_path)
manifest_01b <- read_manifest_local(cfg$module_01b_manifest_path)
manifest_02a <- read_manifest_local(cfg$module_02a_manifest_path)
manifest_02b1 <- read_manifest_local(cfg$module_02b1_manifest_path)
manifest_02b2 <- read_manifest_local(cfg$module_02b2_manifest_path)

raw_rds <- resolve_output_local(manifest_01a, "raw_object")
pre_qc_summary_path <- resolve_output_local(manifest_01b, "sample_qc_summary")
ambient_summary_path <- resolve_output_local(manifest_02a, "ambient_sample_summary")
qc_filter_summary_path <- resolve_output_local(manifest_02b1, "qc_filter_sample_summary")
qc_filter_diagnostics_path <- resolve_output_local(manifest_02b1, "qc_filter_diagnostics")
post_qc_rds <- resolve_output_local(manifest_02b2, "post_qc_object")
doublet_summary_path <- resolve_output_local(manifest_02b2, "doublet_sample_summary")
doublet_cell_path <- resolve_output_local(manifest_02b2, "doublet_cell_diagnostics")
concordance_path <- resolve_output_local(manifest_02b2, "doublet_concordance_summary")
cluster_risk_path <- resolve_output_local(manifest_02b2, "doublet_cluster_risk")

for (required_path in c(
  raw_rds, pre_qc_summary_path, ambient_summary_path, qc_filter_summary_path,
  qc_filter_diagnostics_path, post_qc_rds, doublet_summary_path, doublet_cell_path,
  concordance_path, cluster_risk_path
)) {
  if (!file.exists(required_path)) {
    stop(sprintf("缺少 post-QC EDA 输入: %s", required_path), call. = FALSE)
  }
}

stage_dir <- cfg$post_qc_report_dir
raw_obj <- maybe_join_layers(readRDS(raw_rds))
post_qc_obj <- maybe_join_layers(readRDS(post_qc_rds))
pre_qc_summary <- read.delim(pre_qc_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ambient_summary <- read.delim(ambient_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
qc_filter_summary <- read.delim(qc_filter_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
qc_filter_diagnostics <- read.csv(qc_filter_diagnostics_path, stringsAsFactors = FALSE, check.names = FALSE)
doublet_summary <- read.delim(doublet_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
doublet_cells <- read.csv(doublet_cell_path, stringsAsFactors = FALSE, check.names = FALSE)
concordance_df <- read.delim(concordance_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cluster_risk_df <- read.delim(cluster_risk_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

raw_counts <- raw_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::mutate(sample_id = if ("sample_id" %in% colnames(.)) sample_id else orig.ident) %>%
  dplyr::count(sample_id, name = "raw_cells")

pre_counts <- if (all(c("sample_id", "cells_raw") %in% colnames(pre_qc_summary))) {
  pre_qc_summary[, c("sample_id", "cells_raw"), drop = FALSE] %>%
    dplyr::rename(cells_pre_qc = cells_raw)
} else {
  raw_counts %>% dplyr::rename(cells_pre_qc = raw_cells)
}

post_counts <- doublet_cells %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(cells_post_qc = sum(final_call == "singlet", na.rm = TRUE), .groups = "drop")

post_qc_sample_summary <- doublet_summary %>%
  dplyr::left_join(pre_counts, by = "sample_id") %>%
  dplyr::left_join(post_counts, by = "sample_id") %>%
  dplyr::mutate(
    cells_pre_qc = ifelse(is.na(cells_pre_qc), 0, cells_pre_qc),
    cells_post_qc = ifelse(is.na(cells_post_qc), 0, cells_post_qc),
    retention_rate = ifelse(cells_pre_qc > 0, cells_post_qc / cells_pre_qc, NA_real_)
  )

qc_retention_comparison <- raw_counts %>%
  dplyr::left_join(
    ambient_summary %>%
      dplyr::transmute(sample_id, ambient_corrected = cells_input),
    by = "sample_id"
  ) %>%
  dplyr::left_join(
    qc_filter_summary %>%
      dplyr::transmute(sample_id, qc_filtered = cells_pass),
    by = "sample_id"
  ) %>%
  dplyr::left_join(
    post_counts %>%
      dplyr::transmute(sample_id, singlet = cells_post_qc),
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    ambient_corrected = ifelse(is.na(ambient_corrected), raw_cells, ambient_corrected),
    qc_filtered = ifelse(is.na(qc_filtered), 0, qc_filtered),
    singlet = ifelse(is.na(singlet), 0, singlet)
  )

triage_rows <- list()
append_triage <- function(sample_id, severity, signal_id, evidence) {
  triage_rows[[length(triage_rows) + 1]] <<- make_triage_row(
    sample_id = sample_id,
    severity = severity,
    signal_id = signal_id,
    evidence = evidence
  )
}

as_boolish <- function(x) {
  tolower(normalize_scalar_value(x, "false")) %in% c("true", "yes", "1", "on")
}

for (i in seq_len(nrow(post_qc_sample_summary))) {
  row <- post_qc_sample_summary[i, , drop = FALSE]
  sample_id <- row$sample_id

  if (!is.na(row$retention_rate) && row$retention_rate < 0.30) {
    append_triage(sample_id, "high", "heavy_qc_loss", sprintf("retention_rate=%.3f; cells_pre_qc=%s; cells_post_qc=%s", row$retention_rate, row$cells_pre_qc, row$cells_post_qc))
  }
  if (!is.na(row$observed_primary_rate) && row$observed_primary_rate > 0.15) {
    append_triage(sample_id, "high", "high_primary_doublet_rate", sprintf("observed_primary_rate=%.3f; status=%s", row$observed_primary_rate, row$status))
  }
  discordant <- suppressWarnings(as.numeric(row$primary_only)) + suppressWarnings(as.numeric(row$secondary_only))
  if (is.finite(discordant) && discordant > 0) {
    append_triage(sample_id, "medium", "primary_secondary_discordance", sprintf("primary_only=%s; secondary_only=%s; concordant=%s", row$primary_only, row$secondary_only, row$concordant))
  }
  if (!is.na(row$cells_post_qc) && row$cells_post_qc < 100) {
    append_triage(sample_id, "high", "sample_collapse", sprintf("cells_post_qc=%s", row$cells_post_qc))
  }
}

if (nrow(ambient_summary) > 0 && all(c("sample_id", "recommended_to_replace", "applied_to_main") %in% colnames(ambient_summary))) {
  for (i in seq_len(nrow(ambient_summary))) {
    row <- ambient_summary[i, , drop = FALSE]
    if (as_boolish(row$recommended_to_replace) && !as_boolish(row$applied_to_main)) {
      append_triage(
        row$sample_id,
        "info",
        "ambient_recommended_not_applied",
        sprintf("executed_method=%s; contamination_fraction=%s", normalize_scalar_value(row$executed_method, "none"), normalize_scalar_value(row$contamination_fraction, "NA"))
      )
    }
  }
}

if (nrow(cluster_risk_df) > 0 && "risk_level" %in% colnames(cluster_risk_df)) {
  risky <- cluster_risk_df[cluster_risk_df$risk_level %in% c("high"), , drop = FALSE]
  for (i in seq_len(nrow(risky))) {
    row <- risky[i, , drop = FALSE]
    append_triage(
      row$sample_id,
      "medium",
      "doublet_high_risk_cluster",
      sprintf("cluster=%s; doublet_fraction=%.3f; top_marker=%s; risk_level=%s", row$cluster, suppressWarnings(as.numeric(row$doublet_fraction)), normalize_scalar_value(row$top_marker, "NA"), row$risk_level)
    )
  }
}

triage_df <- if (length(triage_rows) > 0) {
  dplyr::bind_rows(triage_rows)
} else {
  empty_triage_df()
}

metric_long <- function(seu, stage_name) {
  meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      sample_id = if ("sample_id" %in% colnames(.)) sample_id else orig.ident,
      percent.mito = if ("percent.mito" %in% colnames(.)) percent.mito else 0,
      percent.rbc = if ("percent.rbc" %in% colnames(.)) percent.rbc else 0,
      log10GenesPerUMI = if ("log10GenesPerUMI" %in% colnames(.)) log10GenesPerUMI else log10(nFeature_RNA + 1) / log10(nCount_RNA + 1)
    )
  meta$log10GenesPerUMI[!is.finite(meta$log10GenesPerUMI)] <- 0
  meta %>%
    dplyr::select(sample_id, nCount_RNA, nFeature_RNA, percent.mito, percent.rbc, log10GenesPerUMI) %>%
    tidyr::pivot_longer(
      cols = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.rbc", "log10GenesPerUMI"),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::mutate(stage = stage_name)
}

pre_post_long <- dplyr::bind_rows(
  metric_long(raw_obj, "pre_qc"),
  metric_long(post_qc_obj, "post_qc")
)

violin_plot <- ggplot2::ggplot(pre_post_long, ggplot2::aes(x = stage, y = value, fill = stage)) +
  ggplot2::geom_violin(scale = "width", trim = TRUE) +
  ggplot2::facet_grid(metric ~ sample_id, scales = "free_y") +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none") +
  ggplot2::labs(title = "Pre vs post-QC metric distributions", x = NULL, y = NULL)

score_df <- doublet_cells %>%
  dplyr::mutate(primary_score = suppressWarnings(as.numeric(primary_score))) %>%
  dplyr::filter(is.finite(primary_score))
if (nrow(score_df) > 0) {
  score_plot <- ggplot2::ggplot(score_df, ggplot2::aes(x = primary_score, fill = final_call)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.75, position = "identity") +
    ggplot2::facet_wrap(~ sample_id, scales = "free_y") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::labs(title = "Doublet score distribution", x = "primary_score", y = "Cells", fill = "Final call")
} else {
  score_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No primary doublet scores available") +
    ggplot2::theme_void()
}

if (nrow(cluster_risk_df) > 0) {
  heatmap_plot <- ggplot2::ggplot(cluster_risk_df, ggplot2::aes(x = cluster, y = sample_id, fill = doublet_fraction)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#CB181D", labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Cluster doublet heatmap", x = "Cluster", y = NULL, fill = "Doublet fraction")
} else {
  heatmap_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No cluster doublet risk table available") +
    ggplot2::theme_void()
}

post_qc_summary_path <- file.path(stage_dir, "post_qc_sample_summary.tsv")
post_qc_triage_path <- file.path(stage_dir, "post_qc_triage.tsv")
retention_path <- file.path(stage_dir, "qc_retention_comparison.tsv")
report_path <- file.path(stage_dir, "report.md")
violin_png <- file.path(stage_dir, "pre_vs_post_qc_violin.png")
score_png <- file.path(stage_dir, "doublet_score_distribution.png")
heatmap_png <- file.path(stage_dir, "cluster_doublet_heatmap.png")

write_tsv_local(post_qc_sample_summary, post_qc_summary_path)
write_tsv_local(triage_df, post_qc_triage_path)
write_tsv_local(qc_retention_comparison, retention_path)
save_plot_dual(violin_plot, violin_png, width = 12, height = 9)
save_plot_dual(score_plot, score_png, width = 10, height = 6)
save_plot_dual(heatmap_plot, heatmap_png, width = 10, height = 6)

report_lines <- c(
  "# Post-QC EDA Report",
  "",
  sprintf("- 样本数: `%s`", nrow(post_qc_sample_summary)),
  sprintf("- QC/doublet 后总细胞数: `%s`", format(sum(post_qc_sample_summary$cells_post_qc, na.rm = TRUE), big.mark = ",")),
  sprintf("- triage 信号条数: `%s`", nrow(triage_df)),
  "",
  "## Key Files",
  sprintf("- `post_qc_sample_summary.tsv`: `%s`", post_qc_summary_path),
  sprintf("- `post_qc_triage.tsv`: `%s`", post_qc_triage_path),
  sprintf("- `qc_retention_comparison.tsv`: `%s`", retention_path),
  sprintf("- `pre_vs_post_qc_violin.png`: `%s`", violin_png),
  sprintf("- `doublet_score_distribution.png`: `%s`", score_png),
  sprintf("- `cluster_doublet_heatmap.png`: `%s`", heatmap_png),
  "",
  "## Review Focus",
  "- scDblFinder 是 primary caller；DoubletFinder 仅作为 optional secondary/concordance 证据。",
  "- 重点看 retention、primary/secondary discordance、cluster-level doublet risk，以及 ambient correction 是否已应用。"
)

report_lines <- c(
  report_lines,
  "",
  "## Triage Summary",
  render_triage_markdown(triage_df),
  "",
  "## Next Step",
  sprintf("审阅以上报告后，在 `%s` 中将 `post_qc` 的 `status` 改为 `approved`，然后继续 integration 模块。", cfg$eda_gate_file)
)

ensure_dir(dirname(report_path))
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_02c_manifest_path,
  new_outputs = list(
    post_qc_sample_summary = build_output_entry(post_qc_summary_path, "tsv", module_name, "one row per sample post-QC summary", base_dir = cfg$project_root, schema = infer_schema_from_df(post_qc_sample_summary)),
    post_qc_triage = build_output_entry(post_qc_triage_path, "tsv", module_name, "one row per post-QC triage signal", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    qc_retention_comparison = build_output_entry(retention_path, "tsv", module_name, "one row per sample cell count by QC stage", base_dir = cfg$project_root, schema = infer_schema_from_df(qc_retention_comparison)),
    pre_vs_post_qc_violin_png = build_output_entry(violin_png, "png", module_name, "pre-vs-post QC violin plot", base_dir = cfg$project_root),
    pre_vs_post_qc_violin_pdf = build_output_entry(sub("\\.png$", ".pdf", violin_png), "pdf", module_name, "pre-vs-post QC violin plot", base_dir = cfg$project_root),
    doublet_score_distribution_png = build_output_entry(score_png, "png", module_name, "primary doublet score distribution", base_dir = cfg$project_root),
    doublet_score_distribution_pdf = build_output_entry(sub("\\.png$", ".pdf", score_png), "pdf", module_name, "primary doublet score distribution", base_dir = cfg$project_root),
    cluster_doublet_heatmap_png = build_output_entry(heatmap_png, "png", module_name, "cluster doublet risk heatmap", base_dir = cfg$project_root),
    cluster_doublet_heatmap_pdf = build_output_entry(sub("\\.png$", ".pdf", heatmap_png), "pdf", module_name, "cluster doublet risk heatmap", base_dir = cfg$project_root),
    post_qc_report = build_output_entry(report_path, "md", module_name, "human-readable post-QC EDA report", base_dir = cfg$project_root),
    report = build_output_entry(report_path, "md", module_name, "human-readable post-QC EDA report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    manifest_01a = cfg$module_01a_manifest_path,
    manifest_01b = cfg$module_01b_manifest_path,
    manifest_02a = cfg$module_02a_manifest_path,
    manifest_02b1 = cfg$module_02b1_manifest_path,
    manifest_02b2 = cfg$module_02b2_manifest_path,
    raw_object = raw_rds,
    pre_qc_summary = pre_qc_summary_path,
    ambient_sample_summary = ambient_summary_path,
    qc_filter_sample_summary = qc_filter_summary_path,
    qc_filter_diagnostics = qc_filter_diagnostics_path,
    post_qc_object = post_qc_rds,
    doublet_sample_summary = doublet_summary_path,
    doublet_cell_diagnostics = doublet_cell_path,
    doublet_concordance_summary = concordance_path,
    doublet_cluster_risk = cluster_risk_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_01a = cfg$module_01a_manifest_path,
    module_01b = cfg$module_01b_manifest_path,
    module_02a = cfg$module_02a_manifest_path,
    module_02b1 = cfg$module_02b1_manifest_path,
    module_02b2 = cfg$module_02b2_manifest_path
  )
)

message("02c 完成。post-QC EDA 已输出到: ", stage_dir)
