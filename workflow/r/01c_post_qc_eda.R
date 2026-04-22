source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

raw_rds <- file.path(cfg$checkpoint_dir, "00_raw_objects.rds")
qc_rds <- file.path(cfg$checkpoint_dir, "01_after_qc_doublet.rds")
summary_csv <- file.path(cfg$table_dir, "qc_doublet_summary.csv")
diag_csv <- file.path(cfg$table_dir, "qc_doublet_cell_level.csv.gz")
doublet_summary_tsv <- file.path(cfg$table_dir, "doublet_sample_summary.tsv")
concordance_tsv <- file.path(cfg$table_dir, "doublet_concordance_summary.tsv")
cluster_risk_tsv <- file.path(cfg$table_dir, "doublet_cluster_risk.tsv")
ambient_summary_tsv <- file.path(cfg$eda_report_dir, "ambient", "ambient_sample_summary.tsv")

for (required_path in c(raw_rds, qc_rds, summary_csv, diag_csv, doublet_summary_tsv, concordance_tsv, cluster_risk_tsv)) {
  if (!file.exists(required_path)) {
    stop(sprintf("缺少 post-QC EDA 输入: %s", required_path), call. = FALSE)
  }
}

stage_dir <- ensure_eda_stage_dir(cfg, "post_qc")
raw_obj <- readRDS(raw_rds)
raw_obj <- maybe_join_layers(raw_obj)
final_obj <- readRDS(qc_rds)
final_obj <- maybe_join_layers(final_obj)
summary_df <- read.csv(summary_csv, stringsAsFactors = FALSE, check.names = FALSE)
diagnostic_df <- read.csv(gzfile(diag_csv), stringsAsFactors = FALSE, check.names = FALSE)
doublet_summary_df <- read.delim(doublet_summary_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
concordance_df <- read.delim(concordance_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
cluster_risk_df <- read.delim(cluster_risk_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
ambient_summary_df <- if (file.exists(ambient_summary_tsv)) {
  read.delim(ambient_summary_tsv, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(stringsAsFactors = FALSE)
}

retention_df <- summary_df %>%
  mutate(
    qc_retention = ifelse(cells_raw > 0, cells_after_qc / cells_raw, NA_real_),
    singlet_retention = ifelse(cells_after_qc > 0, singlets_after_doublet / cells_after_qc, NA_real_),
    final_retention = ifelse(cells_raw > 0, singlets_after_doublet / cells_raw, NA_real_)
  )

if (nrow(ambient_summary_df) > 0 && "sample_id" %in% colnames(ambient_summary_df)) {
  retention_df <- retention_df %>%
    left_join(
      ambient_summary_df %>%
        select(sample_id, executed_method, contamination_fraction, recommended_to_replace, applied_to_main),
      by = c("sample" = "sample_id")
    )
}

triage_rows <- list()
append_triage <- function(sample_id, severity, signal_id, evidence, recommended_action) {
  triage_rows[[length(triage_rows) + 1]] <<- make_triage_row(
    sample_id = sample_id,
    severity = severity,
    signal_id = signal_id,
    evidence = evidence,
    recommended_action = recommended_action
  )
}

for (i in seq_len(nrow(retention_df))) {
  row <- retention_df[i, , drop = FALSE]
  sample_id <- row$sample

  if (!is.na(row$qc_retention) && row$qc_retention < 0.30) {
    append_triage(
      sample_id,
      "high",
      "heavy_qc_loss",
      sprintf("qc_retention=%.3f", row$qc_retention),
      "该样本在 QC 后保留率偏低，建议复核 pre-QC 阈值和样本级质量。"
    )
  }

  if (!is.na(row$primary_detected_rate) && row$primary_detected_rate > 0.15) {
    append_triage(
      sample_id,
      "high",
      "high_primary_doublet_rate",
      sprintf("primary_detected_rate=%.3f; primary_status=%s", row$primary_detected_rate, row$primary_status),
      "primary caller 检出的 doublet rate 偏高，建议重点审阅 doublet UMAP、score 分布与高风险 cluster。"
    )
  }

  if (!is.na(row$discordant_calls) && row$discordant_calls > 0) {
    append_triage(
      sample_id,
      "medium",
      "primary_secondary_discordance",
      sprintf("discordant_calls=%s; final_filter_caller=%s", row$discordant_calls, row$final_filter_caller),
      "存在 primary/secondary 冲突，建议在保留 primary 为默认主调用器的前提下复核冲突区域。"
    )
  }

  if (row$singlets_after_doublet < 100) {
    append_triage(
      sample_id,
      "high",
      "sample_collapse",
      sprintf("singlets_after_doublet=%s", row$singlets_after_doublet),
      "该样本在 QC/doublet 后剩余细胞过少，后续聚类与整合解释风险较高。"
    )
  }

  if (!is.na(row$contamination_fraction) && row$contamination_fraction >= cfg$ambient_recommend_min_contamination && !isTRUE(as.logical(row$applied_to_main))) {
    append_triage(
      sample_id,
      "info",
      "ambient_recommended_not_applied",
      sprintf("ambient_method=%s; contamination_fraction=%.3f", normalize_scalar_value(row$executed_method, "none"), row$contamination_fraction),
      "ambient 分支建议替换 counts，但当前主对象仍保留原始 counts；进入后续解释前应明确这一点。"
    )
  }
}

if (nrow(cluster_risk_df) > 0) {
  high_risk_rows <- cluster_risk_df[cluster_risk_df$high_risk_cluster == "yes", , drop = FALSE]
  if (nrow(high_risk_rows) > 0) {
    for (i in seq_len(nrow(high_risk_rows))) {
      row <- high_risk_rows[i, , drop = FALSE]
      append_triage(
        row$sample_id,
        "medium",
        "doublet_high_risk_cluster",
        sprintf(
          "cluster=%s; primary_doublet_rate=%.3f; enrichment_vs_sample_baseline=%.3f",
          row$provisional_cluster,
          row$primary_doublet_rate,
          row$enrichment_vs_sample_baseline
        ),
        "该 cluster/区域存在 doublet 富集，建议结合 UMAP 空间与 marker 混合模式复核。"
      )
    }
  }
}

triage_df <- if (length(triage_rows) > 0) bind_rows(triage_rows) else empty_triage_df()

retention_long <- retention_df %>%
  select(sample, qc_retention, singlet_retention, final_retention) %>%
  pivot_longer(cols = -sample, names_to = "stage", values_to = "fraction")

retention_plot <- ggplot(retention_long, aes(x = sample, y = fraction, fill = stage)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Retention after QC and doublet removal", x = NULL, y = "Fraction")

doublet_plot_df <- diagnostic_df %>%
  filter(qc_pass, !is.na(provisional_umap_1), !is.na(provisional_umap_2)) %>%
  mutate(primary_class = ifelse(tolower(scDblFinder.class) == "doublet", "Doublet", "Singlet"))

doublet_plot <- ggplot(doublet_plot_df, aes(x = provisional_umap_1, y = provisional_umap_2, color = primary_class)) +
  geom_point(size = 0.3, alpha = 0.7) +
  facet_wrap(~ sample_id, scales = "free") +
  scale_color_manual(values = c("Singlet" = "#4C78A8", "Doublet" = "#D62728")) +
  theme_classic(base_size = 11) +
  labs(title = "scDblFinder calls in sample-wise UMAP space", x = "UMAP_1", y = "UMAP_2", color = NULL)

score_plot_df <- diagnostic_df %>%
  filter(qc_pass, is.finite(scDblFinder.score))
score_plot <- ggplot(score_plot_df, aes(x = scDblFinder.score, fill = sample_id)) +
  geom_histogram(bins = 30, alpha = 0.75, position = "identity") +
  facet_wrap(~ sample_id, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(title = "scDblFinder score review", x = "scDblFinder.score", y = "Cells", fill = NULL)

cluster_risk_plot_df <- cluster_risk_df %>%
  mutate(provisional_cluster = factor(provisional_cluster, levels = unique(provisional_cluster)))
cluster_risk_plot <- ggplot(cluster_risk_plot_df, aes(x = provisional_cluster, y = primary_doublet_rate, fill = high_risk_cluster)) +
  geom_col() +
  facet_wrap(~ sample_id, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Primary doublet enrichment by provisional cluster", x = "Provisional cluster", y = "Primary doublet rate", fill = "High risk")

final_meta <- final_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(sample_id = if ("sample_id" %in% colnames(.)) sample_id else orig.ident) %>%
  select(sample_id, nCount_RNA, nFeature_RNA, percent.mito, log10GenesPerUMI) %>%
  pivot_longer(
    cols = c("nCount_RNA", "nFeature_RNA", "percent.mito", "log10GenesPerUMI"),
    names_to = "metric",
    values_to = "value"
  )

final_violin <- ggplot(final_meta, aes(x = sample_id, y = value, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Post-QC metric distributions", x = NULL, y = NULL)

save_plot_dual(retention_plot, file.path(stage_dir, "retention_summary.png"), width = 10, height = 5)
save_plot_dual(doublet_plot, file.path(stage_dir, "doublet_review.png"), width = 12, height = 8)
save_plot_dual(score_plot, file.path(stage_dir, "doublet_score_review.png"), width = 12, height = 8)
save_plot_dual(cluster_risk_plot, file.path(stage_dir, "doublet_cluster_risk.png"), width = 12, height = 8)
save_plot_dual(final_violin, file.path(stage_dir, "post_qc_violin.png"), width = 12, height = 8)

write_tsv(retention_df, file.path(stage_dir, "retention_summary.tsv"))
write_tsv(diagnostic_df, file.path(stage_dir, "qc_before_after.tsv"))
write_tsv(triage_df, file.path(stage_dir, "post_qc_triage.tsv"))
write_tsv(doublet_summary_df, file.path(stage_dir, "doublet_sample_summary.tsv"))
write_tsv(concordance_df, file.path(stage_dir, "doublet_concordance_summary.tsv"))
write_tsv(cluster_risk_df, file.path(stage_dir, "doublet_cluster_risk.tsv"))

report_lines <- c(
  "# Post-QC EDA Report",
  "",
  sprintf("- 样本数: `%s`", nrow(retention_df)),
  sprintf("- QC/doublet 后总细胞数: `%s`", format(sum(retention_df$singlets_after_doublet), big.mark = ",")),
  sprintf("- triage 信号条数: `%s`", nrow(triage_df)),
  "",
  "## Key Files",
  sprintf("- `retention_summary.tsv`: `%s`", file.path(stage_dir, "retention_summary.tsv")),
  sprintf("- `doublet_sample_summary.tsv`: `%s`", file.path(stage_dir, "doublet_sample_summary.tsv")),
  sprintf("- `doublet_concordance_summary.tsv`: `%s`", file.path(stage_dir, "doublet_concordance_summary.tsv")),
  sprintf("- `doublet_cluster_risk.tsv`: `%s`", file.path(stage_dir, "doublet_cluster_risk.tsv")),
  sprintf("- `qc_before_after.tsv`: `%s`", file.path(stage_dir, "qc_before_after.tsv")),
  sprintf("- `post_qc_triage.tsv`: `%s`", file.path(stage_dir, "post_qc_triage.tsv")),
  sprintf("- `retention_summary.png`: `%s`", file.path(stage_dir, "retention_summary.png")),
  sprintf("- `doublet_review.png`: `%s`", file.path(stage_dir, "doublet_review.png")),
  sprintf("- `doublet_score_review.png`: `%s`", file.path(stage_dir, "doublet_score_review.png")),
  sprintf("- `doublet_cluster_risk.png`: `%s`", file.path(stage_dir, "doublet_cluster_risk.png")),
  "",
  "## Review Focus",
  "- `scDblFinder` 已是默认 primary caller；`DoubletFinder` 仅作为 secondary/concordance 证据，不再是唯一默认主线。",
  "- 重点看每个样本的 primary detected rate、primary/secondary discordance，以及 provisional cluster 层面的 doublet 富集。",
  "- 若 ambient 报告建议替换 counts 但当前仍未应用，应在进入后续整合和注释前显式记录。"
)

report_lines <- c(report_lines, "", "## Triage Summary", render_triage_markdown(triage_df))

write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("post-QC EDA 已输出到: ", stage_dir)
