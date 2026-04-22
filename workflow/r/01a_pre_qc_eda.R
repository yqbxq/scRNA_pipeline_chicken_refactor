source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "00_raw_objects.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 00_raw_objects.rds", call. = FALSE)
}

stage_dir <- ensure_eda_stage_dir(cfg, "pre_qc")
raw_obj <- readRDS(input_rds)
raw_obj <- maybe_join_layers(raw_obj)

sample_sheet <- main_sample_sheet(cfg)
inventory_df <- read_input_inventory(cfg)
branch_readiness_df <- read_branch_readiness(cfg)
threshold_overrides <- read_qc_threshold_overrides(cfg)

density_peaks <- function(x) {
  x <- stats::na.omit(as.numeric(x))
  if (length(x) < 100 || length(unique(x)) < 20) {
    return(1L)
  }
  dens <- stats::density(x, na.rm = TRUE)
  y <- dens$y
  max(1L, sum(diff(sign(diff(y))) < 0))
}

meta_df <- raw_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    sample_id = if ("sample_id" %in% colnames(.)) sample_id else orig.ident,
    analysis_group = if ("analysis_group" %in% colnames(.)) analysis_group else orig.ident,
    percent.ribo = if ("percent.ribo" %in% colnames(.)) percent.ribo else 0,
    mito_feature_count = if ("mito_feature_count" %in% colnames(.)) mito_feature_count else NA_real_,
    ribo_feature_count = if ("ribo_feature_count" %in% colnames(.)) ribo_feature_count else NA_real_,
    cell_cycle_s_feature_count = if ("cell_cycle_s_feature_count" %in% colnames(.)) cell_cycle_s_feature_count else NA_real_,
    cell_cycle_g2m_feature_count = if ("cell_cycle_g2m_feature_count" %in% colnames(.)) cell_cycle_g2m_feature_count else NA_real_
  )

sample_feature_contracts <- raw_obj@misc$sample_feature_contracts
if (is.null(sample_feature_contracts) || !is.data.frame(sample_feature_contracts)) {
  sample_feature_contracts <- data.frame(sample_id = character(0), stringsAsFactors = FALSE)
}

sample_summary <- bind_rows(lapply(unique(meta_df$sample_id), function(sample_id) {
  sample_meta <- meta_df[meta_df$sample_id == sample_id, , drop = FALSE]
  inventory_row <- resolve_inventory_row(cfg, sample_id, inventory_df = inventory_df)
  readiness_row <- if (nrow(branch_readiness_df) > 0 && "sample_id" %in% colnames(branch_readiness_df)) {
    branch_readiness_df[branch_readiness_df$sample_id == sample_id, , drop = FALSE]
  } else {
    NULL
  }
  thresholds <- get_sample_qc_thresholds(cfg, sample_id, threshold_overrides)
  metrics_path <- resolve_metrics_summary_path(cfg, sample_id, inventory_df = inventory_df)
  metrics <- read_cellranger_metrics(metrics_path)

  data.frame(
    sample_id = sample_id,
    platform = if (!is.null(inventory_row) && "platform_resolved" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$platform_resolved[1]) else "",
    gene_id_type_resolved = if (!is.null(inventory_row) && "gene_id_type_resolved" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$gene_id_type_resolved[1], "unknown") else "unknown",
    feature_name_profile = if (!is.null(inventory_row) && "feature_name_profile" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$feature_name_profile[1], "unknown") else "unknown",
    reference_version = if (!is.null(inventory_row) && "reference_version" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$reference_version[1]) else "",
    filtered_matrix_dir = if (!is.null(inventory_row) && "filtered_matrix_dir" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$filtered_matrix_dir[1]) else "",
    raw_matrix_dir = if (!is.null(inventory_row) && "raw_matrix_dir" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$raw_matrix_dir[1]) else "",
    raw_matrix_kind = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "raw_matrix_kind" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$raw_matrix_kind[1]) else "",
    bam_path = if (!is.null(inventory_row) && "bam_path" %in% colnames(inventory_row)) normalize_scalar_value(inventory_row$bam_path[1]) else "",
    cells_raw = nrow(sample_meta),
    median_ncount = stats::median(sample_meta$nCount_RNA),
    median_nfeature = stats::median(sample_meta$nFeature_RNA),
    median_percent_mito = stats::median(sample_meta$percent.mito),
    median_percent_ribo = stats::median(sample_meta$percent.ribo),
    median_log10GenesPerUMI = stats::median(sample_meta$log10GenesPerUMI),
    mito_feature_count = if ("mito_feature_count" %in% colnames(sample_meta)) dplyr::first(sample_meta$mito_feature_count) else NA_real_,
    ribo_feature_count = if ("ribo_feature_count" %in% colnames(sample_meta)) dplyr::first(sample_meta$ribo_feature_count) else NA_real_,
    cell_cycle_s_feature_count = if ("cell_cycle_s_feature_count" %in% colnames(sample_meta)) dplyr::first(sample_meta$cell_cycle_s_feature_count) else NA_real_,
    cell_cycle_g2m_feature_count = if ("cell_cycle_g2m_feature_count" %in% colnames(sample_meta)) dplyr::first(sample_meta$cell_cycle_g2m_feature_count) else NA_real_,
    frac_below_nfeature_cutoff = mean(sample_meta$nFeature_RNA < thresholds$qc_min_nfeature),
    frac_below_ncount_cutoff = mean(sample_meta$nCount_RNA < thresholds$qc_min_ncount),
    frac_below_log10umi_cutoff = mean(sample_meta$log10GenesPerUMI < thresholds$qc_min_log10umi),
    frac_above_mito_cutoff = mean(sample_meta$percent.mito > thresholds$qc_max_mito_pct),
    nfeature_density_peaks = density_peaks(sample_meta$nFeature_RNA),
    mito_density_peaks = density_peaks(sample_meta$percent.mito),
    qc_min_nfeature = thresholds$qc_min_nfeature,
    qc_min_ncount = thresholds$qc_min_ncount,
    qc_min_log10umi = thresholds$qc_min_log10umi,
    qc_max_mito_pct = thresholds$qc_max_mito_pct,
    estimated_cells = if (!is.null(metrics)) metrics$estimated_cells else NA_real_,
    mean_reads_per_cell = if (!is.null(metrics)) metrics$mean_reads_per_cell else NA_real_,
    median_genes_per_cell = if (!is.null(metrics)) metrics$median_genes_per_cell else NA_real_,
    sequencing_saturation = if (!is.null(metrics)) metrics$sequencing_saturation else NA_real_,
    fraction_reads_in_cells = if (!is.null(metrics)) metrics$fraction_reads_in_cells else NA_real_,
    metrics_summary_path = metrics_path,
    ambient_any_ready = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_any_ready" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_any_ready[1], "false") else "false",
    ambient_soupx_ready = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_soupx_ready" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_soupx_ready[1], "false") else "false",
    ambient_decontx_ready = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_decontx_ready" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_decontx_ready[1], "false") else "false",
    ambient_cellbender_ready = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_cellbender_ready" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_cellbender_ready[1], "false") else "false",
    ambient_preferred_method = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_preferred_method" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_preferred_method[1], "none") else "none",
    ambient_fallback_method = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_fallback_method" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_fallback_method[1], "none") else "none",
    ambient_notes = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "ambient_notes" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$ambient_notes[1]) else "",
    readiness_notes = if (!is.null(readiness_row) && nrow(readiness_row) > 0 && "notes" %in% colnames(readiness_row)) normalize_scalar_value(readiness_row$notes[1]) else "",
    stringsAsFactors = FALSE
  )
})) %>%
  left_join(sample_feature_contracts, by = "sample_id")

triage_rows <- list()
append_triage <- function(sample_id, severity, signal_id, suspected_issue, evidence, recommended_action) {
  triage_rows[[length(triage_rows) + 1]] <<- data.frame(
    sample_id = sample_id,
    severity = severity,
    signal_id = signal_id,
    suspected_issue = suspected_issue,
    evidence = evidence,
    recommended_action = recommended_action,
    manual_review_required = "yes",
    stringsAsFactors = FALSE
  )
}

for (i in seq_len(nrow(sample_summary))) {
  row <- sample_summary[i, , drop = FALSE]
  sample_id <- row$sample_id

  if (row$frac_below_nfeature_cutoff > 0.35 || row$frac_below_ncount_cutoff > 0.35) {
    append_triage(
      sample_id,
      "high",
      "low_complexity_burden",
      "大量细胞落在当前 QC 下限之外",
      sprintf("frac_below_nfeature=%.3f; frac_below_ncount=%.3f", row$frac_below_nfeature_cutoff, row$frac_below_ncount_cutoff),
      "优先审阅样本级原始 QC；必要时为该样本设置 sample-wise QC 阈值，而不是直接沿用全局 cutoff。"
    )
  }

  if (identical(row$mito_detection_method, "failed") || isTRUE(row$mito_feature_count == 0)) {
    append_triage(
      sample_id,
      "high",
      "mito_detection_failed",
      "未能可靠识别线粒体基因，percent.mito 当前不可信",
      sprintf(
        "mito_detection_method=%s; mito_feature_count=%s; detail=%s",
        normalize_scalar_value(row$mito_detection_method, "failed"),
        normalize_scalar_value(row$mito_feature_count, "0"),
        normalize_scalar_value(row$mito_detection_detail, "NA")
      ),
      "优先使用 reference GTF 染色体位点或 config/mito_gene_list.txt 明确线粒体基因列表，再继续解释或依赖 percent.mito。"
    )
  } else {
    if (identical(row$mito_detection_method, "prefix_fallback")) {
      append_triage(
        sample_id,
        "medium",
        "mito_prefix_fallback",
        "当前样本使用 prefix fallback 识别线粒体基因",
        sprintf("detail=%s; genes=%s", normalize_scalar_value(row$mito_detection_detail, "NA"), normalize_scalar_value(row$mito_detected_gene_names, "NA")),
        "结果可用于初步 QC，但建议尽快用 reference GTF 或用户显式列表复核。"
      )
    }

    if (grepl("mito_feature_count_low", normalize_scalar_value(row$mito_warning_codes), fixed = TRUE)) {
      append_triage(
        sample_id,
        "medium",
        "mito_feature_count_low",
        "识别到的线粒体基因数偏少",
        sprintf("mito_feature_count=%s; genes=%s", normalize_scalar_value(row$mito_feature_count, "NA"), normalize_scalar_value(row$mito_detected_gene_names, "NA")),
        "先检查 GTF 染色体命名、feature 列选择和 mito_gene_list.txt 是否需要显式补充。"
      )
    }
  }

  if (row$frac_above_mito_cutoff > 0.25) {
    append_triage(
      sample_id,
      "high",
      "high_mito_burden",
      "当前样本存在较明显的高 mt 负担",
      sprintf("frac_above_mito_cutoff=%.3f; median_percent_mito=%.3f", row$frac_above_mito_cutoff, row$median_percent_mito),
      "先判断这是全样本质量问题还是特定细胞群体的生物信号，再决定是否收紧 mt cutoff。"
    )
  }

  if (row$nfeature_density_peaks > 2 || row$mito_density_peaks > 2) {
    append_triage(
      sample_id,
      "medium",
      "possible_multimodal_qc",
      "QC 指标分布疑似多峰",
      sprintf("nfeature_peaks=%s; mito_peaks=%s", row$nfeature_density_peaks, row$mito_density_peaks),
      "先标记人工复核，避免直接用单一固定 cutoff 清洗。"
    )
  }

  if (!is.na(row$median_genes_per_cell) && !is.na(row$sequencing_saturation) && row$median_genes_per_cell < 800) {
    if (row$sequencing_saturation < 60) {
      append_triage(
        sample_id,
        "medium",
        "low_genes_low_saturation",
        "每细胞基因数偏低且测序饱和度仍偏低",
        sprintf("median_genes_per_cell=%.1f; sequencing_saturation=%.1f", row$median_genes_per_cell, row$sequencing_saturation),
        "追加测序可能仍有收益，同时继续核查上游文库复杂度。"
      )
    } else if (row$sequencing_saturation >= 80) {
      append_triage(
        sample_id,
        "medium",
        "low_genes_high_saturation",
        "每细胞基因数偏低但测序饱和度已接近平稳",
        sprintf("median_genes_per_cell=%.1f; sequencing_saturation=%.1f", row$median_genes_per_cell, row$sequencing_saturation),
        "优先回头检查样本和文库复杂度，不要默认把问题归因于测序深度不足。"
      )
    }
  }

  if (identical(row$ambient_soupx_ready, "true")) {
    append_triage(
      sample_id,
      "info",
      "ambient_soupx_ready",
      "该样本具备 SoupX 主路径条件",
      sprintf("raw_matrix_kind=%s; preferred_method=%s", normalize_scalar_value(row$raw_matrix_kind, "NA"), normalize_scalar_value(row$ambient_preferred_method, "soupx")),
      "进入 ambient branch 时应优先走 SoupX，并输出 contamination summary 与前后 marker leakage 对比。"
    )
  } else if (identical(row$ambient_decontx_ready, "true")) {
    append_triage(
      sample_id,
      "info",
      "ambient_decontx_fallback",
      "该样本缺少 SoupX 主路径条件，但具备 DecontX fallback 条件",
      sprintf("preferred_method=%s; ambient_notes=%s", normalize_scalar_value(row$ambient_preferred_method, "decontx"), normalize_scalar_value(row$ambient_notes, "NA")),
      "进入 ambient branch 时应明确标记为 DecontX fallback，而不是与 SoupX 等价。"
    )
  }

  if (identical(row$gene_id_type_resolved, "unknown")) {
    append_triage(
      sample_id,
      "high",
      "gene_id_type_unrecognized",
      "feature naming profile 未能稳定识别 gene_id_type",
      sprintf("feature_name_profile=%s; mito_feature_count=%s", row$feature_name_profile, row$mito_feature_count),
      "先确认 features.tsv 的命名列和 reference GTF 是否匹配，再继续解释 percent.mito/marker 命中。"
    )
  }

  if (identical(row$ambient_soupx_ready, "false")) {
    append_triage(
      sample_id,
      "info",
      "raw_matrix_missing_for_soupx",
      "该样本不具备 SoupX 原始 droplets 条件",
      sprintf("raw_matrix_dir=%s; preferred_method=%s", ifelse(nzchar(row$raw_matrix_dir), row$raw_matrix_dir, "NA"), normalize_scalar_value(row$ambient_preferred_method, "none")),
      "在 ambient 报告里确认是否降级到 DecontX fallback，且不要把“可做 ambient”与“已完成 correction”混为一谈。"
    )
  }

  if (length(describe_readiness_notes(row$readiness_notes)) > 0) {
    append_triage(
      sample_id,
      "info",
      "input_contract_limits",
      "输入契约限制需要在 pre-QC 审阅中显式查看",
      paste(describe_readiness_notes(row$readiness_notes), collapse = "; "),
      "结合 sample_qc_summary.tsv 一起确认这些限制是否会影响继续推进。"
    )
  }
}

project_branch_row <- if (nrow(branch_readiness_df) > 0 && "sample_id" %in% colnames(branch_readiness_df)) {
  branch_readiness_df[branch_readiness_df$sample_id == "__PROJECT__", , drop = FALSE]
} else {
  data.frame(stringsAsFactors = FALSE)
}
project_notes <- if (nrow(project_branch_row) > 0 && "notes" %in% colnames(project_branch_row)) {
  describe_readiness_notes(project_branch_row$notes[1])
} else {
  character(0)
}
if (length(project_notes) > 0) {
  append_triage(
    "__PROJECT__",
    "info",
    "project_readiness_limits",
    "项目级 intake/readiness 仍存在限制",
    paste(project_notes, collapse = "; "),
    "在进入后续整合和分支模块前，先明确这些限制是否可接受。"
  )
}

triage_df <- if (length(triage_rows) > 0) {
  bind_rows(triage_rows)
} else {
  data.frame(
    sample_id = character(0),
    severity = character(0),
    signal_id = character(0),
    suspected_issue = character(0),
    evidence = character(0),
    recommended_action = character(0),
    manual_review_required = character(0),
    stringsAsFactors = FALSE
  )
}

qc_long <- meta_df %>%
  select(sample_id, nCount_RNA, nFeature_RNA, percent.mito, percent.ribo, log10GenesPerUMI) %>%
  tidyr::pivot_longer(
    cols = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", "log10GenesPerUMI"),
    names_to = "metric",
    values_to = "value"
  )

violin_plot <- ggplot(qc_long, aes(x = sample_id, y = value, fill = sample_id)) +
  geom_violin(scale = "width", trim = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL, title = "Pre-QC metric distributions by sample")

scatter_df <- meta_df
scatter_plot <- ggplot(scatter_df, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mito)) +
  geom_point(size = 0.25, alpha = 0.6) +
  facet_wrap(~ sample_id, scales = "free") +
  scale_color_gradientn(colors = paper_feature_palette()) +
  theme_classic(base_size = 11) +
  labs(title = "nCount vs nFeature before QC", x = "nCount_RNA", y = "nFeature_RNA", color = "percent.mito")

cutoff_df <- sample_summary %>%
  select(sample_id, frac_below_nfeature_cutoff, frac_below_ncount_cutoff, frac_below_log10umi_cutoff, frac_above_mito_cutoff) %>%
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "metric",
    values_to = "fraction"
  )

cutoff_plot <- ggplot(cutoff_df, aes(x = sample_id, y = fraction, fill = metric)) +
  geom_col(position = "dodge") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Fraction of cells outside current QC cutoffs", x = NULL, y = "Fraction")

save_plot_dual(violin_plot, file.path(stage_dir, "pre_qc_violin.png"), width = 12, height = 8)
save_plot_dual(scatter_plot, file.path(stage_dir, "pre_qc_scatter.png"), width = 12, height = 8)
save_plot_dual(cutoff_plot, file.path(stage_dir, "pre_qc_cutoff_burden.png"), width = 10, height = 5)

write_tsv(sample_summary, file.path(stage_dir, "sample_qc_summary.tsv"))
write_tsv(triage_df, file.path(stage_dir, "triage.tsv"))

report_lines <- c(
  "# Pre-QC EDA Report",
  "",
  sprintf("- 样本数: `%s`", length(unique(sample_summary$sample_id))),
  sprintf("- 细胞总数: `%s`", format(sum(sample_summary$cells_raw), big.mark = ",")),
  sprintf("- triage 信号条数: `%s`", nrow(triage_df)),
  "",
  "## Key Files",
  sprintf("- `sample_qc_summary.tsv`: `%s`", file.path(stage_dir, "sample_qc_summary.tsv")),
  sprintf("- `triage.tsv`: `%s`", file.path(stage_dir, "triage.tsv")),
  sprintf("- `pre_qc_violin.png`: `%s`", file.path(stage_dir, "pre_qc_violin.png")),
  sprintf("- `pre_qc_scatter.png`: `%s`", file.path(stage_dir, "pre_qc_scatter.png")),
  sprintf("- `pre_qc_cutoff_burden.png`: `%s`", file.path(stage_dir, "pre_qc_cutoff_burden.png")),
  "",
  "## Review Focus",
  "- 先确认是否存在样本级整体质量问题，再决定是否调整 sample-wise QC 阈值。",
  "- 如果 `triage.tsv` 出现 `possible_multimodal_qc`，优先人工审阅，不要直接套全局固定 cutoff。",
  sprintf("- ambient branch 报告位置: `%s`", file.path(cfg$eda_report_dir, "ambient", "report.md")),
  "- 如果 `mito_detection_method` 不是 `gtf` 或 `user_list`，先确认 mito 基因识别是否可信，再解释 `percent.mito`。"
)

if (length(project_notes) > 0) {
  report_lines <- c(report_lines, "", "## Readiness Limits")
  for (note in project_notes) {
    report_lines <- c(report_lines, sprintf("- %s", note))
  }
}

report_lines <- c(report_lines, "", "## Sample Contracts")
for (i in seq_len(nrow(sample_summary))) {
  row <- sample_summary[i, , drop = FALSE]
  report_lines <- c(
    report_lines,
    sprintf("- `%s`: platform=`%s`; profile=`%s`; gene_id_type=`%s`; reference_version=`%s`", row$sample_id, row$platform, row$feature_name_profile, row$gene_id_type_resolved, ifelse(nzchar(row$reference_version), row$reference_version, "NA")),
    sprintf("  filtered_matrix_dir=`%s`", ifelse(nzchar(row$filtered_matrix_dir), row$filtered_matrix_dir, "NA")),
    sprintf("  raw_matrix_dir=`%s`", ifelse(nzchar(row$raw_matrix_dir), row$raw_matrix_dir, "NA")),
    sprintf("  ambient_preferred_method=`%s`; ambient_fallback_method=`%s`; soupx_ready=`%s`; decontx_ready=`%s`; cellbender_ready=`%s`", row$ambient_preferred_method, row$ambient_fallback_method, row$ambient_soupx_ready, row$ambient_decontx_ready, row$ambient_cellbender_ready),
    sprintf("  metrics_path=`%s`", ifelse(nzchar(row$metrics_summary_path), row$metrics_summary_path, "NA")),
    sprintf("  bam_path=`%s`", ifelse(nzchar(row$bam_path), row$bam_path, "NA")),
    sprintf("  mito_detection_method=`%s`; detail=`%s`; species_guess=`%s`", normalize_scalar_value(row$mito_detection_method, "failed"), ifelse(nzchar(normalize_scalar_value(row$mito_detection_detail)), row$mito_detection_detail, "NA"), normalize_scalar_value(row$species_guess, "unknown")),
    sprintf("  mito_reference_seqnames=`%s`", ifelse(nzchar(normalize_scalar_value(row$mito_reference_seqnames)), row$mito_reference_seqnames, "NA")),
    sprintf("  mito_feature_count=`%s`; ribo_feature_count=`%s`; cell_cycle_s_features=`%s`; cell_cycle_g2m_features=`%s`", row$mito_feature_count, row$ribo_feature_count, row$cell_cycle_s_feature_count, row$cell_cycle_g2m_feature_count),
    sprintf("  mito_gene_list=`%s`", ifelse(nzchar(normalize_scalar_value(row$mito_detected_gene_names)), row$mito_detected_gene_names, "NA")),
    sprintf("  mito_warnings=`%s`", ifelse(nzchar(normalize_scalar_value(row$mito_warning_messages)), row$mito_warning_messages, "NA")),
    sprintf("  ambient_notes=`%s`", ifelse(nzchar(normalize_scalar_value(row$ambient_notes)), row$ambient_notes, "NA"))
  )
}

if (nrow(triage_df) > 0) {
  report_lines <- c(report_lines, "", "## Triage Summary")
  for (i in seq_len(nrow(triage_df))) {
    row <- triage_df[i, , drop = FALSE]
    report_lines <- c(
      report_lines,
      sprintf("- `%s` [%s] `%s`: %s", row$sample_id, row$severity, row$signal_id, row$recommended_action)
    )
  }
}

write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("pre-QC EDA 已输出到: ", stage_dir)
