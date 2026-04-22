source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "annotation_stats_helpers.R"))

suppressPackageStartupMessages({
  library(dplyr)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

manifest_path <- file.path(cfg$table_dir, "annotation", "layer_annotation_manifest.tsv")
if (!file.exists(manifest_path)) {
  stop(sprintf("缺少 annotation manifest: %s", manifest_path), call. = FALSE)
}

manifest_df <- read_tsv_optional(manifest_path)
if (nrow(manifest_df) == 0) {
  stop("annotation manifest 为空。", call. = FALSE)
}

stage_dir <- ensure_eda_stage_dir(cfg, "annotation")
summary_rows <- list()

for (i in seq_len(nrow(manifest_df))) {
  row <- manifest_df[i, , drop = FALSE]
  annotation_path <- row$annotation_table_tsv[[1]]
  evidence_path <- file.path(dirname(annotation_path), "annotation_evidence.tsv")
  ann_df <- read_tsv_optional(annotation_path)
  evidence_df <- read_tsv_optional(evidence_path)

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    layer_id = row$layer_id[[1]],
    layer_role = row$layer_role[[1]],
    cluster_count = nrow(ann_df),
    determined_n = sum(ann_df$confidence == "确定"),
    tentative_n = sum(ann_df$confidence == "暂定"),
    undetermined_n = sum(ann_df$confidence == "未定"),
    panel_validated_clusters = sum(ann_df$evidence_relation %in% c("一致", "部分一致", "冲突")),
    data_only_clusters = sum(ann_df$evidence_relation %in% c("仅数据驱动", "无文献验证")),
    marker_panel_available = if (nrow(ann_df) > 0) unique(ann_df$marker_panel_available)[1] else "no",
    conflict_n = sum(ann_df$evidence_relation == "冲突"),
    annotation_table_tsv = annotation_path,
    evidence_table_tsv = evidence_path,
    report_md = row$report_md[[1]],
    stringsAsFactors = FALSE
  )
}

summary_df <- bind_rows(summary_rows)
summary_path <- file.path(cfg$table_dir, "annotation", "layer_annotation_summary.tsv")
write_tsv(summary_df, summary_path)

report_lines <- c(
  "# Annotation EDA Summary",
  "",
  sprintf("- annotation manifest: `%s`", manifest_path),
  sprintf("- annotation summary: `%s`", summary_path),
  "",
  "## Layer Overview"
)

for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, , drop = FALSE]
  report_lines <- c(
    report_lines,
    sprintf(
      "- `%s` (%s): clusters=%s, 确定=%s, 暂定=%s, 未定=%s, conflicts=%s, marker_panel_available=%s",
      row$layer_id,
      row$layer_role,
      row$cluster_count,
      row$determined_n,
      row$tentative_n,
      row$undetermined_n,
      row$conflict_n,
      row$marker_panel_available
    ),
    sprintf("  annotation_table: %s", row$annotation_table_tsv),
    sprintf("  evidence_table: %s", row$evidence_table_tsv),
    sprintf("  layer_report: %s", row$report_md)
  )
}

data_only_layers <- summary_df$layer_id[summary_df$data_only_clusters > 0]
if (length(data_only_layers) > 0) {
  report_lines <- c(
    report_lines,
    "",
    "## Data-Driven Only Layers",
    sprintf("- `%s`", paste(data_only_layers, collapse = ", "))
  )
}

write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("annotation EDA summary 已输出到: ", stage_dir)
