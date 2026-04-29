#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source_utf8 <- function(path) source(path, encoding = "UTF-8")

source_utf8(file.path(.script_dir, "helpers", "metadata_dsl_parser.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_comparison.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_communication.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_regulation.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_trajectory.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_enrichment.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_fanout_st.R"))

env_or <- function(name, default = "") {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) default else value
}

project_root <- env_or("PROJECT_ROOT", getwd())
metadata_dir <- env_or("METADATA_DIR", file.path(project_root, "metadata"))
results_dir <- env_or("RESULTS_DIR", file.path(project_root, "results"))
question_path <- env_or("ANALYSIS_QUESTIONS_FILE", file.path(metadata_dir, "analysis_questions.tsv"))
summary_path <- env_or("METADATA_GENERATOR_SUMMARY", file.path(results_dir, "00_validation", "metadata_generator_summary.tsv"))

comparison_axes <- c("cluster_marker", "directional_DEG", "pairwise", "contrast_only", "composition")
communication_axes <- c("bidirectional", "directional", "sequential", "symmetric", "pairwise_comm")
regulation_axes <- c("regulation_per", "regulation_pair", "regulation_stage")
trajectory_axes <- c("lineage", "velocity")
st_axes <- c("spatial_clustering", "spatial_neighbor", "spatial_overlay", "deconv_pair", "SVG", "SVG_diff", "deconv_validation")

append_df <- function(items, df) {
  if (nrow(df) > 0) {
    items[[length(items) + 1L]] <- df
  }
  items
}

bind_or_empty <- function(items, cols) {
  if (length(items) == 0) {
    return(m3_empty_df(cols))
  }
  df <- do.call(rbind, items)
  rownames(df) <- NULL
  for (col in cols) {
    if (!col %in% colnames(df)) df[[col]] <- ""
  }
  df[, cols, drop = FALSE]
}

questions <- m3_read_tsv(question_path)
if (nrow(questions) == 0) {
  stop(sprintf("analysis questions table is missing or empty: %s", question_path), call. = FALSE)
}
required <- c("question_id", "scope", "sender_groups", "receiver_groups", "condition_split", "contrast_axis", "tools_to_run", "status")
missing <- setdiff(required, colnames(questions))
if (length(missing) > 0) {
  stop(sprintf("analysis questions table is missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

active <- questions[questions$status == "active", , drop = FALSE]

comparison_items <- list()
communication_items <- list()
scenic_items <- list()
trajectory_items <- list()

for (idx in seq_len(nrow(active))) {
  q <- active[idx, , drop = FALSE]
  axis <- q$contrast_axis[[1]]
  if (axis %in% comparison_axes) {
    comparison_items <- append_df(comparison_items, m3_fanout_comparison(q))
  } else if (axis %in% communication_axes) {
    communication_items <- append_df(communication_items, m3_fanout_communication(q))
  } else if (axis %in% regulation_axes) {
    scenic_items <- append_df(scenic_items, m3_fanout_regulation(q))
  } else if (axis %in% trajectory_axes) {
    trajectory_items <- append_df(trajectory_items, m3_fanout_trajectory(q))
  } else if (axis %in% st_axes) {
    m3_fanout_st(q)
  } else if (identical(axis, "enrichment_target")) {
    next
  } else {
    stop(sprintf("No M3 fan-out handler for contrast_axis=%s question_id=%s", axis, q$question_id[[1]]), call. = FALSE)
  }
}

comparisons <- bind_or_empty(comparison_items, m3_comparison_cols)
communication_pairs <- bind_or_empty(communication_items, m3_communication_cols)
trajectory_pairs <- bind_or_empty(trajectory_items, m3_trajectory_cols)
scenic_targets <- bind_or_empty(scenic_items, m3_scenic_cols)
enrichment_targets <- m3_fanout_enrichment(questions, comparisons)
deconv_pairs <- m3_empty_df(m3_deconv_cols)
spatial_pairs <- m3_empty_df(m3_spatial_cols)

m3_assert_unique(comparisons, "comparison_id", "comparisons.tsv")
m3_assert_unique(communication_pairs, "pair_id", "communication_pairs.tsv")
m3_assert_unique(trajectory_pairs, "trajectory_id", "trajectory_pairs.tsv")
m3_assert_unique(scenic_targets, "target_id", "scenic_targets.tsv")
m3_assert_unique(enrichment_targets, "target_id", "enrichment_targets.tsv")

m3_write_generated_tsv(comparisons, file.path(metadata_dir, "comparisons.tsv"), m3_comparison_cols)
m3_write_generated_tsv(communication_pairs, file.path(metadata_dir, "communication_pairs.tsv"), m3_communication_cols)
m3_write_generated_tsv(trajectory_pairs, file.path(metadata_dir, "trajectory_pairs.tsv"), m3_trajectory_cols)
m3_write_generated_tsv(scenic_targets, file.path(metadata_dir, "scenic_targets.tsv"), m3_scenic_cols)
m3_write_generated_tsv(enrichment_targets, file.path(metadata_dir, "enrichment_targets.tsv"), m3_enrichment_cols)
m3_write_generated_tsv(deconv_pairs, file.path(metadata_dir, "deconv_pairs.tsv"), m3_deconv_cols)
m3_write_generated_tsv(spatial_pairs, file.path(metadata_dir, "spatial_pairs.tsv"), m3_spatial_cols)

summary <- data.frame(
  table = c("comparisons.tsv", "communication_pairs.tsv", "trajectory_pairs.tsv", "scenic_targets.tsv", "enrichment_targets.tsv", "deconv_pairs.tsv", "spatial_pairs.tsv"),
  rows = c(nrow(comparisons), nrow(communication_pairs), nrow(trajectory_pairs), nrow(scenic_targets), nrow(enrichment_targets), nrow(deconv_pairs), nrow(spatial_pairs)),
  stringsAsFactors = FALSE
)
dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(summary, summary_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat("metadata generator summary\n")
for (idx in seq_len(nrow(summary))) {
  cat(sprintf("%s\t%s\n", summary$table[[idx]], summary$rows[[idx]]))
}
cat(sprintf("SUMMARY=%s\n", summary_path))
