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
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03b_integration_eda"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03a2 <- read_manifest_local(cfg$module_03a2_manifest_path)
candidate_index_tsv <- resolve_output_local(manifest_03a2, "candidate_index_tsv")
candidate_index <- read_tsv_optional(candidate_index_tsv)
if (nrow(candidate_index) == 0) {
  stop(sprintf("候选索引为空: %s", candidate_index_tsv), call. = FALSE)
}

r2_by_factor <- function(values, group) {
  group <- as.factor(group)
  if (length(unique(stats::na.omit(group))) < 2) {
    return(NA_real_)
  }
  fit <- stats::lm(values ~ group)
  summary(fit)$r.squared
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  max(x)
}

sample_rows <- function(n, max_n = 1500L) {
  if (n <= max_n) {
    return(seq_len(n))
  }
  sort(sample(seq_len(n), max_n))
}

mean_silhouette_by_group <- function(emb, group) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    return(NA_real_)
  }
  group <- as.character(group)
  keep <- !is.na(group) & nzchar(group)
  if (sum(keep) < 3 || length(unique(group[keep])) < 2) {
    return(NA_real_)
  }
  emb <- emb[keep, , drop = FALSE]
  group <- group[keep]
  idx <- sample_rows(nrow(emb))
  emb <- emb[idx, , drop = FALSE]
  group <- group[idx]
  if (length(unique(group)) < 2) {
    return(NA_real_)
  }
  sil <- cluster::silhouette(as.integer(factor(group)), stats::dist(emb))
  mean(sil[, "sil_width"], na.rm = TRUE)
}

same_sample_knn_fraction <- function(emb, group, k = 5L) {
  group <- as.character(group)
  keep <- !is.na(group) & nzchar(group)
  if (sum(keep) <= k || length(unique(group[keep])) < 2) {
    return(NA_real_)
  }
  emb <- emb[keep, , drop = FALSE]
  group <- group[keep]
  idx <- sample_rows(nrow(emb))
  emb <- emb[idx, , drop = FALSE]
  group <- group[idx]
  d <- as.matrix(stats::dist(emb))
  diag(d) <- Inf
  k <- min(k, nrow(d) - 1L)
  nn <- t(apply(d, 1, function(x) order(x)[seq_len(k)]))
  mean(vapply(seq_len(nrow(nn)), function(i) mean(group[nn[i, ]] == group[i]), numeric(1)), na.rm = TRUE)
}

candidate_metric_rows <- list()
umap_rows <- list()

for (idx in seq_len(nrow(candidate_index))) {
  row <- candidate_index[idx, , drop = FALSE]
  seu <- readRDS(row$out_rds[[1]])
  seu <- maybe_join_layers(seu)
  candidate_id <- sprintf("%s__%s", row$normalization[[1]], row$integration[[1]])
  reduction_name <- row$reduction_name[[1]]
  umap_name <- row$umap_name[[1]]
  emb <- Embeddings(seu, reduction = reduction_name)
  pc_ids <- seq_len(min(20L, ncol(emb)))
  sample_r2 <- vapply(pc_ids, function(i) r2_by_factor(emb[, i], seu$orig.ident), numeric(1))
  group_var <- if ("analysis_group" %in% colnames(seu@meta.data)) "analysis_group" else if ("condition" %in% colnames(seu@meta.data)) "condition" else ""
  group_r2 <- if (nzchar(group_var)) vapply(pc_ids, function(i) r2_by_factor(emb[, i], seu@meta.data[[group_var]]), numeric(1)) else rep(NA_real_, length(pc_ids))

  metric_emb <- if (umap_name %in% Reductions(seu)) Embeddings(seu, reduction = umap_name) else emb[, seq_len(min(2L, ncol(emb))), drop = FALSE]
  candidate_metric_rows[[length(candidate_metric_rows) + 1]] <- data.frame(
    candidate_id = candidate_id,
    layer_id = row$layer_id[[1]],
    normalization = row$normalization[[1]],
    integration = row$integration[[1]],
    reduction_name = reduction_name,
    umap_name = umap_name,
    max_sample_r2 = safe_max(sample_r2),
    max_group_r2 = safe_max(group_r2),
    silhouette_by_orig_ident = mean_silhouette_by_group(metric_emb, seu$orig.ident),
    same_sample_knn_fraction = same_sample_knn_fraction(metric_emb, seu$orig.ident),
    runtime_sec = suppressWarnings(as.numeric(row$runtime_sec[[1]])),
    downgrade_reason = normalize_scalar_value(row$downgrade_reason[[1]]),
    stringsAsFactors = FALSE
  )

  if (ncol(metric_emb) >= 2) {
    plot_df <- as.data.frame(metric_emb[, 1:2, drop = FALSE])
    colnames(plot_df)[1:2] <- c("UMAP_1", "UMAP_2")
    plot_df$candidate_id <- candidate_id
    plot_df$orig.ident <- as.character(seu$orig.ident)
    plot_df$group <- if (nzchar(group_var)) as.character(seu@meta.data[[group_var]]) else ""
    umap_rows[[length(umap_rows) + 1]] <- plot_df
  }
}

compare_df <- dplyr::bind_rows(candidate_metric_rows)
rank_metric <- function(x, decreasing = FALSE) {
  if (all(is.na(x))) {
    return(rep(1, length(x)))
  }
  value <- x
  fill <- if (decreasing) min(value, na.rm = TRUE) - 1 else max(value, na.rm = TRUE) + 1
  value[is.na(value)] <- fill
  rank(if (decreasing) -value else value, ties.method = "min")
}
compare_df$rank_score <- rank_metric(compare_df$max_sample_r2) +
  rank_metric(compare_df$silhouette_by_orig_ident) +
  rank_metric(compare_df$same_sample_knn_fraction)
compare_df <- dplyr::arrange(compare_df, rank_score, max_sample_r2, silhouette_by_orig_ident)
recommended <- compare_df[1, , drop = FALSE]
recommended_value <- recommended$candidate_id[[1]]

read_selected_value <- function(path) {
  if (!file.exists(path)) {
    return("")
  }
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) == 0) "" else lines[[1]]
}

selected_value <- read_selected_value(cfg$selected_integration_file)
if (nzchar(selected_value)) {
  if (!selected_value %in% compare_df$candidate_id) {
    stop(sprintf("selected_integration.txt 中的选择不合法: %s。合法值: %s", selected_value, paste(compare_df$candidate_id, collapse = ", ")), call. = FALSE)
  }
} else {
  ensure_dir(dirname(cfg$selected_integration_file))
  writeLines(
    c(
      "# automatic fallback recommendation from 03b_integration_eda",
      "# review the 03b report and approve the integration gate before 03c",
      recommended_value
    ),
    cfg$selected_integration_file,
    useBytes = TRUE
  )
  selected_value <- recommended_value
}

triage_rows <- list()
append_triage <- function(severity, signal_id, evidence) {
  triage_rows[[length(triage_rows) + 1]] <<- make_triage_row(
    sample_id = panorama_spec$layer_id,
    severity = severity,
    signal_id = signal_id,
    evidence = evidence
  )
}
for (i in seq_len(nrow(compare_df))) {
  if (nzchar(compare_df$downgrade_reason[[i]])) {
    append_triage("medium", "integration_method_unavailable", sprintf("%s: %s", compare_df$candidate_id[[i]], compare_df$downgrade_reason[[i]]))
  }
}
if (is.finite(recommended$max_sample_r2[[1]]) && recommended$max_sample_r2[[1]] >= 0.30) {
  append_triage("medium", "integration_undercorrection_suspected", sprintf("recommended=%s; max_sample_r2=%.3f", recommended_value, recommended$max_sample_r2[[1]]))
}
none_rows <- compare_df[compare_df$integration == "none", , drop = FALSE]
if (nrow(none_rows) > 0 && is.finite(none_rows$max_group_r2[[1]]) && is.finite(recommended$max_group_r2[[1]]) &&
    none_rows$max_group_r2[[1]] >= 0.15 && recommended$max_group_r2[[1]] < none_rows$max_group_r2[[1]] * 0.5) {
  append_triage("medium", "integration_overcorrection_suspected", sprintf("none max_group_r2=%.3f; recommended max_group_r2=%.3f", none_rows$max_group_r2[[1]], recommended$max_group_r2[[1]]))
}
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

compare_tsv <- file.path(cfg$integration_report_dir_layer, "integration_compare.tsv")
triage_tsv <- file.path(cfg$integration_report_dir_layer, "integration_triage.tsv")
umap_grid_png <- file.path(cfg$integration_report_dir_layer, "candidate_umap_grid.png")
r2_compare_png <- file.path(cfg$integration_report_dir_layer, "candidate_r2_compare.png")
report_path <- file.path(cfg$integration_report_dir_layer, "report.md")

write_tsv_local(compare_df, compare_tsv)
write_tsv_local(triage_df, triage_tsv)

umap_df <- dplyr::bind_rows(umap_rows)
if (nrow(umap_df) > 0) {
  umap_plot <- ggplot2::ggplot(umap_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = orig.ident)) +
    ggplot2::geom_point(size = 0.2, alpha = 0.65) +
    ggplot2::facet_wrap(~ candidate_id, scales = "free") +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::labs(title = "Integration candidates by sample", color = NULL)
} else {
  umap_plot <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No candidate UMAP available")
}
save_plot_local(umap_plot, umap_grid_png, width = 10, height = 7)

r2_long <- rbind(
  data.frame(candidate_id = compare_df$candidate_id, metric = "max_sample_r2", value = compare_df$max_sample_r2),
  data.frame(candidate_id = compare_df$candidate_id, metric = "max_group_r2", value = compare_df$max_group_r2)
)
r2_plot <- ggplot2::ggplot(r2_long, ggplot2::aes(x = candidate_id, y = value, fill = metric)) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::coord_flip() +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::labs(title = "Candidate PCA association", x = NULL, y = "R2", fill = NULL)
save_plot_local(r2_plot, r2_compare_png, width = 8, height = 5)

report_lines <- c(
  "# 03b Integration EDA Report",
  "",
  sprintf("- layer_id: `%s`", panorama_spec$layer_id),
  sprintf("- candidate_index: `%s`", candidate_index_tsv),
  sprintf("- selected_integration_file: `%s`", cfg$selected_integration_file),
  sprintf("- current selected value: `%s`", selected_value),
  sprintf("- automatic fallback recommendation: `%s`", recommended_value),
  "",
  "## Candidate Ranking",
  render_markdown_table_local(compare_df),
  "",
  "## Triage",
  render_triage_markdown(triage_df),
  "",
  "## Review Action",
  sprintf("Edit `%s` so the first non-comment line is the desired `<normalization>__<integration>` value, then approve the `integration` gate.", cfg$selected_integration_file)
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_03b_manifest_path,
  new_outputs = list(
    report = build_output_entry(report_path, "md", module_name, "human-readable integration EDA report", base_dir = cfg$project_root),
    integration_compare_tsv = build_output_entry(compare_tsv, "tsv", module_name, "one row per integration candidate diagnostic summary", base_dir = cfg$project_root, schema = infer_schema_from_df(compare_df)),
    integration_triage_tsv = build_output_entry(triage_tsv, "tsv", module_name, "one row per integration triage signal", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
    umap_grid_png = build_output_entry(umap_grid_png, "png", module_name, "faceted candidate UMAP by sample", base_dir = cfg$project_root),
    r2_compare_png = build_output_entry(r2_compare_png, "png", module_name, "candidate sample/group R2 comparison", base_dir = cfg$project_root),
    selected_integration_txt = build_output_entry(cfg$selected_integration_file, "txt", module_name, "first non-comment line selects normalization__integration for 03c", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(module_03a2_manifest = cfg$module_03a2_manifest_path, candidate_index_tsv = candidate_index_tsv),
  version = cfg$module_version,
  depends_on = list(module_03a2 = cfg$module_03a2_manifest_path)
)

message("03b 完成。推荐: ", recommended_value, "；报告: ", report_path)
