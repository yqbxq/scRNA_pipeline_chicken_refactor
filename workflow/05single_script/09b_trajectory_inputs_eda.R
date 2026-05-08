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

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_04.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_05.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_09.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "jsonlite"))

cfg <- get_single_script_config_09()
module_name <- "09b_trajectory_inputs_eda"
prepare_dirs_09(cfg)

empty_df_09b <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols), stringsAsFactors = FALSE)
}

resolve_09a_output <- function(manifest, key) {
  path <- resolve_output_local(manifest, key)
  if (!file.exists(path)) {
    stop(sprintf("09a manifest output missing on disk: %s -> %s", key, path), call. = FALSE)
  }
  path
}

input_status_cols_09b <- c(
  "pair_id", "source_question_id", "layer_scope", "split_var", "split_value",
  "cell_n_before", "outlier_n", "cell_n_after", "retention_fraction",
  "root_group_present", "terminal_group_present", "cc_status", "pca_npcs",
  "input_rds", "status", "reason"
)

pre_post_cols_09b <- c(
  "pair_id", "split_value", "metric", "stage", "n_cells",
  "median", "p05", "p95"
)

split_balance_cols_09b <- c(
  "pair_id", "split_var", "split_value", "label_var", "label_value",
  "n_cells", "fraction", "total_cells", "imbalance_fraction",
  "warn_imbalanced_split", "warn_low_label_count"
)

umap_index_cols_09b <- c(
  "pair_id", "split_value", "color_var", "figure_path", "status", "reason"
)

manifest_09a <- read_manifest_local(cfg$module_09a_manifest_path)
input_index_path <- resolve_09a_output(manifest_09a, "trajectory_input_index")
outlier_cells_path <- resolve_09a_output(manifest_09a, "outlier_qc_cells")
split_index_path <- resolve_09a_output(manifest_09a, "split_index")
cc_summary_path <- resolve_09a_output(manifest_09a, "cc_score_summary")
pca_summary_path <- resolve_09a_output(manifest_09a, "pca_summary")

input_index <- read_tsv_optional(input_index_path)
outlier_cells <- read_tsv_optional(outlier_cells_path)
split_index <- read_tsv_optional(split_index_path)
cc_summary <- read_tsv_optional(cc_summary_path)
pca_summary <- read_tsv_optional(pca_summary_path)

if (nrow(input_index) == 0) {
  input_status <- empty_df_09b(input_status_cols_09b)
} else {
  for (col in input_status_cols_09b) {
    if (!col %in% colnames(input_index)) {
      input_index[[col]] <- ""
    }
  }
  input_status <- input_index[, input_status_cols_09b, drop = FALSE]
}
write_tsv_local(input_status, cfg$trajectory_input_status_tsv)

summarize_metric_stage_09b <- function(df, metric, stage) {
  values <- suppressWarnings(as.numeric(df[[metric]]))
  data.frame(
    metric = metric,
    stage = stage,
    n_cells = sum(is.finite(values)),
    median = safe_median(values),
    p05 = safe_quantile(values, 0.05),
    p95 = safe_quantile(values, 0.95),
    stringsAsFactors = FALSE
  )
}

build_outlier_pre_post_09b <- function(outlier_cells) {
  if (nrow(outlier_cells) == 0) {
    return(empty_df_09b(pre_post_cols_09b))
  }
  metrics <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "neighbor_distance_z"), colnames(outlier_cells))
  if (length(metrics) == 0) {
    return(empty_df_09b(pre_post_cols_09b))
  }
  groups <- split(outlier_cells, paste(outlier_cells$pair_id, outlier_cells$split_value, sep = "\r"))
  rows <- list()
  for (group in groups) {
    pair_id <- group$pair_id[[1]]
    split_value <- group$split_value[[1]]
    retained <- tolower(as.character(group$retained)) %in% c("true", "yes", "1")
    for (metric in metrics) {
      pre <- summarize_metric_stage_09b(group, metric, "pre")
      post <- summarize_metric_stage_09b(group[retained, , drop = FALSE], metric, "post")
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(pair_id = pair_id, split_value = split_value, stringsAsFactors = FALSE),
        pre
      )
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(pair_id = pair_id, split_value = split_value, stringsAsFactors = FALSE),
        post
      )
    }
  }
  dplyr::bind_rows(rows)
}

outlier_pre_post <- build_outlier_pre_post_09b(outlier_cells)
write_tsv_local(outlier_pre_post, cfg$trajectory_outlier_pre_post_summary_tsv)

plot_outlier_pre_post_09b <- function(outlier_cells, path) {
  ensure_dir(dirname(path))
  metrics <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "neighbor_distance_z"), colnames(outlier_cells))
  if (nrow(outlier_cells) == 0 || length(metrics) == 0) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No outlier QC data") +
      ggplot2::theme_void()
    save_plot_local(plot_obj, path, width = 6, height = 3)
    return(invisible(path))
  }
  retained <- tolower(as.character(outlier_cells$retained)) %in% c("true", "yes", "1")
  pre <- outlier_cells
  pre$stage <- "pre"
  post <- outlier_cells[retained, , drop = FALSE]
  post$stage <- "post"
  combined <- rbind(pre, post)
  rows <- list()
  for (metric in metrics) {
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = combined$pair_id,
      split_value = combined$split_value,
      stage = combined$stage,
      metric = metric,
      value = suppressWarnings(as.numeric(combined[[metric]])),
      stringsAsFactors = FALSE
    )
  }
  long <- dplyr::bind_rows(rows)
  long <- long[is.finite(long$value), , drop = FALSE]
  if (nrow(long) == 0) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No finite QC metrics") +
      ggplot2::theme_void()
  } else {
    long$input_id <- paste(long$pair_id, long$split_value, sep = " / ")
    plot_obj <- ggplot2::ggplot(long, ggplot2::aes(x = stage, y = value, fill = stage)) +
      ggplot2::geom_violin(scale = "width", trim = TRUE, na.rm = TRUE) +
      ggplot2::geom_boxplot(width = 0.14, outlier.size = 0.2, na.rm = TRUE) +
      ggplot2::facet_grid(metric ~ input_id, scales = "free_y") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
      ggplot2::labs(x = NULL, y = NULL)
  }
  save_plot_local(plot_obj, path, width = 10, height = 7)
  invisible(path)
}

plot_outlier_pre_post_09b(outlier_cells, cfg$trajectory_inputs_qc_violin_png)

triage_rows <- list()
if (nrow(input_status) > 0) {
  failed <- input_status[input_status$status != "ok", , drop = FALSE]
  if (nrow(failed) > 0) {
    for (i in seq_len(nrow(failed))) {
      triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
        sample_id = failed$pair_id[[i]],
        severity = if (identical(failed$status[[i]], "failed")) "error" else "warning",
        signal_id = "trajectory_input_not_ok",
        evidence = sprintf(
          "split=%s status=%s reason=%s",
          failed$split_value[[i]],
          failed$status[[i]],
          failed$reason[[i]]
        ),
        recommended_action = "检查对应 layer_scope、split 变量、label 列和 Seurat 输入对象后重跑 09a。"
      )
    }
  }
}

build_split_balance_09b <- function(input_index, split_index) {
  if (nrow(input_index) == 0 || nrow(split_index) == 0) {
    return(empty_df_09b(split_balance_cols_09b))
  }
  ok_inputs <- input_index[input_index$status == "ok", , drop = FALSE]
  if (nrow(ok_inputs) == 0) {
    return(empty_df_09b(split_balance_cols_09b))
  }
  rows <- list()
  for (pair_id in unique(ok_inputs$pair_id)) {
    pair_inputs <- ok_inputs[ok_inputs$pair_id == pair_id, , drop = FALSE]
    pair_splits <- split_index[split_index$pair_id == pair_id, , drop = FALSE]
    if (nrow(pair_splits) == 0) {
      next
    }
    coarse_var <- pair_inputs$coarse_label_var[[1]]
    pair_splits <- pair_splits[pair_splits$label_var == coarse_var, , drop = FALSE]
    if (nrow(pair_splits) == 0) {
      next
    }
    totals <- unique(pair_inputs[, c("split_value", "cell_n_after"), drop = FALSE])
    totals$cell_n_after <- suppressWarnings(as.numeric(totals$cell_n_after))
    imbalance <- NA_real_
    warn_imbalanced <- FALSE
    split_var <- normalize_scalar_value(pair_inputs$split_var[[1]])
    if (nzchar(split_var) && nrow(totals) > 1 && max(totals$cell_n_after, na.rm = TRUE) > 0) {
      imbalance <- (max(totals$cell_n_after, na.rm = TRUE) - min(totals$cell_n_after, na.rm = TRUE)) /
        max(totals$cell_n_after, na.rm = TRUE)
      warn_imbalanced <- is.finite(imbalance) && imbalance >= cfg$trajectory_balance_warn_fraction
      if (warn_imbalanced) {
        triage_rows[[length(triage_rows) + 1L]] <<- make_triage_row(
          sample_id = pair_id,
          severity = "warning",
          signal_id = "trajectory_split_imbalanced",
          evidence = sprintf("split cell-count imbalance %.3f >= %.3f", imbalance, cfg$trajectory_balance_warn_fraction),
          recommended_action = "审阅 split_balance 表；必要时改用 pooled trajectory 或下调 split 分析优先级。"
        )
      }
    }
    total_lookup <- setNames(totals$cell_n_after, totals$split_value)
    pair_splits$n_cells <- suppressWarnings(as.numeric(pair_splits$n_cells))
    for (i in seq_len(nrow(pair_splits))) {
      total_cells <- unname(total_lookup[pair_splits$split_value[[i]]])
      low_label <- is.finite(pair_splits$n_cells[[i]]) && pair_splits$n_cells[[i]] < cfg$trajectory_split_min_cells
      if (low_label) {
        triage_rows[[length(triage_rows) + 1L]] <<- make_triage_row(
          sample_id = pair_id,
          severity = "warning",
          signal_id = "trajectory_split_label_low_cells",
          evidence = sprintf(
            "%s=%s has %s cells in split %s; threshold=%s",
            pair_splits$label_var[[i]],
            pair_splits$label_value[[i]],
            pair_splits$n_cells[[i]],
            pair_splits$split_value[[i]],
            cfg$trajectory_split_min_cells
          ),
          recommended_action = "低细胞数 label 的拟时序结果只作探索；必要时合并相邻 label 或改用 pooled 输入。"
        )
      }
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = pair_id,
        split_var = split_var,
        split_value = pair_splits$split_value[[i]],
        label_var = pair_splits$label_var[[i]],
        label_value = pair_splits$label_value[[i]],
        n_cells = pair_splits$n_cells[[i]],
        fraction = suppressWarnings(as.numeric(pair_splits$fraction[[i]])),
        total_cells = total_cells,
        imbalance_fraction = imbalance,
        warn_imbalanced_split = warn_imbalanced,
        warn_low_label_count = low_label,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(empty_df_09b(split_balance_cols_09b))
  }
  dplyr::bind_rows(rows)
}

split_balance <- build_split_balance_09b(input_index, split_index)
write_tsv_local(split_balance, cfg$trajectory_split_balance_summary_tsv)
if (nrow(split_balance) > 0) {
  for (pair_id in unique(split_balance$pair_id)) {
    write_tsv_local(
      split_balance[split_balance$pair_id == pair_id, , drop = FALSE],
      trajectory_split_balance_path_09(cfg, pair_id)
    )
  }
}

plot_umap_for_input_09b <- function(input_row) {
  pair_id <- input_row$pair_id[[1]]
  split_value <- input_row$split_value[[1]]
  input_rds <- normalize_scalar_value(input_row$input_rds[[1]])
  if (!nzchar(input_rds) || !file.exists(input_rds)) {
    return(data.frame(
      pair_id = pair_id,
      split_value = split_value,
      color_var = "",
      figure_path = "",
      status = "missing_rds",
      reason = input_rds,
      stringsAsFactors = FALSE
    ))
  }

  seu <- readRDS(input_rds)
  emb <- tryCatch(Seurat::Embeddings(seu, reduction = "umap"), error = function(e) NULL)
  if (is.null(emb) || ncol(emb) < 2) {
    return(data.frame(
      pair_id = pair_id,
      split_value = split_value,
      color_var = "",
      figure_path = "",
      status = "missing_umap",
      reason = "input object has no umap reduction",
      stringsAsFactors = FALSE
    ))
  }
  color_vars <- unique(c(
    normalize_scalar_value(input_row$coarse_label_var[[1]]),
    "group_id",
    normalize_scalar_value(input_row$split_var[[1]]),
    "Phase"
  ))
  color_vars <- color_vars[nzchar(color_vars) & color_vars %in% colnames(seu@meta.data)]
  if (length(color_vars) == 0) {
    color_vars <- "trajectory_pair_id"
  }
  fig_dir <- trajectory_pair_figure_dir_09(cfg, pair_id, if (identical(split_value, "pooled")) "" else split_value)
  ensure_dir(fig_dir)
  rows <- list()
  for (color_var in color_vars) {
    plot_df <- data.frame(
      UMAP_1 = emb[, 1],
      UMAP_2 = emb[, 2],
      color = as.character(seu@meta.data[rownames(emb), color_var]),
      stringsAsFactors = FALSE
    )
    out_png <- file.path(fig_dir, sprintf("umap__%s.png", safe_id_09(color_var)))
    plot_obj <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = color)) +
      ggplot2::geom_point(size = 0.25, alpha = 0.75) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(legend.position = "right") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
      ggplot2::labs(color = color_var, title = paste(pair_id, split_value), x = "UMAP 1", y = "UMAP 2")
    save_plot_local(plot_obj, out_png, width = 6, height = 4.8)
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = pair_id,
      split_value = split_value,
      color_var = color_var,
      figure_path = out_png,
      status = "ok",
      reason = "",
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

umap_rows <- list()
ok_inputs <- input_index[input_index$status == "ok", , drop = FALSE]
if (nrow(ok_inputs) > 0) {
  for (i in seq_len(nrow(ok_inputs))) {
    umap_rows[[length(umap_rows) + 1L]] <- tryCatch(
      plot_umap_for_input_09b(ok_inputs[i, , drop = FALSE]),
      error = function(e) data.frame(
        pair_id = ok_inputs$pair_id[[i]],
        split_value = ok_inputs$split_value[[i]],
        color_var = "",
        figure_path = "",
        status = "failed",
        reason = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    )
  }
}
umap_index <- if (length(umap_rows) > 0) dplyr::bind_rows(umap_rows) else empty_df_09b(umap_index_cols_09b)
write_tsv_local(umap_index, cfg$trajectory_umap_figure_index_tsv)

if (nrow(umap_index) > 0) {
  not_ok <- umap_index[umap_index$status != "ok", , drop = FALSE]
  for (i in seq_len(nrow(not_ok))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = not_ok$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_umap_plot_failed",
      evidence = sprintf("split=%s status=%s reason=%s", not_ok$split_value[[i]], not_ok$status[[i]], not_ok$reason[[i]]),
      recommended_action = "确认 09a 输出对象包含 umap reduction；必要时检查 PCA dims 和 cell 数。"
    )
  }
}

triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
write_tsv_local(triage_df, cfg$trajectory_inputs_triage_tsv)

status_counts <- if (nrow(input_status) > 0) {
  as.data.frame(table(input_status$status), stringsAsFactors = FALSE)
} else {
  data.frame(Var1 = character(0), Freq = integer(0), stringsAsFactors = FALSE)
}
colnames(status_counts) <- c("status", "n")
ok_n <- sum(input_status$status == "ok")
failed_n <- sum(input_status$status != "ok")
removed_n <- sum(suppressWarnings(as.numeric(input_status$outlier_n)), na.rm = TRUE)

warning_balance <- split_balance[split_balance$warn_imbalanced_split == TRUE | split_balance$warn_low_label_count == TRUE, , drop = FALSE]
cc_display <- cc_summary
if (nrow(cc_display) > 12) {
  cc_display <- utils::head(cc_display, 12)
}
pca_display <- pca_summary[pca_summary$pc %in% c("1", "2", "3", 1, 2, 3), , drop = FALSE]
if (nrow(pca_display) > 12) {
  pca_display <- utils::head(pca_display, 12)
}

report_lines <- build_report_lines_v04(
  title = "09b Trajectory Input EDA",
  header_bullets = c(
    sprintf("Trajectory input rows: %s ok, %s not-ok.", ok_n, failed_n),
    sprintf("Outlier cells removed before trajectory modeling: %s.", removed_n),
    sprintf("Split balance warning threshold: %.0f%%; label minimum: %s cells.", cfg$trajectory_balance_warn_fraction * 100, cfg$trajectory_split_min_cells)
  ),
  key_files = list(
    input_status = cfg$trajectory_input_status_tsv,
    outlier_pre_post = cfg$trajectory_outlier_pre_post_summary_tsv,
    split_balance = cfg$trajectory_split_balance_summary_tsv,
    triage = cfg$trajectory_inputs_triage_tsv,
    outlier_violin = cfg$trajectory_inputs_qc_violin_png,
    umap_index = cfg$trajectory_umap_figure_index_tsv
  ),
  review_focus = c(
    "Approve trajectory_inputs only after every required pair/split has status ok or an accepted skipped/failed reason.",
    "For split trajectories, confirm syf/f5 balance and low-cell labels before interpreting downstream pseudotime."
  ),
  extra_sections = list(
    "Input Status Counts" = render_markdown_table_local(status_counts),
    "Input Status" = render_markdown_table_local(utils::head(input_status, 20)),
    "Split Balance Warnings" = render_markdown_table_local(utils::head(warning_balance, 20)),
    "Cell Cycle Summary" = render_markdown_table_local(cc_display),
    "PCA Summary Preview" = render_markdown_table_local(pca_display)
  ),
  triage_df = triage_df
)
ensure_dir(dirname(cfg$trajectory_inputs_report_md))
write_markdown_local(report_lines, cfg$trajectory_inputs_report_md)

if (file.exists(cfg$module_09b_manifest_path)) {
  unlink(cfg$module_09b_manifest_path)
}
outputs <- list(
  report_md = build_output_entry(cfg$trajectory_inputs_report_md, "md", module_name, "trajectory input EDA gate report", base_dir = cfg$project_root),
  input_status_tsv = build_output_entry(cfg$trajectory_input_status_tsv, "tsv", module_name, "trajectory input status table", base_dir = cfg$project_root, schema = infer_schema_from_df(input_status)),
  outlier_pre_post_summary_tsv = build_output_entry(cfg$trajectory_outlier_pre_post_summary_tsv, "tsv", module_name, "pre/post outlier metric summary", base_dir = cfg$project_root, schema = infer_schema_from_df(outlier_pre_post)),
  split_balance_summary_tsv = build_output_entry(cfg$trajectory_split_balance_summary_tsv, "tsv", module_name, "split balance diagnostics by pair/split/label", base_dir = cfg$project_root, schema = infer_schema_from_df(split_balance)),
  triage_tsv = build_output_entry(cfg$trajectory_inputs_triage_tsv, "tsv", module_name, "trajectory input EDA triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df)),
  outlier_pre_post_violin_png = build_output_entry(cfg$trajectory_inputs_qc_violin_png, "png", module_name, "pre/post outlier QC violin plot", base_dir = cfg$project_root),
  umap_figure_index_tsv = build_output_entry(cfg$trajectory_umap_figure_index_tsv, "tsv", module_name, "UMAP figure index for trajectory inputs", base_dir = cfg$project_root, schema = infer_schema_from_df(umap_index))
)
write_manifest_local(
  manifest_path = cfg$module_09b_manifest_path,
  new_outputs = outputs,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_09a = cfg$module_09a_manifest_path,
    trajectory_pairs = cfg$trajectory_pairs_sheet
  ),
  version = cfg$module_version,
  depends_on = list(module_09a = cfg$module_09a_manifest_path)
)

message("09b completed. report: ", cfg$trajectory_inputs_report_md)
