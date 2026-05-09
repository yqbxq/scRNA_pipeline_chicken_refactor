#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_10.R"), encoding = "UTF-8")

load_required_packages(c("dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_10()
module_name <- "10g_velocity_consistency"
prepare_dirs_10(cfg)

empty_metric_10g <- function(pair_id, split_value, metric_id, status, reason) {
  data.frame(
    pair_id = pair_id,
    split_value = display_scalar_value(split_value, "pooled"),
    metric_id = metric_id,
    group_a = "",
    group_b = "",
    n = 0L,
    value = NA_real_,
    p_value = NA_real_,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

cosine_10g <- function(x1, y1, x2, y2) {
  x1 <- suppressWarnings(as.numeric(x1))
  y1 <- suppressWarnings(as.numeric(y1))
  x2 <- suppressWarnings(as.numeric(x2))
  y2 <- suppressWarnings(as.numeric(y2))
  denom <- sqrt(x1^2 + y1^2) * sqrt(x2^2 + y2^2)
  out <- (x1 * x2 + y1 * y2) / denom
  out[!is.finite(out)] <- NA_real_
  out
}

spearman_metric_10g <- function(x, y) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) {
    return(list(n = sum(keep), value = NA_real_, p_value = NA_real_, status = "not_evaluated", reason = "fewer than 3 finite pairs"))
  }
  test <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  list(n = sum(keep), value = unname(test$estimate), p_value = test$p.value, status = "ok", reason = "")
}

metric_row_10g <- function(pair_id, split_value, metric_id, value, n, status = "ok", reason = "", p_value = NA_real_, group_a = "", group_b = "") {
  data.frame(
    pair_id = pair_id,
    split_value = display_scalar_value(split_value, "pooled"),
    metric_id = metric_id,
    group_a = group_a,
    group_b = group_b,
    n = n,
    value = value,
    p_value = p_value,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

read_delim_optional_10g <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read_tsv_optional(path)
  }
}

read_slingshot_index_10g <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09d_manifest_path,
    "slingshot_index_tsv",
    cfg$trajectory_slingshot_index_tsv,
    cols = c(trajectory_method_index_cols_09, "fine_output_path", "lineage_count", "root_label", "label_var")
  )
}

slingshot_row_for_velocity_10g <- function(slingshot_index, pair_id, split_value) {
  split_value <- display_scalar_value(split_value, "pooled")
  ok <- slingshot_index[
    slingshot_index$status == "ok" &
      slingshot_index$pair_id == pair_id &
      slingshot_index$split_value == split_value,
    ,
    drop = FALSE
  ]
  if (nrow(ok) == 0 && identical(split_value, "pooled")) {
    ok <- slingshot_index[
      slingshot_index$status == "ok" &
        slingshot_index$pair_id == pair_id,
      ,
      drop = FALSE
    ]
  }
  if (nrow(ok) == 0) {
    return(NULL)
  }
  ok[1, , drop = FALSE]
}

plot_consistency_10g <- function(metrics, path, title) {
  ensure_dir(dirname(path))
  ok <- metrics[metrics$status == "ok" & is.finite(suppressWarnings(as.numeric(metrics$value))), , drop = FALSE]
  if (nrow(ok) == 0) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No velocity consistency metrics") +
      ggplot2::theme_void()
  } else {
    ok$value <- suppressWarnings(as.numeric(ok$value))
    plot_obj <- ggplot2::ggplot(ok, ggplot2::aes(x = metric_id, y = value)) +
      ggplot2::geom_col(fill = "#4C78A8") +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = title, x = NULL, y = "metric value")
  }
  save_plot_local(plot_obj, path, width = 7, height = 4.8)
}

scvelo_rows <- velocity_read_scvelo_index_10(cfg)
steady_index <- read_tsv_optional(cfg$velocity_velocyto_steady_index_tsv)
cellrank_index <- velocity_read_cellrank_index_10(cfg, include_not_ok = TRUE)
driver_overlap <- read_tsv_optional(cfg$velocity_driver_overlap_tsv)
slingshot_index <- read_slingshot_index_10g(cfg)

index_rows <- list()
metric_rows <- list()
split_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(scvelo_rows))) {
  unit <- scvelo_rows[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  out_tsv <- velocity_consistency_path_10(cfg, pair_id, split_value)
  fig_png <- velocity_consistency_figure_path_10(cfg, pair_id, split_value)
  started <- proc.time()[["elapsed"]]
  rows <- list()

  vectors <- read_delim_optional_10g(unit$vector_tsv[[1]])
  if (nrow(vectors) == 0) {
    rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "dynamical_vs_stochastic_cosine", "not_evaluated", "missing 10c vector_tsv")
  } else {
    required <- c("velocity_umap_1", "velocity_umap_2", "velocity_umap_1_stochastic", "velocity_umap_2_stochastic")
    if (all(required %in% colnames(vectors))) {
      cos <- cosine_10g(vectors$velocity_umap_1, vectors$velocity_umap_2, vectors$velocity_umap_1_stochastic, vectors$velocity_umap_2_stochastic)
      finite_cos <- is.finite(cos)
      rows[[length(rows) + 1L]] <- metric_row_10g(
        pair_id,
        split_value,
        "dynamical_vs_stochastic_cosine",
        if (any(finite_cos)) mean(cos[finite_cos]) else NA_real_,
        sum(finite_cos),
        status = ifelse(any(finite_cos), "ok", "not_evaluated"),
        reason = ifelse(any(finite_cos), "", "no finite cosine values")
      )
    } else {
      rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "dynamical_vs_stochastic_cosine", "not_evaluated", "vector_tsv lacks stochastic velocity columns")
    }

    sling_row <- slingshot_row_for_velocity_10g(slingshot_index, pair_id, split_value)
    if (!is.null(sling_row) && file.exists(sling_row$output_path[[1]]) && "latent_time" %in% colnames(vectors)) {
      sling <- read_delim_optional_10g(sling_row$output_path[[1]])
      if (all(c("cell_id", "pseudotime") %in% colnames(sling))) {
        merged <- merge(vectors[, c("cell_id", "latent_time"), drop = FALSE], sling[, c("cell_id", "pseudotime"), drop = FALSE], by = "cell_id")
        cor <- spearman_metric_10g(merged$latent_time, merged$pseudotime)
        rows[[length(rows) + 1L]] <- metric_row_10g(pair_id, split_value, "latent_time_vs_slingshot_spearman", cor$value, cor$n, cor$status, cor$reason, cor$p_value)
      } else {
        rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "latent_time_vs_slingshot_spearman", "not_evaluated", "Slingshot pseudotime table lacks cell_id/pseudotime")
      }
    } else {
      rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "latent_time_vs_slingshot_spearman", "not_evaluated", "missing Slingshot match or latent_time")
    }
  }

  steady_hit <- steady_index[steady_index$pair_id == pair_id & steady_index$split_value == split_value & steady_index$status == "ok", , drop = FALSE]
  if (nrow(steady_hit) > 0 && nrow(vectors) > 0) {
    steady <- read_delim_optional_10g(steady_hit$direction_tsv[[1]])
    if (all(c("cell_id", "velocity_length") %in% colnames(steady)) && "velocity_length" %in% colnames(vectors)) {
      merged <- merge(vectors[, c("cell_id", "velocity_length"), drop = FALSE], steady[, c("cell_id", "velocity_length"), drop = FALSE], by = "cell_id", suffixes = c("_scvelo", "_velocyto"))
      cor <- spearman_metric_10g(merged$velocity_length_scvelo, merged$velocity_length_velocyto)
      rows[[length(rows) + 1L]] <- metric_row_10g(pair_id, split_value, "dynamical_vs_velocyto_length_spearman", cor$value, cor$n, cor$status, cor$reason, cor$p_value)
    } else {
      rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "dynamical_vs_velocyto_length_spearman", "not_evaluated", "velocyto direction table lacks length values")
    }
  } else {
    rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "dynamical_vs_velocyto_length_spearman", "not_evaluated", "missing successful 10d steady-state row")
  }

  cr_hit <- cellrank_index[cellrank_index$pair_id == pair_id & cellrank_index$split_value == split_value & cellrank_index$status == "ok", , drop = FALSE]
  if (nrow(cr_hit) > 0 && nrow(vectors) > 0) {
    fate <- read_delim_optional_10g(cr_hit$fate_csv[[1]])
    sling_row <- slingshot_row_for_velocity_10g(slingshot_index, pair_id, split_value)
    if (nrow(fate) > 0 && !is.null(sling_row) && file.exists(sling_row$output_path[[1]])) {
      fate_cols <- setdiff(colnames(fate), c("pair_id", "split_value", "cell_id"))
      fate$max_fate <- if (length(fate_cols) > 0) fate_cols[max.col(as.matrix(fate[, fate_cols, drop = FALSE]), ties.method = "first")] else ""
      sling <- read_delim_optional_10g(sling_row$output_path[[1]])
      if (all(c("cell_id", "label") %in% colnames(sling))) {
        conf <- merge(fate[, c("cell_id", "max_fate"), drop = FALSE], sling[, c("cell_id", "label"), drop = FALSE], by = "cell_id")
        if (nrow(conf) > 0) {
          tab <- as.data.frame(table(conf$label, conf$max_fate), stringsAsFactors = FALSE)
          for (j in seq_len(nrow(tab))) {
            rows[[length(rows) + 1L]] <- metric_row_10g(pair_id, split_value, "cellrank_fate_vs_slingshot_lineage_count", tab$Freq[[j]], tab$Freq[[j]], "ok", "", group_a = as.character(tab$Var1[[j]]), group_b = as.character(tab$Var2[[j]]))
          }
        } else {
          rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "cellrank_fate_vs_slingshot_lineage_count", "not_evaluated", "no overlapping CellRank/Slingshot cells")
        }
      }
    } else {
      rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "cellrank_fate_vs_slingshot_lineage_count", "not_evaluated", "missing CellRank fate or Slingshot table")
    }
  } else {
    rows[[length(rows) + 1L]] <- empty_metric_10g(pair_id, split_value, "cellrank_fate_vs_slingshot_lineage_count", "not_evaluated", "missing successful 10f CellRank row")
  }

  metrics <- dplyr::bind_rows(rows)
  write_tsv_local(metrics, out_tsv)
  plot_consistency_10g(metrics, fig_png, sprintf("Velocity consistency %s %s", pair_id, split_value))
  metric_rows[[length(metric_rows) + 1L]] <- metrics
  index_rows[[length(index_rows) + 1L]] <- data.frame(
    pair_id = pair_id,
    split_value = split_value,
    method = "velocity_consistency",
    methods_enabled = "yes",
    input_path = unit$output_path[[1]],
    output_path = out_tsv,
    extra_path = "",
    figure_path = fig_png,
    n_cells = suppressWarnings(as.integer(unit$n_cells[[1]])),
    status = "ok",
    reason = "",
    runtime_s = round(proc.time()[["elapsed"]] - started, 3),
    split_compare_path = "",
    stringsAsFactors = FALSE
  )
  dynamic_outputs[[sprintf("velocity_consistency__%s", velocity_unit_file_id_10(pair_id, split_value))]] <- build_output_entry(out_tsv, "tsv", module_name, "velocity consistency metrics for one unit", base_dir = cfg$project_root)
}

for (pair_id in unique(scvelo_rows$pair_id)) {
  pair_units <- scvelo_rows[scvelo_rows$pair_id == pair_id & scvelo_rows$split_value != "pooled", , drop = FALSE]
  if (nrow(pair_units) < 2L) {
    next
  }
  split_out <- velocity_split_compare_path_10(cfg, pair_id)
  split_fig <- velocity_split_compare_figure_path_10(cfg, pair_id)
  rows <- list()
  vecs <- lapply(seq_len(nrow(pair_units)), function(i) read_delim_optional_10g(pair_units$vector_tsv[[i]]))
  names(vecs) <- pair_units$split_value
  split_names <- names(vecs)
  for (a in seq_along(split_names)) {
    for (b in seq_along(split_names)) {
      if (b <= a) next
      va <- vecs[[a]]
      vb <- vecs[[b]]
      if (nrow(va) == 0 || nrow(vb) == 0 || !"cluster" %in% colnames(va) || !"cluster" %in% colnames(vb)) next
      for (cluster in intersect(unique(va$cluster), unique(vb$cluster))) {
        ca <- va[va$cluster == cluster, , drop = FALSE]
        cb <- vb[vb$cluster == cluster, , drop = FALSE]
        angle_a <- atan2(suppressWarnings(as.numeric(ca$velocity_umap_2)), suppressWarnings(as.numeric(ca$velocity_umap_1)))
        angle_b <- atan2(suppressWarnings(as.numeric(cb$velocity_umap_2)), suppressWarnings(as.numeric(cb$velocity_umap_1)))
        rows[[length(rows) + 1L]] <- data.frame(
          pair_id = pair_id,
          comparison_id = sprintf("%s_vs_%s", split_names[[a]], split_names[[b]]),
          metric_id = "cluster_velocity_angle_delta",
          group_a = cluster,
          group_b = cluster,
          n_a = sum(is.finite(angle_a)),
          n_b = sum(is.finite(angle_b)),
          value = mean(angle_a, na.rm = TRUE) - mean(angle_b, na.rm = TRUE),
          p_value = NA_real_,
          status = "ok",
          reason = "",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (nrow(driver_overlap) > 0) {
    hit <- driver_overlap[driver_overlap$pair_id == pair_id, , drop = FALSE]
    for (j in seq_len(nrow(hit))) {
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = pair_id,
        comparison_id = sprintf("%s_vs_%s", hit$split_a[[j]], hit$split_b[[j]]),
        metric_id = "driver_overlap_jaccard",
        group_a = hit$split_a[[j]],
        group_b = hit$split_b[[j]],
        n_a = suppressWarnings(as.integer(hit$driver_n_a[[j]])),
        n_b = suppressWarnings(as.integer(hit$driver_n_b[[j]])),
        value = suppressWarnings(as.numeric(hit$jaccard[[j]])),
        p_value = NA_real_,
        status = "ok",
        reason = "",
        stringsAsFactors = FALSE
      )
    }
  }
  split_df <- if (length(rows) > 0) dplyr::bind_rows(rows) else velocity_empty_df_10(velocity_split_compare_cols_10)
  write_tsv_local(split_df, split_out)
  plot_consistency_10g(transform(split_df, split_value = "split"), split_fig, sprintf("Velocity split comparison %s", pair_id))
  split_rows[[length(split_rows) + 1L]] <- split_df
  dynamic_outputs[[sprintf("velocity_split_compare__%s", safe_id_09(pair_id))]] <- build_output_entry(split_out, "tsv", module_name, "velocity split comparison metrics", base_dir = cfg$project_root)
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else velocity_empty_df_10(velocity_consistency_index_cols_10)
summary_df <- if (length(metric_rows) > 0) dplyr::bind_rows(metric_rows) else velocity_empty_df_10(velocity_consistency_summary_cols_10)
split_df <- if (length(split_rows) > 0) dplyr::bind_rows(split_rows) else velocity_empty_df_10(velocity_split_compare_cols_10)

write_tsv_local(index_df, cfg$velocity_consistency_index_tsv)
write_tsv_local(summary_df, cfg$velocity_consistency_summary_tsv)
write_tsv_local(split_df, cfg$velocity_split_compare_tsv)

velocity_manifest_10(
  cfg,
  cfg$module_10g_manifest_path,
  module_name,
  outputs = c(
    list(
      velocity_consistency_index_tsv = build_output_entry(cfg$velocity_consistency_index_tsv, "tsv", module_name, "velocity consistency status by unit", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
      velocity_consistency_summary_tsv = build_output_entry(cfg$velocity_consistency_summary_tsv, "tsv", module_name, "velocity consistency metrics by unit", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
      velocity_split_compare_tsv = build_output_entry(cfg$velocity_split_compare_tsv, "tsv", module_name, "velocity split comparison metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(split_df))
    ),
    dynamic_outputs
  ),
  inputs = list(
    module_10c = cfg$module_10c_manifest_path,
    module_10d = cfg$module_10d_manifest_path,
    module_10e = cfg$module_10e_manifest_path,
    module_10f = cfg$module_10f_manifest_path,
    module_09d = cfg$module_09d_manifest_path
  ),
  depends_on = list(
    module_10c = cfg$module_10c_manifest_path,
    module_10d = cfg$module_10d_manifest_path,
    module_10e = cfg$module_10e_manifest_path,
    module_10f = cfg$module_10f_manifest_path,
    module_09d = cfg$module_09d_manifest_path
  )
)

message("10g completed. velocity consistency index: ", cfg$velocity_consistency_index_tsv)
