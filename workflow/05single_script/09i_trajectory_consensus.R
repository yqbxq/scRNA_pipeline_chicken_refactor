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
source_utf8(file.path(.script_dir, "helpers", "trajectory_utils.R"))

load_required_packages(c("dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09i_trajectory_consensus"
prepare_dirs_09(cfg)

consensus_index_cols <- c(
  trajectory_method_index_cols_09,
  "method_count", "methods_used", "correlation_path", "jaccard_path",
  "consensus_pseudotime_path", "conflict_count"
)
correlation_cols <- c("pair_id", "split_value", "method_a", "method_b", "spearman_rho", "kendall_tau", "p_value", "n_cells")
jaccard_cols <- c("pair_id", "split_value", "method_a", "method_b", "label_jaccard", "n_shared_labels")
split_agreement_cols <- c("pair_id", "split_a", "split_b", "method_count_a", "method_count_b", "ks_statistic", "ks_p_value", "median_consensus_a", "median_consensus_b", "status", "reason")

method_specs_09i <- list(
  slingshot = list(manifest = cfg$module_09d_manifest_path, key = "slingshot_index_tsv", fallback = cfg$trajectory_slingshot_index_tsv),
  monocle3 = list(manifest = cfg$module_09e_manifest_path, key = "monocle3_index_tsv", fallback = cfg$trajectory_monocle3_index_tsv),
  paga_dpt = list(manifest = cfg$module_09f_manifest_path, key = "paga_dpt_index_tsv", fallback = cfg$trajectory_paga_index_tsv),
  palantir = list(manifest = cfg$module_09g_manifest_path, key = "palantir_index_tsv", fallback = cfg$trajectory_palantir_index_tsv),
  monocle2 = list(manifest = cfg$module_09e2_manifest_path, key = "monocle2_index_tsv", fallback = cfg$trajectory_monocle2_index_tsv)
)

read_method_indexes_09i <- function() {
  out <- list()
  for (method in names(method_specs_09i)) {
    spec <- method_specs_09i[[method]]
    out[[method]] <- trajectory_read_method_index_09(spec$manifest, spec$key, spec$fallback)
  }
  out
}

method_row_for_unit_09i <- function(index_df, pair_id, split_value) {
  if (nrow(index_df) == 0) {
    return(NULL)
  }
  split_value <- display_scalar_value(split_value, "pooled")
  hit <- index_df[index_df$pair_id == pair_id & index_df$split_value == split_value, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(NULL)
  }
  hit[1, , drop = FALSE]
}

read_unit_method_pt_09i <- function(indexes, pair_id, split_value) {
  rows <- list()
  for (method in names(indexes)) {
    row <- method_row_for_unit_09i(indexes[[method]], pair_id, split_value)
    if (is.null(row) || row$status[[1]] != "ok" || !file.exists(row$output_path[[1]])) {
      next
    }
    pt <- trajectory_read_pseudotime_csv_09(row$output_path[[1]], method)
    if (nrow(pt) == 0) {
      next
    }
    pt$pseudotime <- trajectory_scale01_09(pt$pseudotime)
    rows[[method]] <- pt
  }
  rows
}

label_jaccard_09i <- function(a, b) {
  if (!"label" %in% colnames(a) || !"label" %in% colnames(b)) {
    return(c(label_jaccard = NA_real_, n_shared_labels = 0))
  }
  labels <- intersect(unique(as.character(a$label)), unique(as.character(b$label)))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0) {
    return(c(label_jaccard = NA_real_, n_shared_labels = 0))
  }
  scores <- vapply(labels, function(label) {
    a_cells <- a$cell_id[as.character(a$label) == label]
    b_cells <- b$cell_id[as.character(b$label) == label]
    denom <- length(union(a_cells, b_cells))
    if (denom == 0) NA_real_ else length(intersect(a_cells, b_cells)) / denom
  }, numeric(1))
  c(label_jaccard = mean(scores, na.rm = TRUE), n_shared_labels = length(labels))
}

consensus_for_unit_09i <- function(unit, indexes) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  unit_id <- trajectory_unit_file_id_09(pair_id, split_value)
  corr_tsv <- trajectory_method_output_path_09(cfg, "consensus", pair_id, split_value, "consensus_pseudotime_correlation", "tsv")
  jaccard_tsv <- trajectory_method_output_path_09(cfg, "consensus", pair_id, split_value, "consensus_jaccard", "tsv")
  consensus_csv <- trajectory_method_output_path_09(cfg, "consensus", pair_id, split_value, "consensus_pseudotime", "csv")
  fig_png <- trajectory_method_figure_path_09(cfg, "consensus", pair_id, split_value, "Figure_ConsensusCorrelation")
  ensure_dir(dirname(corr_tsv))
  ensure_dir(dirname(fig_png))
  started <- proc.time()[["elapsed"]]

  pt_list <- read_unit_method_pt_09i(indexes, pair_id, split_value)
  methods <- names(pt_list)
  if (length(methods) < cfg$trajectory_consensus_min_methods) {
    write_tsv_local(trajectory_empty_df_09(correlation_cols), corr_tsv)
    write_tsv_local(trajectory_empty_df_09(jaccard_cols), jaccard_tsv)
    write.csv(data.frame(cell_id = character(0), consensus_pseudotime = numeric(0)), consensus_csv, row.names = FALSE)
    return(list(
      index = data.frame(
        pair_id = pair_id, split_value = split_value, method = "consensus", methods_enabled = "yes",
        input_rds = unit$input_rds[[1]], output_path = corr_tsv, extra_path = jaccard_tsv,
        figure_path = fig_png, n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
        status = "skipped_insufficient_methods",
        reason = sprintf("need at least %s successful methods; observed %s", cfg$trajectory_consensus_min_methods, length(methods)),
        runtime_s = round(proc.time()[["elapsed"]] - started, 3), method_count = length(methods),
        methods_used = paste(methods, collapse = ","), correlation_path = corr_tsv,
        jaccard_path = jaccard_tsv, consensus_pseudotime_path = consensus_csv,
        conflict_count = 0L, stringsAsFactors = FALSE
      ),
      correlations = trajectory_empty_df_09(correlation_cols),
      jaccards = trajectory_empty_df_09(jaccard_cols),
      consensus = data.frame(cell_id = character(0), consensus_pseudotime = numeric(0), stringsAsFactors = FALSE)
    ))
  }

  common <- Reduce(intersect, lapply(pt_list, function(x) x$cell_id[is.finite(x$pseudotime)]))
  if (length(common) < 3L) {
    write_tsv_local(trajectory_empty_df_09(correlation_cols), corr_tsv)
    write_tsv_local(trajectory_empty_df_09(jaccard_cols), jaccard_tsv)
    write.csv(data.frame(cell_id = character(0), consensus_pseudotime = numeric(0)), consensus_csv, row.names = FALSE)
    return(list(
      index = data.frame(
        pair_id = pair_id, split_value = split_value, method = "consensus", methods_enabled = "yes",
        input_rds = unit$input_rds[[1]], output_path = corr_tsv, extra_path = jaccard_tsv,
        figure_path = fig_png, n_cells = length(common), status = "skipped_low_common_cells",
        reason = "fewer than 3 common finite cells across methods",
        runtime_s = round(proc.time()[["elapsed"]] - started, 3), method_count = length(methods),
        methods_used = paste(methods, collapse = ","), correlation_path = corr_tsv,
        jaccard_path = jaccard_tsv, consensus_pseudotime_path = consensus_csv,
        conflict_count = 0L, stringsAsFactors = FALSE
      ),
      correlations = trajectory_empty_df_09(correlation_cols),
      jaccards = trajectory_empty_df_09(jaccard_cols),
      consensus = data.frame(cell_id = character(0), consensus_pseudotime = numeric(0), stringsAsFactors = FALSE)
    ))
  }

  mat <- do.call(cbind, lapply(pt_list, function(df) {
    values <- df$pseudotime[match(common, df$cell_id)]
    trajectory_scale01_09(values)
  }))
  colnames(mat) <- methods
  rownames(mat) <- common
  pairs <- utils::combn(methods, 2, simplify = FALSE)
  corr_rows <- lapply(pairs, function(pair) {
    x <- mat[, pair[[1]]]
    y <- mat[, pair[[2]]]
    sp <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
    kd <- suppressWarnings(stats::cor.test(x, y, method = "kendall", exact = FALSE))
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method_a = pair[[1]],
      method_b = pair[[2]],
      spearman_rho = unname(sp$estimate),
      kendall_tau = unname(kd$estimate),
      p_value = sp$p.value,
      n_cells = length(common),
      stringsAsFactors = FALSE
    )
  })
  corr_df <- dplyr::bind_rows(corr_rows)

  jaccard_rows <- lapply(pairs, function(pair) {
    score <- label_jaccard_09i(pt_list[[pair[[1]]]], pt_list[[pair[[2]]]])
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method_a = pair[[1]],
      method_b = pair[[2]],
      label_jaccard = score[["label_jaccard"]],
      n_shared_labels = score[["n_shared_labels"]],
      stringsAsFactors = FALSE
    )
  })
  jaccard_df <- dplyr::bind_rows(jaccard_rows)
  consensus_df <- data.frame(
    cell_id = common,
    consensus_pseudotime = rowMeans(mat, na.rm = TRUE),
    method_count = ncol(mat),
    stringsAsFactors = FALSE
  )
  write_tsv_local(corr_df, corr_tsv)
  write_tsv_local(jaccard_df, jaccard_tsv)
  write.csv(consensus_df, consensus_csv, row.names = FALSE)

  plot_df <- corr_df
  plot_df$method_pair <- paste(plot_df$method_a, plot_df$method_b, sep = " vs ")
  plot_obj <- ggplot2::ggplot(plot_df, ggplot2::aes(x = method_pair, y = spearman_rho, fill = spearman_rho < cfg$trajectory_consensus_conflict_rho)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::scale_fill_manual(values = c("FALSE" = "#4C78A8", "TRUE" = "#D95F02"), guide = "none") +
    ggplot2::labs(title = sprintf("Consensus correlation %s %s", pair_id, split_value), x = NULL, y = "Spearman rho")
  save_plot_local(plot_obj, fig_png, width = 6.5, height = 4.5)

  conflict_n <- sum(corr_df$spearman_rho < cfg$trajectory_consensus_conflict_rho, na.rm = TRUE)
  list(
    index = data.frame(
      pair_id = pair_id, split_value = split_value, method = "consensus", methods_enabled = "yes",
      input_rds = unit$input_rds[[1]], output_path = corr_tsv, extra_path = jaccard_tsv,
      figure_path = fig_png, n_cells = length(common), status = "ok", reason = "",
      runtime_s = round(proc.time()[["elapsed"]] - started, 3), method_count = length(methods),
      methods_used = paste(methods, collapse = ","), correlation_path = corr_tsv,
      jaccard_path = jaccard_tsv, consensus_pseudotime_path = consensus_csv,
      conflict_count = conflict_n, stringsAsFactors = FALSE
    ),
    correlations = corr_df,
    jaccards = jaccard_df,
    consensus = consensus_df
  )
}

split_agreement_09i <- function(index_df) {
  ok <- index_df[index_df$status == "ok" & nzchar(index_df$consensus_pseudotime_path), , drop = FALSE]
  rows <- list()
  for (pair_id in unique(ok$pair_id)) {
    pair_rows <- ok[ok$pair_id == pair_id, , drop = FALSE]
    split_values <- setdiff(unique(pair_rows$split_value), "pooled")
    if (length(split_values) < 2) {
      next
    }
    for (pair in utils::combn(split_values, 2, simplify = FALSE)) {
      a <- pair_rows[pair_rows$split_value == pair[[1]], , drop = FALSE][1, , drop = FALSE]
      b <- pair_rows[pair_rows$split_value == pair[[2]], , drop = FALSE][1, , drop = FALSE]
      pa <- read.csv(a$consensus_pseudotime_path[[1]], stringsAsFactors = FALSE)
      pb <- read.csv(b$consensus_pseudotime_path[[1]], stringsAsFactors = FALSE)
      if (nrow(pa) < 3 || nrow(pb) < 3) {
        rows[[length(rows) + 1L]] <- data.frame(pair_id = pair_id, split_a = pair[[1]], split_b = pair[[2]], method_count_a = a$method_count[[1]], method_count_b = b$method_count[[1]], ks_statistic = NA_real_, ks_p_value = NA_real_, median_consensus_a = NA_real_, median_consensus_b = NA_real_, status = "skipped_low_cells", reason = "too few consensus cells", stringsAsFactors = FALSE)
        next
      }
      ks <- suppressWarnings(stats::ks.test(pa$consensus_pseudotime, pb$consensus_pseudotime))
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = pair_id,
        split_a = pair[[1]],
        split_b = pair[[2]],
        method_count_a = a$method_count[[1]],
        method_count_b = b$method_count[[1]],
        ks_statistic = unname(ks$statistic),
        ks_p_value = ks$p.value,
        median_consensus_a = median(pa$consensus_pseudotime, na.rm = TRUE),
        median_consensus_b = median(pb$consensus_pseudotime, na.rm = TRUE),
        status = "ok",
        reason = "",
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) trajectory_empty_df_09(split_agreement_cols) else dplyr::bind_rows(rows)
}

units <- trajectory_execution_units_09(cfg)
indexes <- read_method_indexes_09i()
index_rows <- list()
corr_rows <- list()
jaccard_rows <- list()
triage_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  result <- consensus_for_unit_09i(units[i, , drop = FALSE], indexes)
  index_rows[[length(index_rows) + 1L]] <- result$index
  corr_rows[[length(corr_rows) + 1L]] <- result$correlations
  jaccard_rows[[length(jaccard_rows) + 1L]] <- result$jaccards
  unit_id <- trajectory_unit_file_id_09(result$index$pair_id[[1]], result$index$split_value[[1]])
  dynamic_outputs[[sprintf("consensus_correlation_%s", unit_id)]] <- build_output_entry(result$index$correlation_path[[1]], "tsv", module_name, "pairwise method pseudotime correlations", base_dir = cfg$project_root)
  dynamic_outputs[[sprintf("consensus_pseudotime_%s", unit_id)]] <- build_output_entry(result$index$consensus_pseudotime_path[[1]], "csv", module_name, "mean scaled consensus pseudotime", base_dir = cfg$project_root)
  if (nrow(result$correlations) > 0) {
    conflicts <- result$correlations[result$correlations$spearman_rho < cfg$trajectory_consensus_conflict_rho, , drop = FALSE]
    for (j in seq_len(nrow(conflicts))) {
      triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
        sample_id = conflicts$pair_id[[j]],
        severity = "warning",
        signal_id = "trajectory_consensus_low_correlation",
        evidence = sprintf("%s split=%s %s vs %s Spearman rho=%.3f", conflicts$pair_id[[j]], conflicts$split_value[[j]], conflicts$method_a[[j]], conflicts$method_b[[j]], conflicts$spearman_rho[[j]]),
        recommended_action = "审阅各方法 pseudotime UMAP 和 root 方向；必要时调整 root 或将冲突方法降级为探索。"
      )
    }
  }
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(consensus_index_cols)
corr_df <- if (length(corr_rows) > 0) dplyr::bind_rows(corr_rows) else trajectory_empty_df_09(correlation_cols)
jaccard_df <- if (length(jaccard_rows) > 0) dplyr::bind_rows(jaccard_rows) else trajectory_empty_df_09(jaccard_cols)
split_agreement_df <- split_agreement_09i(index_df)
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

write_tsv_local(index_df, cfg$trajectory_consensus_index_tsv)
write_tsv_local(corr_df, cfg$trajectory_consensus_correlation_summary_tsv)
write_tsv_local(jaccard_df, cfg$trajectory_consensus_jaccard_summary_tsv)
write_tsv_local(split_agreement_df, cfg$trajectory_consensus_split_agreement_tsv)
write_tsv_local(triage_df, cfg$trajectory_consensus_triage_tsv)

outputs <- c(
  list(
    consensus_index_tsv = build_output_entry(cfg$trajectory_consensus_index_tsv, "tsv", module_name, "consensus status by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
    consensus_correlation_summary_tsv = build_output_entry(cfg$trajectory_consensus_correlation_summary_tsv, "tsv", module_name, "all pairwise method pseudotime correlations", base_dir = cfg$project_root, schema = infer_schema_from_df(corr_df)),
    consensus_jaccard_summary_tsv = build_output_entry(cfg$trajectory_consensus_jaccard_summary_tsv, "tsv", module_name, "all pairwise label Jaccard summaries", base_dir = cfg$project_root, schema = infer_schema_from_df(jaccard_df)),
    consensus_split_agreement_tsv = build_output_entry(cfg$trajectory_consensus_split_agreement_tsv, "tsv", module_name, "split-mode consensus distribution agreement", base_dir = cfg$project_root, schema = infer_schema_from_df(split_agreement_df)),
    consensus_triage_tsv = build_output_entry(cfg$trajectory_consensus_triage_tsv, "tsv", module_name, "consensus conflict triage rows", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
  ),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09i_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09d = cfg$module_09d_manifest_path, module_09e = cfg$module_09e_manifest_path, module_09f = cfg$module_09f_manifest_path, module_09g = cfg$module_09g_manifest_path, module_09e2 = cfg$module_09e2_manifest_path),
  depends_on = list(module_09d = cfg$module_09d_manifest_path, module_09e = cfg$module_09e_manifest_path, module_09f = cfg$module_09f_manifest_path, module_09g = cfg$module_09g_manifest_path)
)

message("09i completed. consensus index: ", cfg$trajectory_consensus_index_tsv)
