#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")

load_required_packages(c("dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09m_trajectory_split_compare"
prepare_dirs_09(cfg)

pseudotime_cols_09m <- c(
  "pair_id", "split_a", "split_b", "label_var", "label_value",
  "n_cells_a", "n_cells_b", "median_pseudotime_a", "median_pseudotime_b",
  "wilcox_p_value", "ks_statistic", "ks_p_value", "status", "reason"
)

lineage_cols_09m <- c(
  "pair_id", "split_a", "split_b", "label_var", "lineage_count_a",
  "lineage_count_b", "label_n_a", "label_n_b", "label_jaccard",
  "mean_pseudotime_spearman", "status", "reason"
)

driver_cols_09m <- c(
  "pair_id", "split_a", "split_b", "option", "driver_n_a",
  "driver_n_b", "overlap_n", "jaccard", "top_overlap_genes",
  "status", "reason"
)

option_cols_09m <- c(
  "pair_id", "split_values", "option_a_driver_n", "option_b_driver_n",
  "overlap_n", "jaccard", "status", "reason"
)

entropy_cols_09m <- c(
  "pair_id", "split_a", "split_b", "label_var", "label_value",
  "n_cells_a", "n_cells_b", "median_entropy_a", "median_entropy_b",
  "wilcox_p_value", "ks_statistic", "ks_p_value", "status", "reason"
)

split_index_cols_09m <- c(
  "pair_id", "split_values", "pseudotime_distribution_path",
  "lineage_backbone_path", "driver_overlap_path",
  "option_a_vs_b_consistency_path", "entropy_compare_path",
  "figure_path", "report_path", "status", "reason"
)

tradeseq_index_cols_09m <- c(
  trajectory_method_index_cols_09,
  "option", "association_path", "pattern_path", "driver_path", "conditiontest_path"
)

slingshot_index_cols_09m <- c(
  trajectory_method_index_cols_09,
  "fine_output_path", "lineage_count", "root_label", "label_var"
)

palantir_index_cols_09m <- c(trajectory_method_index_cols_09, "fate_path")

read_slingshot_index_09m <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09d_manifest_path,
    "slingshot_index_tsv",
    cfg$trajectory_slingshot_index_tsv,
    cols = slingshot_index_cols_09m
  )
}

read_tradeseq_index_09m <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09h_manifest_path,
    "tradeseq_index_tsv",
    cfg$trajectory_tradeseq_index_tsv,
    cols = tradeseq_index_cols_09m
  )
}

read_palantir_index_09m <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09g_manifest_path,
    "palantir_index_tsv",
    cfg$trajectory_palantir_index_tsv,
    cols = palantir_index_cols_09m
  )
}

file_available_09m <- function(path) {
  path <- normalize_scalar_value(path)
  nzchar(path) && file.exists(path) && !isTRUE(file.info(path)$size == 0)
}

split_pair_ids_09m <- function(cfg) {
  pairs <- trajectory_active_pairs_09(cfg)
  if (nrow(pairs) == 0) {
    return(character(0))
  }
  if (!"condition_split_var" %in% colnames(pairs)) {
    pairs$condition_split_var <- ""
  }
  ids <- pairs$trajectory_id[nzchar(vapply(pairs$condition_split_var, normalize_scalar_value, character(1)))]
  unique(ids[nzchar(ids)])
}

row_for_unit_09m <- function(index_df, pair_id, split_value, option = NULL) {
  if (nrow(index_df) == 0) {
    return(NULL)
  }
  split_value <- display_scalar_value(split_value, "pooled")
  hit <- index_df[index_df$pair_id == pair_id & index_df$split_value == split_value, , drop = FALSE]
  if (!is.null(option) && "option" %in% colnames(hit)) {
    hit <- hit[hit$option == option, , drop = FALSE]
  }
  if (nrow(hit) == 0) NULL else hit[1, , drop = FALSE]
}

read_slingshot_pt_09m <- function(row) {
  if (is.null(row) || row$status[[1]] != "ok" || !file_available_09m(row$output_path[[1]])) {
    return(data.frame(cell_id = character(0), pseudotime = numeric(0), label = character(0), stringsAsFactors = FALSE))
  }
  df <- trajectory_read_pseudotime_csv_09(row$output_path[[1]], "slingshot")
  if (!"label" %in% colnames(df)) {
    df$label <- "__all__"
  }
  df$label <- as.character(df$label)
  df$label[is.na(df$label) | !nzchar(df$label)] <- "unknown"
  df$pseudotime <- trajectory_scale01_09(df$pseudotime)
  df[, c("cell_id", "pseudotime", "label"), drop = FALSE]
}

compare_numeric_09m <- function(a, b) {
  a <- suppressWarnings(as.numeric(a))
  b <- suppressWarnings(as.numeric(b))
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 3L || length(b) < 3L) {
    return(list(wilcox_p = NA_real_, ks_stat = NA_real_, ks_p = NA_real_, status = "skipped_low_cells", reason = "fewer than 3 finite values in at least one split"))
  }
  wx <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE))
  ks <- suppressWarnings(stats::ks.test(a, b))
  list(
    wilcox_p = wx$p.value,
    ks_stat = unname(ks$statistic),
    ks_p = ks$p.value,
    status = "ok",
    reason = ""
  )
}

compare_pseudotime_pair_09m <- function(pair_id, split_a, split_b, sling_a, sling_b) {
  pt_a <- read_slingshot_pt_09m(sling_a)
  pt_b <- read_slingshot_pt_09m(sling_b)
  label_var <- if (!is.null(sling_a)) normalize_scalar_value(sling_a$label_var[[1]], "label") else "label"
  if (nrow(pt_a) == 0 || nrow(pt_b) == 0) {
    return(data.frame(
      pair_id = pair_id, split_a = split_a, split_b = split_b, label_var = label_var,
      label_value = "__all__", n_cells_a = nrow(pt_a), n_cells_b = nrow(pt_b),
      median_pseudotime_a = NA_real_, median_pseudotime_b = NA_real_,
      wilcox_p_value = NA_real_, ks_statistic = NA_real_, ks_p_value = NA_real_,
      status = "skipped_no_slingshot", reason = "missing successful Slingshot pseudotime",
      stringsAsFactors = FALSE
    ))
  }
  labels <- intersect(unique(pt_a$label), unique(pt_b$label))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0) {
    labels <- "__all__"
    pt_a$label <- "__all__"
    pt_b$label <- "__all__"
  }
  rows <- lapply(labels, function(label) {
    va <- pt_a$pseudotime[pt_a$label == label]
    vb <- pt_b$pseudotime[pt_b$label == label]
    cmp <- compare_numeric_09m(va, vb)
    data.frame(
      pair_id = pair_id,
      split_a = split_a,
      split_b = split_b,
      label_var = label_var,
      label_value = label,
      n_cells_a = sum(is.finite(va)),
      n_cells_b = sum(is.finite(vb)),
      median_pseudotime_a = safe_median(va),
      median_pseudotime_b = safe_median(vb),
      wilcox_p_value = cmp$wilcox_p,
      ks_statistic = cmp$ks_stat,
      ks_p_value = cmp$ks_p,
      status = cmp$status,
      reason = cmp$reason,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

label_jaccard_09m <- function(a, b) {
  union_n <- length(union(a, b))
  if (union_n == 0L) {
    return(NA_real_)
  }
  length(intersect(a, b)) / union_n
}

compare_lineage_pair_09m <- function(pair_id, split_a, split_b, sling_a, sling_b) {
  pt_a <- read_slingshot_pt_09m(sling_a)
  pt_b <- read_slingshot_pt_09m(sling_b)
  label_var <- if (!is.null(sling_a)) normalize_scalar_value(sling_a$label_var[[1]], "label") else "label"
  if (is.null(sling_a) || is.null(sling_b) || nrow(pt_a) == 0 || nrow(pt_b) == 0) {
    return(data.frame(
      pair_id = pair_id, split_a = split_a, split_b = split_b, label_var = label_var,
      lineage_count_a = 0L, lineage_count_b = 0L, label_n_a = dplyr::n_distinct(pt_a$label),
      label_n_b = dplyr::n_distinct(pt_b$label), label_jaccard = NA_real_,
      mean_pseudotime_spearman = NA_real_, status = "skipped_no_slingshot",
      reason = "missing successful Slingshot pseudotime", stringsAsFactors = FALSE
    ))
  }
  labels_a <- unique(pt_a$label[nzchar(pt_a$label)])
  labels_b <- unique(pt_b$label[nzchar(pt_b$label)])
  common_labels <- intersect(labels_a, labels_b)
  rho <- NA_real_
  if (length(common_labels) >= 3L) {
    mean_a <- stats::aggregate(pseudotime ~ label, data = pt_a[pt_a$label %in% common_labels, , drop = FALSE], FUN = mean)
    mean_b <- stats::aggregate(pseudotime ~ label, data = pt_b[pt_b$label %in% common_labels, , drop = FALSE], FUN = mean)
    merged <- merge(mean_a, mean_b, by = "label", suffixes = c("_a", "_b"), sort = FALSE)
    if (nrow(merged) >= 3L) {
      rho <- suppressWarnings(stats::cor(merged$pseudotime_a, merged$pseudotime_b, method = "spearman", use = "complete.obs"))
    }
  }
  data.frame(
    pair_id = pair_id,
    split_a = split_a,
    split_b = split_b,
    label_var = label_var,
    lineage_count_a = suppressWarnings(as.integer(sling_a$lineage_count[[1]])),
    lineage_count_b = suppressWarnings(as.integer(sling_b$lineage_count[[1]])),
    label_n_a = length(labels_a),
    label_n_b = length(labels_b),
    label_jaccard = label_jaccard_09m(labels_a, labels_b),
    mean_pseudotime_spearman = rho,
    status = "ok",
    reason = "",
    stringsAsFactors = FALSE
  )
}

driver_genes_09m <- function(row) {
  if (is.null(row) || row$status[[1]] != "ok" || !file_available_09m(row$driver_path[[1]])) {
    return(character(0))
  }
  df <- read_tsv_optional(row$driver_path[[1]])
  gene_col <- intersect(c("gene_id", "gene", "gene_name", "symbol"), colnames(df))
  if (length(gene_col) == 0) {
    return(character(0))
  }
  genes <- trimws(as.character(df[[gene_col[[1]]]]))
  unique(genes[!is.na(genes) & nzchar(genes)])
}

compare_driver_pair_09m <- function(pair_id, split_a, split_b, tr_a, tr_b) {
  genes_a <- driver_genes_09m(tr_a)
  genes_b <- driver_genes_09m(tr_b)
  overlap <- intersect(genes_a, genes_b)
  union_genes <- union(genes_a, genes_b)
  status <- if (length(genes_a) == 0 || length(genes_b) == 0) "skipped_no_drivers" else "ok"
  reason <- if (identical(status, "ok")) "" else "missing option A driver table for at least one split"
  data.frame(
    pair_id = pair_id,
    split_a = split_a,
    split_b = split_b,
    option = "option_a",
    driver_n_a = length(genes_a),
    driver_n_b = length(genes_b),
    overlap_n = length(overlap),
    jaccard = if (length(union_genes) > 0) length(overlap) / length(union_genes) else NA_real_,
    top_overlap_genes = paste(utils::head(overlap, 30), collapse = ","),
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

compare_option_a_b_09m <- function(pair_id, split_values, tradeseq_index) {
  option_a_rows <- tradeseq_index[
    tradeseq_index$pair_id == pair_id &
      tradeseq_index$split_value %in% split_values &
      tradeseq_index$option == "option_a",
    ,
    drop = FALSE
  ]
  option_a_genes <- unique(unlist(lapply(seq_len(nrow(option_a_rows)), function(i) driver_genes_09m(option_a_rows[i, , drop = FALSE])), use.names = FALSE))
  option_b <- row_for_unit_09m(tradeseq_index, pair_id, "pooled", option = "option_b")
  option_b_genes <- driver_genes_09m(option_b)
  overlap <- intersect(option_a_genes, option_b_genes)
  union_genes <- union(option_a_genes, option_b_genes)
  jaccard <- if (length(union_genes) > 0) length(overlap) / length(union_genes) else NA_real_
  status <- if (length(option_a_genes) == 0 || length(option_b_genes) == 0) {
    "skipped_no_drivers"
  } else if (is.finite(jaccard) && jaccard < 0.10) {
    "warning_low_overlap"
  } else {
    "ok"
  }
  reason <- if (identical(status, "skipped_no_drivers")) {
    "option A union or option B condition driver table is empty"
  } else if (identical(status, "warning_low_overlap")) {
    "option A split drivers and option B condition drivers have low overlap"
  } else {
    ""
  }
  data.frame(
    pair_id = pair_id,
    split_values = paste(split_values, collapse = ","),
    option_a_driver_n = length(option_a_genes),
    option_b_driver_n = length(option_b_genes),
    overlap_n = length(overlap),
    jaccard = jaccard,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

read_entropy_09m <- function(pal_row, sling_row) {
  if (is.null(pal_row) || pal_row$status[[1]] != "ok" || !file_available_09m(pal_row$extra_path[[1]])) {
    return(data.frame(cell_id = character(0), entropy = numeric(0), label = character(0), stringsAsFactors = FALSE))
  }
  entropy <- read.csv(pal_row$extra_path[[1]], stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% colnames(entropy) || !"entropy" %in% colnames(entropy)) {
    return(data.frame(cell_id = character(0), entropy = numeric(0), label = character(0), stringsAsFactors = FALSE))
  }
  entropy <- entropy[, c("cell_id", "entropy"), drop = FALSE]
  entropy$cell_id <- as.character(entropy$cell_id)
  entropy$entropy <- suppressWarnings(as.numeric(entropy$entropy))
  pt <- read_slingshot_pt_09m(sling_row)
  if (nrow(pt) == 0) {
    entropy$label <- "__all__"
  } else {
    entropy <- merge(entropy, pt[, c("cell_id", "label"), drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
    entropy$label[is.na(entropy$label) | !nzchar(entropy$label)] <- "__all__"
  }
  entropy[is.finite(entropy$entropy), c("cell_id", "entropy", "label"), drop = FALSE]
}

compare_entropy_pair_09m <- function(pair_id, split_a, split_b, pal_a, pal_b, sling_a, sling_b) {
  ent_a <- read_entropy_09m(pal_a, sling_a)
  ent_b <- read_entropy_09m(pal_b, sling_b)
  label_var <- if (!is.null(sling_a)) normalize_scalar_value(sling_a$label_var[[1]], "label") else "label"
  if (nrow(ent_a) == 0 || nrow(ent_b) == 0) {
    return(data.frame(
      pair_id = pair_id, split_a = split_a, split_b = split_b, label_var = label_var,
      label_value = "__all__", n_cells_a = nrow(ent_a), n_cells_b = nrow(ent_b),
      median_entropy_a = NA_real_, median_entropy_b = NA_real_,
      wilcox_p_value = NA_real_, ks_statistic = NA_real_, ks_p_value = NA_real_,
      status = "skipped_no_entropy", reason = "missing successful Palantir entropy output",
      stringsAsFactors = FALSE
    ))
  }
  labels <- intersect(unique(ent_a$label), unique(ent_b$label))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0) {
    labels <- "__all__"
    ent_a$label <- "__all__"
    ent_b$label <- "__all__"
  }
  rows <- lapply(labels, function(label) {
    va <- ent_a$entropy[ent_a$label == label]
    vb <- ent_b$entropy[ent_b$label == label]
    cmp <- compare_numeric_09m(va, vb)
    data.frame(
      pair_id = pair_id,
      split_a = split_a,
      split_b = split_b,
      label_var = label_var,
      label_value = label,
      n_cells_a = sum(is.finite(va)),
      n_cells_b = sum(is.finite(vb)),
      median_entropy_a = safe_median(va),
      median_entropy_b = safe_median(vb),
      wilcox_p_value = cmp$wilcox_p,
      ks_statistic = cmp$ks_stat,
      ks_p_value = cmp$ks_p,
      status = cmp$status,
      reason = cmp$reason,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

pair_table_path_09m <- function(kind, pair_id) {
  file.path(cfg$trajectory_split_compare_table_dir, sprintf("%s_%s.tsv", kind, safe_id_09(pair_id)))
}

figure_path_09m <- function(pair_id) {
  out_dir <- trajectory_method_figure_dir_09(cfg, "split_compare", pair_id, "")
  ensure_dir(out_dir)
  file.path(out_dir, sprintf("Figure_Trajectory_SplitCompare_%s.png", safe_id_09(pair_id)))
}

plot_pair_summary_09m <- function(pair_id, pseudotime_df, lineage_df, driver_df, option_df, entropy_df, path) {
  rows <- list()
  if (nrow(pseudotime_df) > 0 && "ks_statistic" %in% colnames(pseudotime_df)) {
    rows[[length(rows) + 1L]] <- data.frame(metric = "pseudotime_ks", value = mean(pseudotime_df$ks_statistic, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  if (nrow(lineage_df) > 0 && "label_jaccard" %in% colnames(lineage_df)) {
    rows[[length(rows) + 1L]] <- data.frame(metric = "lineage_label_jaccard", value = mean(lineage_df$label_jaccard, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  if (nrow(driver_df) > 0 && "jaccard" %in% colnames(driver_df)) {
    rows[[length(rows) + 1L]] <- data.frame(metric = "driver_jaccard", value = mean(driver_df$jaccard, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  if (nrow(option_df) > 0 && "jaccard" %in% colnames(option_df)) {
    rows[[length(rows) + 1L]] <- data.frame(metric = "option_a_b_jaccard", value = mean(option_df$jaccard, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  if (nrow(entropy_df) > 0 && "ks_statistic" %in% colnames(entropy_df)) {
    rows[[length(rows) + 1L]] <- data.frame(metric = "entropy_ks", value = mean(entropy_df$ks_statistic, na.rm = TRUE), stringsAsFactors = FALSE)
  }
  plot_df <- if (length(rows) > 0) dplyr::bind_rows(rows) else data.frame(metric = character(0), value = numeric(0), stringsAsFactors = FALSE)
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  if (nrow(plot_df) == 0) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No split comparison metrics") +
      ggplot2::theme_void() +
      ggplot2::labs(title = sprintf("Split comparison %s", pair_id))
  } else {
    plot_df$metric <- factor(plot_df$metric, levels = plot_df$metric)
    plot_obj <- ggplot2::ggplot(plot_df, ggplot2::aes(x = metric, y = value, fill = metric)) +
      ggplot2::geom_col(width = 0.72) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = sprintf("Split comparison %s", pair_id), x = NULL, y = "summary value")
  }
  save_plot_dual(plot_obj, path, width = 6.5, height = 4.2)
  invisible(path)
}

units <- trajectory_execution_units_09(cfg)
split_pair_ids <- split_pair_ids_09m(cfg)
slingshot_index <- read_slingshot_index_09m(cfg)
tradeseq_index <- read_tradeseq_index_09m(cfg)
palantir_index <- read_palantir_index_09m(cfg)

pseudotime_rows <- list()
lineage_rows <- list()
driver_rows <- list()
option_rows <- list()
entropy_rows <- list()
index_rows <- list()
dynamic_outputs <- list()

for (pair_id in split_pair_ids) {
  pair_units <- units[units$pair_id == pair_id, , drop = FALSE]
  split_values <- unique(vapply(pair_units$split_value, display_scalar_value, character(1), default = "pooled"))
  split_values <- setdiff(split_values, "pooled")
  pair_pseudo <- trajectory_empty_df_09(pseudotime_cols_09m)
  pair_lineage <- trajectory_empty_df_09(lineage_cols_09m)
  pair_driver <- trajectory_empty_df_09(driver_cols_09m)
  pair_option <- trajectory_empty_df_09(option_cols_09m)
  pair_entropy <- trajectory_empty_df_09(entropy_cols_09m)
  fig_path <- figure_path_09m(pair_id)

  if (length(split_values) >= 2L) {
    combos <- utils::combn(split_values, 2, simplify = FALSE)
    pair_pseudo_rows <- list()
    pair_lineage_rows <- list()
    pair_driver_rows <- list()
    pair_entropy_rows <- list()
    for (combo in combos) {
      split_a <- combo[[1]]
      split_b <- combo[[2]]
      sling_a <- row_for_unit_09m(slingshot_index, pair_id, split_a)
      sling_b <- row_for_unit_09m(slingshot_index, pair_id, split_b)
      tr_a <- row_for_unit_09m(tradeseq_index, pair_id, split_a, option = "option_a")
      tr_b <- row_for_unit_09m(tradeseq_index, pair_id, split_b, option = "option_a")
      pal_a <- row_for_unit_09m(palantir_index, pair_id, split_a)
      pal_b <- row_for_unit_09m(palantir_index, pair_id, split_b)

      pair_pseudo_rows[[length(pair_pseudo_rows) + 1L]] <- compare_pseudotime_pair_09m(pair_id, split_a, split_b, sling_a, sling_b)
      pair_lineage_rows[[length(pair_lineage_rows) + 1L]] <- compare_lineage_pair_09m(pair_id, split_a, split_b, sling_a, sling_b)
      pair_driver_rows[[length(pair_driver_rows) + 1L]] <- compare_driver_pair_09m(pair_id, split_a, split_b, tr_a, tr_b)
      pair_entropy_rows[[length(pair_entropy_rows) + 1L]] <- compare_entropy_pair_09m(pair_id, split_a, split_b, pal_a, pal_b, sling_a, sling_b)
    }
    pair_pseudo <- dplyr::bind_rows(pair_pseudo_rows)
    pair_lineage <- dplyr::bind_rows(pair_lineage_rows)
    pair_driver <- dplyr::bind_rows(pair_driver_rows)
    pair_option <- compare_option_a_b_09m(pair_id, split_values, tradeseq_index)
    pair_entropy <- dplyr::bind_rows(pair_entropy_rows)
  }

  pseudo_path <- pair_table_path_09m("pseudotime_distribution_compare", pair_id)
  lineage_path <- pair_table_path_09m("lineage_backbone_match", pair_id)
  driver_path <- pair_table_path_09m("driver_overlap", pair_id)
  option_path <- pair_table_path_09m("option_a_vs_b_consistency", pair_id)
  entropy_path <- pair_table_path_09m("entropy_compare", pair_id)
  write_tsv_local(pair_pseudo, pseudo_path)
  write_tsv_local(pair_lineage, lineage_path)
  write_tsv_local(pair_driver, driver_path)
  write_tsv_local(pair_option, option_path)
  write_tsv_local(pair_entropy, entropy_path)
  plot_pair_summary_09m(pair_id, pair_pseudo, pair_lineage, pair_driver, pair_option, pair_entropy, fig_path)

  pseudotime_rows[[length(pseudotime_rows) + 1L]] <- pair_pseudo
  lineage_rows[[length(lineage_rows) + 1L]] <- pair_lineage
  driver_rows[[length(driver_rows) + 1L]] <- pair_driver
  option_rows[[length(option_rows) + 1L]] <- pair_option
  entropy_rows[[length(entropy_rows) + 1L]] <- pair_entropy

  all_status <- c(pair_pseudo$status, pair_lineage$status, pair_driver$status, pair_option$status, pair_entropy$status)
  all_reason <- unique(c(pair_pseudo$reason, pair_lineage$reason, pair_driver$reason, pair_option$reason, pair_entropy$reason))
  all_reason <- all_reason[nzchar(all_reason)]
  status <- if (length(split_values) < 2L) {
    "skipped_not_split"
  } else if (any(grepl("^warning|failed", all_status))) {
    "warning"
  } else if (all(all_status %in% c("ok", "skipped_no_entropy", "skipped_no_drivers", "skipped_low_cells", "skipped_no_slingshot"))) {
    "ok_partial"
  } else {
    "ok"
  }
  if (length(split_values) >= 2L && any(all_status == "ok")) {
    status <- if (any(grepl("^warning", all_status))) "warning" else "ok"
  }
  index_rows[[length(index_rows) + 1L]] <- data.frame(
    pair_id = pair_id,
    split_values = paste(split_values, collapse = ","),
    pseudotime_distribution_path = pseudo_path,
    lineage_backbone_path = lineage_path,
    driver_overlap_path = driver_path,
    option_a_vs_b_consistency_path = option_path,
    entropy_compare_path = entropy_path,
    figure_path = fig_path,
    report_path = cfg$trajectory_split_compare_report_md,
    status = status,
    reason = paste(all_reason, collapse = "; "),
    stringsAsFactors = FALSE
  )

  dynamic_outputs[[sprintf("trajectory_split_compare_%s", safe_id_09(pair_id))]] <- build_output_entry(
    fig_path,
    "png",
    module_name,
    "split-mode trajectory comparison figure",
    base_dir = cfg$project_root
  )
}

pseudotime_df <- if (length(pseudotime_rows) > 0) dplyr::bind_rows(pseudotime_rows) else trajectory_empty_df_09(pseudotime_cols_09m)
lineage_df <- if (length(lineage_rows) > 0) dplyr::bind_rows(lineage_rows) else trajectory_empty_df_09(lineage_cols_09m)
driver_df <- if (length(driver_rows) > 0) dplyr::bind_rows(driver_rows) else trajectory_empty_df_09(driver_cols_09m)
option_df <- if (length(option_rows) > 0) dplyr::bind_rows(option_rows) else trajectory_empty_df_09(option_cols_09m)
entropy_df <- if (length(entropy_rows) > 0) dplyr::bind_rows(entropy_rows) else trajectory_empty_df_09(entropy_cols_09m)
split_index <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(split_index_cols_09m)

write_tsv_local(pseudotime_df, cfg$trajectory_pseudotime_distribution_compare_tsv)
write_tsv_local(lineage_df, cfg$trajectory_lineage_backbone_match_tsv)
write_tsv_local(driver_df, cfg$trajectory_driver_overlap_tsv)
write_tsv_local(option_df, cfg$trajectory_option_a_b_consistency_tsv)
write_tsv_local(entropy_df, cfg$trajectory_entropy_compare_tsv)
write_tsv_local(split_index, cfg$trajectory_split_compare_index_tsv)

status_counts <- if (nrow(split_index) > 0) {
  as.data.frame(table(split_index$status), stringsAsFactors = FALSE)
} else {
  data.frame(status = character(0), n = integer(0), stringsAsFactors = FALSE)
}
if (nrow(status_counts) > 0) {
  colnames(status_counts) <- c("status", "n")
}

report_lines <- build_report_lines_v04(
  title = "09m Trajectory Split Compare",
  header_bullets = c(
    sprintf("Split-mode trajectory pairs compared: %s.", nrow(split_index)),
    "Comparisons include pseudotime distributions, Slingshot label backbones, tradeSeq driver overlap, option A/B consistency, and Palantir entropy."
  ),
  key_files = list(
    split_compare_index = cfg$trajectory_split_compare_index_tsv,
    pseudotime_distribution_compare = cfg$trajectory_pseudotime_distribution_compare_tsv,
    lineage_backbone_match = cfg$trajectory_lineage_backbone_match_tsv,
    driver_overlap = cfg$trajectory_driver_overlap_tsv,
    option_a_vs_b_consistency = cfg$trajectory_option_a_b_consistency_tsv,
    entropy_compare = cfg$trajectory_entropy_compare_tsv
  ),
  review_focus = c(
    "Treat split differences as descriptive unless supported by tradeSeq option B and upstream biological replication.",
    "Low option A/B overlap is a consistency warning, not an automatic failure."
  ),
  extra_sections = list(
    "Status Counts" = render_markdown_table_local(status_counts),
    "Split Compare Index" = render_markdown_table_local(split_index),
    "Pseudotime Distribution Preview" = render_markdown_table_local(utils::head(pseudotime_df, 30)),
    "Lineage Backbone Preview" = render_markdown_table_local(utils::head(lineage_df, 30)),
    "Driver Overlap Preview" = render_markdown_table_local(utils::head(driver_df, 30)),
    "Option A vs B Consistency" = render_markdown_table_local(option_df),
    "Entropy Compare Preview" = render_markdown_table_local(utils::head(entropy_df, 30))
  )
)
ensure_dir(dirname(cfg$trajectory_split_compare_report_md))
write_markdown_local(report_lines, cfg$trajectory_split_compare_report_md)

outputs <- c(
  list(
    report_md = build_output_entry(cfg$trajectory_split_compare_report_md, "md", module_name, "split-mode trajectory comparison report", base_dir = cfg$project_root),
    split_compare_index_tsv = build_output_entry(cfg$trajectory_split_compare_index_tsv, "tsv", module_name, "split comparison index by trajectory pair", base_dir = cfg$project_root, schema = infer_schema_from_df(split_index)),
    pseudotime_distribution_compare_tsv = build_output_entry(cfg$trajectory_pseudotime_distribution_compare_tsv, "tsv", module_name, "pseudotime distribution tests across split values", base_dir = cfg$project_root, schema = infer_schema_from_df(pseudotime_df)),
    lineage_backbone_match_tsv = build_output_entry(cfg$trajectory_lineage_backbone_match_tsv, "tsv", module_name, "Slingshot lineage backbone agreement across split values", base_dir = cfg$project_root, schema = infer_schema_from_df(lineage_df)),
    driver_overlap_tsv = build_output_entry(cfg$trajectory_driver_overlap_tsv, "tsv", module_name, "tradeSeq option A driver overlap across split values", base_dir = cfg$project_root, schema = infer_schema_from_df(driver_df)),
    option_a_vs_b_consistency_tsv = build_output_entry(cfg$trajectory_option_a_b_consistency_tsv, "tsv", module_name, "tradeSeq option A/B driver consistency", base_dir = cfg$project_root, schema = infer_schema_from_df(option_df)),
    entropy_compare_tsv = build_output_entry(cfg$trajectory_entropy_compare_tsv, "tsv", module_name, "Palantir entropy distribution tests across split values", base_dir = cfg$project_root, schema = infer_schema_from_df(entropy_df))
  ),
  dynamic_outputs
)

trajectory_method_manifest_09(
  cfg,
  cfg$module_09m_manifest_path,
  module_name,
  outputs,
  inputs = list(
    module_09a = cfg$module_09a_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09g = cfg$module_09g_manifest_path,
    module_09h = cfg$module_09h_manifest_path,
    trajectory_pairs = cfg$trajectory_pairs_sheet
  ),
  depends_on = list(
    module_09a = cfg$module_09a_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09h = cfg$module_09h_manifest_path
  )
)

message("09m completed. split compare index: ", cfg$trajectory_split_compare_index_tsv)
