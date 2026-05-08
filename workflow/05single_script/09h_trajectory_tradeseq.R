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
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_09.R"))
source_utf8(file.path(.script_dir, "helpers", "trajectory_utils.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09h_trajectory_tradeseq"
prepare_dirs_09(cfg)

tradeseq_index_cols <- c(trajectory_method_index_cols_09, "option", "association_path", "pattern_path", "driver_path", "conditiontest_path")
tradeseq_driver_cols <- c("pair_id", "split_value", "option", "driver_path", "n_drivers", "status", "reason")

empty_test_df_09h <- function() {
  data.frame(gene_id = character(0), waldStat = numeric(0), df = numeric(0), pvalue = numeric(0), padj = numeric(0), stringsAsFactors = FALSE)
}

write_empty_test_09h <- function(path) {
  ensure_dir(dirname(path))
  write_tsv_local(empty_test_df_09h(), path)
}

normalize_test_df_09h <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return(empty_test_df_09h())
  }
  df <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  df$gene_id <- rownames(df)
  rownames(df) <- NULL
  p_cols <- grep("^p(value)?$|pval|p\\.value", colnames(df), ignore.case = TRUE, value = TRUE)
  if (!"pvalue" %in% colnames(df)) {
    df$pvalue <- if (length(p_cols) > 0) suppressWarnings(as.numeric(df[[p_cols[[1]]]])) else NA_real_
  }
  if (!"padj" %in% colnames(df)) {
    df$padj <- stats::p.adjust(df$pvalue, method = "BH")
  }
  df <- df[, unique(c("gene_id", setdiff(colnames(df), "gene_id"))), drop = FALSE]
  df[order(df$padj, df$pvalue, na.last = TRUE), , drop = FALSE]
}

read_slingshot_index_09h <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09d_manifest_path,
    "slingshot_index_tsv",
    cfg$trajectory_slingshot_index_tsv,
    cols = c(trajectory_method_index_cols_09, "fine_output_path", "lineage_count", "root_label", "label_var")
  )
}

slingshot_row_for_unit_09h <- function(slingshot_index, pair_id, split_value) {
  hit <- slingshot_index[
    slingshot_index$pair_id == pair_id &
      slingshot_index$split_value == display_scalar_value(split_value, "pooled"),
    ,
    drop = FALSE
  ]
  if (nrow(hit) == 0) {
    return(NULL)
  }
  hit[1, , drop = FALSE]
}

tradeseq_inputs_09h <- function(seu, pseudotime_csv) {
  pt <- trajectory_read_pseudotime_csv_09(pseudotime_csv, "slingshot")
  counts <- trajectory_counts_matrix_09(seu)
  common <- intersect(colnames(counts), pt$cell_id)
  if (length(common) < 10L) {
    stop("fewer than 10 cells overlap Seurat counts and Slingshot pseudotime", call. = FALSE)
  }
  counts <- counts[, common, drop = FALSE]
  pt <- pt[match(common, pt$cell_id), , drop = FALSE]
  lineage_cols <- grep("^lineage_", colnames(pt), value = TRUE)
  if (length(lineage_cols) == 0) {
    pseudotime <- matrix(pt$pseudotime, ncol = 1)
  } else {
    pseudotime <- as.matrix(pt[, lineage_cols, drop = FALSE])
  }
  pseudotime <- apply(pseudotime, 2, trajectory_scale01_09)
  if (is.null(dim(pseudotime))) {
    pseudotime <- matrix(pseudotime, ncol = 1)
  }
  cell_weights <- ifelse(is.finite(pseudotime), 1, 0)
  pseudotime[!is.finite(pseudotime)] <- 0
  rownames(pseudotime) <- common
  rownames(cell_weights) <- common
  list(counts = counts, pseudotime = pseudotime, cell_weights = cell_weights, cells = common)
}

fit_tradeseq_09h <- function(counts, pseudotime, cell_weights, conditions = NULL) {
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("missing R package: tradeSeq", call. = FALSE)
  }
  args <- list(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    nknots = cfg$tradeseq_knots,
    verbose = FALSE
  )
  if (!is.null(conditions)) {
    args$conditions <- conditions
  }
  do.call(tradeSeq::fitGAM, args)
}

plot_driver_scores_09h <- function(df, path, title) {
  ensure_dir(dirname(path))
  if (nrow(df) == 0 || !"padj" %in% colnames(df)) {
    plot_obj <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No tradeSeq driver data") +
      ggplot2::theme_void()
  } else {
    show <- utils::head(df[order(df$padj, na.last = TRUE), , drop = FALSE], 30)
    show$gene_id <- factor(show$gene_id, levels = rev(show$gene_id))
    show$score <- -log10(pmax(suppressWarnings(as.numeric(show$padj)), .Machine$double.xmin))
    plot_obj <- ggplot2::ggplot(show, ggplot2::aes(x = score, y = gene_id)) +
      ggplot2::geom_col(fill = "#4C78A8") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(title = title, x = "-log10(BH adjusted p)", y = NULL)
  }
  save_plot_local(plot_obj, path, width = 6.5, height = 6)
}

driver_table_09h <- function(assoc_df, pattern_df, pair_id, split_value, option) {
  rows <- list()
  if (nrow(assoc_df) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(gene_id = assoc_df$gene_id, assoc_padj = assoc_df$padj, stringsAsFactors = FALSE)
  }
  if (nrow(pattern_df) > 0) {
    rows[[length(rows) + 1L]] <- data.frame(gene_id = pattern_df$gene_id, pattern_padj = pattern_df$padj, stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) {
    return(data.frame(pair_id = character(0), split_value = character(0), option = character(0), gene_id = character(0), driver_score = numeric(0), stringsAsFactors = FALSE))
  }
  merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE, sort = FALSE), rows)
  p_cols <- setdiff(colnames(merged), "gene_id")
  merged$driver_padj <- apply(merged[, p_cols, drop = FALSE], 1, function(x) min(suppressWarnings(as.numeric(x)), na.rm = TRUE))
  merged$driver_padj[!is.finite(merged$driver_padj)] <- NA_real_
  merged$driver_score <- -log10(pmax(merged$driver_padj, .Machine$double.xmin))
  merged <- merged[order(merged$driver_padj, na.last = TRUE), , drop = FALSE]
  merged <- utils::head(merged, 200)
  cbind(
    data.frame(pair_id = pair_id, split_value = split_value, option = option, stringsAsFactors = FALSE),
    merged
  )
}

run_option_a_09h <- function(unit, sling_row) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  assoc_tsv <- trajectory_method_output_path_09(cfg, "tradeseq", pair_id, split_value, "tradeseq_assoc", "tsv")
  pattern_tsv <- trajectory_method_output_path_09(cfg, "tradeseq", pair_id, split_value, "tradeseq_pattern", "tsv")
  driver_tsv <- trajectory_method_output_path_09(cfg, "tradeseq", pair_id, split_value, "tradeseq_drivers_top200", "tsv")
  fig_png <- trajectory_method_figure_path_09(cfg, "tradeseq", pair_id, split_value, "Figure_5G_DriverHeatmap")
  started <- proc.time()[["elapsed"]]

  row <- tryCatch({
    if (is.null(sling_row) || sling_row$status[[1]] != "ok" || !file.exists(sling_row$output_path[[1]])) {
      stop("missing successful Slingshot pseudotime", call. = FALSE)
    }
    seu <- readRDS(unit$input_rds[[1]])
    inp <- tradeseq_inputs_09h(seu, sling_row$output_path[[1]])
    fit <- fit_tradeseq_09h(inp$counts, inp$pseudotime, inp$cell_weights)
    assoc <- normalize_test_df_09h(tradeSeq::associationTest(fit))
    pattern <- tryCatch(normalize_test_df_09h(tradeSeq::patternTest(fit)), error = function(e) empty_test_df_09h())
    drivers <- driver_table_09h(assoc, pattern, pair_id, split_value, "option_a")
    write_tsv_local(assoc, assoc_tsv)
    write_tsv_local(pattern, pattern_tsv)
    write_tsv_local(drivers, driver_tsv)
    plot_driver_scores_09h(drivers, fig_png, sprintf("tradeSeq drivers %s %s", pair_id, split_value))
    data.frame(
      pair_id = pair_id, split_value = split_value, method = "tradeseq", methods_enabled = "yes",
      input_rds = unit$input_rds[[1]], output_path = assoc_tsv, extra_path = pattern_tsv,
      figure_path = fig_png, n_cells = ncol(seu), status = "ok", reason = "",
      runtime_s = round(proc.time()[["elapsed"]] - started, 3), option = "option_a",
      association_path = assoc_tsv, pattern_path = pattern_tsv, driver_path = driver_tsv,
      conditiontest_path = "", stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write_empty_test_09h(assoc_tsv)
    write_empty_test_09h(pattern_tsv)
    write_tsv_local(data.frame(pair_id = character(0), split_value = character(0), option = character(0), gene_id = character(0), driver_score = numeric(0)), driver_tsv)
    status <- if (grepl("missing R package", conditionMessage(e))) "skipped_package_missing" else "skipped_no_slingshot"
    if (!grepl("missing R package|Slingshot", conditionMessage(e))) {
      status <- "failed"
    }
    data.frame(
      pair_id = pair_id, split_value = split_value, method = "tradeseq", methods_enabled = "yes",
      input_rds = unit$input_rds[[1]], output_path = assoc_tsv, extra_path = pattern_tsv,
      figure_path = fig_png, n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = status, reason = conditionMessage(e), runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      option = "option_a", association_path = assoc_tsv, pattern_path = pattern_tsv,
      driver_path = driver_tsv, conditiontest_path = "", stringsAsFactors = FALSE
    )
  })
  row
}

run_option_b_09h <- function(pair_id, units_pair, sling_index) {
  split_values <- unique(vapply(units_pair$split_value, display_scalar_value, character(1), default = "pooled"))
  condition_tsv <- trajectory_method_output_path_09(cfg, "tradeseq", pair_id, "", "tradeseq_conditiontest", "tsv")
  driver_tsv <- trajectory_method_output_path_09(cfg, "tradeseq", pair_id, "", "tradeseq_drivers_condition_top200", "tsv")
  fig_png <- trajectory_method_figure_path_09(cfg, "tradeseq", pair_id, "", "Figure_5G_ConditionTest")
  started <- proc.time()[["elapsed"]]
  tryCatch({
    if (length(split_values) < 2 || all(split_values == "pooled")) {
      stop("pair is not split-mode; conditionTest not applicable", call. = FALSE)
    }
    counts_list <- list()
    pt_list <- list()
    cond <- character(0)
    for (i in seq_len(nrow(units_pair))) {
      unit <- units_pair[i, , drop = FALSE]
      split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
      sling <- slingshot_row_for_unit_09h(sling_index, pair_id, split_value)
      if (is.null(sling) || sling$status[[1]] != "ok") {
        next
      }
      seu <- readRDS(unit$input_rds[[1]])
      inp <- tradeseq_inputs_09h(seu, sling$output_path[[1]])
      counts_list[[split_value]] <- inp$counts
      pt_list[[split_value]] <- inp
      cond <- c(cond, rep(split_value, length(inp$cells)))
    }
    if (length(counts_list) < 2) {
      stop("fewer than two split values have successful Slingshot pseudotime", call. = FALSE)
    }
    common_genes <- Reduce(intersect, lapply(counts_list, rownames))
    counts <- do.call(cbind, lapply(counts_list, function(x) x[common_genes, , drop = FALSE]))
    pseudotime <- do.call(rbind, lapply(pt_list, function(x) x$pseudotime[, 1, drop = FALSE]))
    cell_weights <- do.call(rbind, lapply(pt_list, function(x) x$cell_weights[, 1, drop = FALSE]))
    fit <- fit_tradeseq_09h(counts, pseudotime, cell_weights, conditions = factor(cond))
    condition_df <- normalize_test_df_09h(tradeSeq::conditionTest(fit))
    drivers <- cbind(
      data.frame(pair_id = pair_id, split_value = "pooled", option = "option_b", stringsAsFactors = FALSE),
      utils::head(condition_df, 200)
    )
    write_tsv_local(condition_df, condition_tsv)
    write_tsv_local(drivers, driver_tsv)
    plot_driver_scores_09h(condition_df, fig_png, sprintf("tradeSeq conditionTest %s", pair_id))
    data.frame(
      pair_id = pair_id, split_value = "pooled", method = "tradeseq", methods_enabled = "yes",
      input_rds = paste(units_pair$input_rds, collapse = ";"), output_path = condition_tsv,
      extra_path = driver_tsv, figure_path = fig_png, n_cells = ncol(counts), status = "ok",
      reason = "", runtime_s = round(proc.time()[["elapsed"]] - started, 3), option = "option_b",
      association_path = "", pattern_path = "", driver_path = driver_tsv,
      conditiontest_path = condition_tsv, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write_empty_test_09h(condition_tsv)
    write_tsv_local(data.frame(pair_id = character(0), split_value = character(0), option = character(0), gene_id = character(0), driver_score = numeric(0)), driver_tsv)
    status <- if (grepl("not split-mode", conditionMessage(e))) "skipped_not_split" else if (grepl("missing R package", conditionMessage(e))) "skipped_package_missing" else "skipped_no_slingshot"
    data.frame(
      pair_id = pair_id, split_value = "pooled", method = "tradeseq", methods_enabled = "yes",
      input_rds = paste(units_pair$input_rds, collapse = ";"), output_path = condition_tsv,
      extra_path = driver_tsv, figure_path = fig_png, n_cells = sum(suppressWarnings(as.integer(units_pair$cell_n_after)), na.rm = TRUE),
      status = status, reason = conditionMessage(e), runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      option = "option_b", association_path = "", pattern_path = "", driver_path = driver_tsv,
      conditiontest_path = condition_tsv, stringsAsFactors = FALSE
    )
  })
}

units <- trajectory_execution_units_09(cfg)
slingshot_index <- read_slingshot_index_09h(cfg)
index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  row <- run_option_a_09h(unit, slingshot_row_for_unit_09h(slingshot_index, unit$pair_id[[1]], split_value))
  index_rows[[length(index_rows) + 1L]] <- row
  unit_id <- trajectory_unit_file_id_09(unit$pair_id[[1]], split_value)
  dynamic_outputs[[sprintf("tradeseq_assoc_%s", unit_id)]] <- build_output_entry(row$association_path[[1]], "tsv", module_name, "tradeSeq associationTest output", base_dir = cfg$project_root)
  dynamic_outputs[[sprintf("tradeseq_pattern_%s", unit_id)]] <- build_output_entry(row$pattern_path[[1]], "tsv", module_name, "tradeSeq patternTest output", base_dir = cfg$project_root)
  dynamic_outputs[[sprintf("tradeseq_drivers_%s_top200", unit_id)]] <- build_output_entry(row$driver_path[[1]], "tsv", module_name, "top tradeSeq trajectory drivers", base_dir = cfg$project_root)
}

if (nrow(units) > 0) {
  for (pair_id in unique(units$pair_id)) {
    units_pair <- units[units$pair_id == pair_id, , drop = FALSE]
    if (length(unique(vapply(units_pair$split_value, display_scalar_value, character(1), default = "pooled"))) > 1) {
      row <- run_option_b_09h(pair_id, units_pair, slingshot_index)
      index_rows[[length(index_rows) + 1L]] <- row
      dynamic_outputs[[sprintf("tradeseq_conditiontest_%s", safe_id_09(pair_id))]] <- build_output_entry(row$conditiontest_path[[1]], "tsv", module_name, "tradeSeq conditionTest output for split-mode pair", base_dir = cfg$project_root)
      dynamic_outputs[[sprintf("tradeseq_condition_drivers_%s_top200", safe_id_09(pair_id))]] <- build_output_entry(row$driver_path[[1]], "tsv", module_name, "top condition-dependent trajectory drivers", base_dir = cfg$project_root)
    }
  }
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(tradeseq_index_cols)
write_tsv_local(index_df, cfg$trajectory_tradeseq_index_tsv)

condition_index <- index_df[index_df$option == "option_b", , drop = FALSE]
driver_index <- index_df[nzchar(index_df$driver_path), c("pair_id", "split_value", "option", "driver_path", "status", "reason"), drop = FALSE]
driver_index$n_drivers <- vapply(driver_index$driver_path, function(path) nrow(read_tsv_optional(path)), integer(1))
driver_index <- driver_index[, tradeseq_driver_cols, drop = FALSE]
write_tsv_local(condition_index, cfg$trajectory_tradeseq_condition_index_tsv)
write_tsv_local(driver_index, cfg$trajectory_tradeseq_driver_index_tsv)

outputs <- c(
  list(
    tradeseq_index_tsv = build_output_entry(cfg$trajectory_tradeseq_index_tsv, "tsv", module_name, "tradeSeq status by pair/split and option", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
    tradeseq_condition_index_tsv = build_output_entry(cfg$trajectory_tradeseq_condition_index_tsv, "tsv", module_name, "tradeSeq split-mode conditionTest status", base_dir = cfg$project_root, schema = infer_schema_from_df(condition_index)),
    tradeseq_driver_index_tsv = build_output_entry(cfg$trajectory_tradeseq_driver_index_tsv, "tsv", module_name, "tradeSeq top-driver file index", base_dir = cfg$project_root, schema = infer_schema_from_df(driver_index))
  ),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09h_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09d = cfg$module_09d_manifest_path, trajectory_pairs = cfg$trajectory_pairs_sheet),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09d = cfg$module_09d_manifest_path)
)

message("09h completed. tradeSeq index: ", cfg$trajectory_tradeseq_index_tsv)
