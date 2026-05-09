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

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09e2_trajectory_monocle2"
prepare_dirs_09(cfg)

monocle2_index_cols <- c(trajectory_method_index_cols_09, "branches_path", "root_state")

state_of_root_09e2 <- function(cds, label_var, root_label) {
  pheno <- Biobase::pData(cds)
  if (!"State" %in% colnames(pheno) || !label_var %in% colnames(pheno)) {
    return(NA)
  }
  root_rows <- as.character(pheno[[label_var]]) == root_label
  if (!any(root_rows)) {
    return(NA)
  }
  tab <- sort(table(pheno$State[root_rows]), decreasing = TRUE)
  suppressWarnings(as.numeric(names(tab)[[1]]))
}

run_monocle2_09e2 <- function(seu, label_var, root_label, out_csv, tree_rds, fig_png, branches_tsv) {
  if (!requireNamespace("monocle", quietly = TRUE)) {
    stop("missing R package: monocle", call. = FALSE)
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("missing R package: Biobase", call. = FALSE)
  }
  counts <- trajectory_counts_matrix_09(seu)
  pd <- Biobase::AnnotatedDataFrame(data = seu@meta.data)
  fd <- Biobase::AnnotatedDataFrame(data = data.frame(gene_short_name = rownames(counts), row.names = rownames(counts)))
  cds <- monocle::newCellDataSet(
    as(counts, "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    expressionFamily = monocle::negbinomial.size()
  )
  cds <- monocle::estimateSizeFactors(cds)
  cds <- monocle::estimateDispersions(cds)
  cds <- monocle::detectGenes(cds, min_expr = 0.1)
  expressed <- rownames(Biobase::fData(cds))[Biobase::fData(cds)$num_cells_expressed >= max(10L, ceiling(ncol(seu) * 0.02))]
  if (length(expressed) < 20L) {
    expressed <- rownames(counts)
  }
  n_cores <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
  n_cores <- if (length(n_cores) == 0 || is.na(n_cores)) 1L else max(1L, min(4L, n_cores))
  day_genes <- tryCatch(
    monocle::differentialGeneTest(
      cds[expressed, ],
      fullModelFormulaStr = paste0("~", label_var),
      cores = n_cores
    ),
    error = function(e) {
      stop(sprintf("Monocle 2 differentialGeneTest failed; refusing q=1 ordering fallback: %s", conditionMessage(e)), call. = FALSE)
    }
  )
  if (!"qval" %in% colnames(day_genes)) {
    stop("Monocle 2 differentialGeneTest output lacks qval; refusing ordering fallback", call. = FALSE)
  }
  day_genes$qval <- suppressWarnings(as.numeric(day_genes$qval))
  if (!any(is.finite(day_genes$qval))) {
    stop("Monocle 2 differentialGeneTest returned no finite qval; refusing ordering fallback", call. = FALSE)
  }
  ordering_genes <- rownames(day_genes)[order(day_genes$qval)]
  ordering_genes <- utils::head(ordering_genes[nzchar(ordering_genes)], min(2000L, length(ordering_genes)))
  if (length(ordering_genes) < 20L) {
    stop("Monocle 2 differentialGeneTest returned fewer than 20 ordering genes", call. = FALSE)
  }
  cds <- monocle::setOrderingFilter(cds, ordering_genes = ordering_genes)
  cds <- monocle::reduceDimension(cds, method = "DDRTree", max_components = 2, norm_method = "none", pseudo_expr = 0)
  root_state <- state_of_root_09e2(cds, label_var, root_label)
  cds <- if (is.finite(root_state)) monocle::orderCells(cds, root_state = root_state) else monocle::orderCells(cds)
  pheno <- Biobase::pData(cds)
  pst <- trajectory_scale01_09(pheno$Pseudotime)
  names(pst) <- rownames(pheno)
  extra <- data.frame(State = pheno$State, stringsAsFactors = FALSE)
  trajectory_write_pseudotime_09(out_csv, names(pst), pst, label = pheno[[label_var]], extra = extra)
  branches <- as.data.frame(table(State = pheno$State), stringsAsFactors = FALSE)
  write_tsv_local(branches, branches_tsv)
  saveRDS(cds, tree_rds)
  png(fig_png, width = 1800, height = 1400, res = 220)
  print(monocle::plot_cell_trajectory(cds, color_by = label_var))
  dev.off()
  root_state
}

units_all <- trajectory_execution_units_09(cfg)
enabled_mask <- if (nrow(units_all) > 0) {
  vapply(seq_len(nrow(units_all)), function(i) {
    trajectory_method_enabled_09(units_all[i, , drop = FALSE], "monocle2")
  }, logical(1))
} else {
  logical(0)
}
units <- units_all[enabled_mask, , drop = FALSE]

index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  out_csv <- trajectory_method_output_path_09(cfg, "monocle2", pair_id, split_value, "monocle2_pseudotime")
  tree_rds <- trajectory_method_output_path_09(cfg, "monocle2", pair_id, split_value, "monocle2_tree", "rds")
  branches_tsv <- trajectory_method_output_path_09(cfg, "monocle2", pair_id, split_value, "monocle2_branches", "tsv")
  fig_png <- trajectory_method_figure_path_09(cfg, "monocle2", pair_id, split_value, "Figure_M2")
  ensure_dir(dirname(out_csv))
  ensure_dir(dirname(fig_png))
  started <- proc.time()[["elapsed"]]
  row <- tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])
    root_state <- run_monocle2_09e2(seu, unit$coarse_label_var[[1]], root, out_csv, tree_rds, fig_png, branches_tsv)
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "monocle2",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = tree_rds,
      figure_path = fig_png,
      n_cells = ncol(seu),
      status = "ok",
      reason = "",
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      branches_path = branches_tsv,
      root_state = root_state,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), out_csv, row.names = FALSE)
    write_tsv_local(data.frame(State = character(0), Freq = integer(0)), branches_tsv)
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "monocle2",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = tree_rds,
      figure_path = fig_png,
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = if (grepl("missing R package", conditionMessage(e))) "skipped_package_missing" else "failed",
      reason = conditionMessage(e),
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      branches_path = branches_tsv,
      root_state = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  index_rows[[length(index_rows) + 1L]] <- row
  dynamic_outputs[[sprintf("monocle2_%s_pseudotime", trajectory_unit_file_id_09(pair_id, split_value))]] <- build_output_entry(out_csv, "csv", module_name, "Monocle 2 pseudotime", base_dir = cfg$project_root)
}

if (nrow(units_all) > 0 && any(!enabled_mask)) {
  disabled <- units_all[!enabled_mask, , drop = FALSE]
  disabled_rows <- lapply(seq_len(nrow(disabled)), function(i) {
    unit <- disabled[i, , drop = FALSE]
    pair_id <- unit$pair_id[[1]]
    split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "monocle2",
      methods_enabled = "no",
      input_rds = unit$input_rds[[1]],
      output_path = "",
      extra_path = "",
      figure_path = "",
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = "skipped_disabled",
      reason = trajectory_method_disabled_reason_09(unit, "monocle2"),
      runtime_s = 0,
      branches_path = "",
      root_state = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  index_rows <- c(index_rows, disabled_rows)
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(monocle2_index_cols)
write_tsv_local(index_df, cfg$trajectory_monocle2_index_tsv)

outputs <- c(
  list(monocle2_index_tsv = build_output_entry(cfg$trajectory_monocle2_index_tsv, "tsv", module_name, "Monocle 2 opt-in status by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09e2_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path, trajectory_pairs = cfg$trajectory_pairs_sheet),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path)
)

message("09e2 completed. Monocle 2 index: ", cfg$trajectory_monocle2_index_tsv)
