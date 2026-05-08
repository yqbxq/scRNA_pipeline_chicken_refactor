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
module_name <- "09e_trajectory_monocle3"
prepare_dirs_09(cfg)

monocle3_index_cols <- c(trajectory_method_index_cols_09, "root_cell_n", "root_label")

run_monocle3_09e <- function(seu, label_var, root_label, out_csv, graph_rds, fig_png, pair_id, split_value) {
  if (!requireNamespace("monocle3", quietly = TRUE)) {
    stop("missing R package: monocle3", call. = FALSE)
  }
  counts <- trajectory_counts_matrix_09(seu)
  cell_meta <- seu@meta.data
  cell_meta$cell_id <- rownames(cell_meta)
  gene_meta <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle3::new_cell_data_set(counts, cell_metadata = cell_meta, gene_metadata = gene_meta)
  cds <- monocle3::preprocess_cds(cds, num_dim = min(50L, ncol(seu) - 1L), verbose = FALSE)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", verbose = FALSE)
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE, verbose = FALSE)
  root_cells <- trajectory_root_cells_09(seu, label_var, root_label)
  if (length(root_cells) > 0) {
    cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
  } else {
    cds <- monocle3::order_cells(cds, reduction_method = "UMAP")
  }
  pst <- monocle3::pseudotime(cds)
  pst <- trajectory_scale01_09(pst)
  names(pst) <- colnames(cds)
  label <- if (label_var %in% colnames(seu@meta.data)) as.character(seu@meta.data[names(pst), label_var]) else NULL
  trajectory_write_pseudotime_09(out_csv, names(pst), pst, label = label)
  saveRDS(cds, graph_rds)
  trajectory_plot_pseudotime_09(seu, pst, fig_png, sprintf("Monocle 3 %s %s", pair_id, split_value))
  length(root_cells)
}

units <- trajectory_execution_units_09(cfg)
index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  out_csv <- trajectory_method_output_path_09(cfg, "monocle3", pair_id, split_value, "monocle3_pseudotime")
  graph_rds <- trajectory_method_output_path_09(cfg, "monocle3", pair_id, split_value, "monocle3_graph", "rds")
  fig_png <- trajectory_method_figure_path_09(cfg, "monocle3", pair_id, split_value, "Figure_Monocle3")
  ensure_dir(dirname(out_csv))
  ensure_dir(dirname(fig_png))
  started <- proc.time()[["elapsed"]]

  row <- tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])
    root_n <- run_monocle3_09e(seu, unit$coarse_label_var[[1]], root, out_csv, graph_rds, fig_png, pair_id, split_value)
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "monocle3",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = graph_rds,
      figure_path = fig_png,
      n_cells = ncol(seu),
      status = "ok",
      reason = "",
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      root_cell_n = root_n,
      root_label = root,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), out_csv, row.names = FALSE)
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "monocle3",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = graph_rds,
      figure_path = fig_png,
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = "failed",
      reason = conditionMessage(e),
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      root_cell_n = 0L,
      root_label = "",
      stringsAsFactors = FALSE
    )
  })
  index_rows[[length(index_rows) + 1L]] <- row
  dynamic_outputs[[sprintf("monocle3_%s_pseudotime", trajectory_unit_file_id_09(pair_id, split_value))]] <- build_output_entry(out_csv, "csv", module_name, "Monocle 3 pseudotime", base_dir = cfg$project_root)
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(monocle3_index_cols)
write_tsv_local(index_df, cfg$trajectory_monocle3_index_tsv)

outputs <- c(
  list(monocle3_index_tsv = build_output_entry(cfg$trajectory_monocle3_index_tsv, "tsv", module_name, "Monocle 3 method status by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09e_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path)
)

message("09e completed. Monocle 3 index: ", cfg$trajectory_monocle3_index_tsv)
