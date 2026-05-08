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

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2", "SingleCellExperiment"))

cfg <- get_single_script_config_09()
module_name <- "09d_trajectory_slingshot"
prepare_dirs_09(cfg)

slingshot_index_cols <- c(
  trajectory_method_index_cols_09,
  "fine_output_path", "lineage_count", "root_label", "label_var"
)

run_slingshot_label_09d <- function(seu, label_var, root_label) {
  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("missing R package: slingshot", call. = FALSE)
  }
  emb <- trajectory_umap_09(seu)
  if (is.null(emb)) {
    stop("input object has no UMAP reduction", call. = FALSE)
  }
  labels <- as.character(seu@meta.data[[label_var]])
  labels[is.na(labels) | !nzchar(labels)] <- "unknown"
  if (length(unique(labels)) < 2L) {
    stop(sprintf("label_var=%s has fewer than 2 labels", label_var), call. = FALSE)
  }
  sce <- Seurat::as.SingleCellExperiment(seu)
  SingleCellExperiment::reducedDim(sce, "UMAP") <- emb[colnames(sce), , drop = FALSE]
  SummarizedExperiment::colData(sce)[[label_var]] <- labels
  start.clus <- normalize_scalar_value(root_label)
  if (!nzchar(start.clus) || !start.clus %in% unique(labels)) {
    start.clus <- NULL
  }
  sce <- slingshot::slingshot(
    sce,
    clusterLabels = label_var,
    reducedDim = "UMAP",
    start.clus = start.clus
  )
  pst <- slingshot::slingPseudotime(sce)
  pst_vec <- rowMeans(pst, na.rm = TRUE)
  pst_vec[!is.finite(pst_vec)] <- NA_real_
  names(pst_vec) <- rownames(pst)
  list(sce = sce, pseudotime = pst, pseudotime_vector = pst_vec, labels = labels)
}

fine_root_label_09d <- function(seu, coarse_var, fine_var, root_label) {
  root_cells <- trajectory_root_cells_09(seu, coarse_var, root_label)
  if (length(root_cells) == 0 || !fine_var %in% colnames(seu@meta.data)) {
    return("")
  }
  tab <- sort(table(as.character(seu@meta.data[root_cells, fine_var])), decreasing = TRUE)
  if (length(tab) == 0) "" else names(tab)[[1]]
}

units <- trajectory_execution_units_09(cfg)
index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  method_dir <- trajectory_method_table_dir_09(cfg, "slingshot", pair_id, if (identical(split_value, "pooled")) "" else split_value)
  ensure_dir(method_dir)
  out_csv <- trajectory_method_output_path_09(cfg, "slingshot", pair_id, split_value, "slingshot_pseudotime")
  fine_csv <- trajectory_method_output_path_09(cfg, "slingshot", pair_id, split_value, "slingshot_pseudotime", "seurat_clusters.csv")
  curve_rds <- trajectory_method_output_path_09(cfg, "slingshot", pair_id, split_value, "slingshot_lineage_curves", "rds")
  fig_png <- trajectory_method_figure_path_09(cfg, "slingshot", pair_id, split_value, "Figure_Slingshot")
  started <- proc.time()[["elapsed"]]

  row <- tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    coarse <- normalize_scalar_value(unit$coarse_label_var[[1]])
    fine <- normalize_scalar_value(unit$fine_label_var[[1]])
    root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])

    coarse_fit <- run_slingshot_label_09d(seu, coarse, root)
    extra <- as.data.frame(coarse_fit$pseudotime, stringsAsFactors = FALSE)
    colnames(extra) <- sprintf("lineage_%s", seq_len(ncol(extra)))
    trajectory_write_pseudotime_09(out_csv, names(coarse_fit$pseudotime_vector), coarse_fit$pseudotime_vector, label = coarse_fit$labels, extra = extra)
    saveRDS(coarse_fit$sce, curve_rds)
    trajectory_plot_pseudotime_09(seu, coarse_fit$pseudotime_vector, fig_png, sprintf("Slingshot %s %s", pair_id, split_value))

    if (nzchar(fine) && fine %in% colnames(seu@meta.data) && length(unique(as.character(seu@meta.data[[fine]]))) >= 2L) {
      fine_root <- fine_root_label_09d(seu, coarse, fine, root)
      fine_fit <- run_slingshot_label_09d(seu, fine, fine_root)
      fine_extra <- as.data.frame(fine_fit$pseudotime, stringsAsFactors = FALSE)
      colnames(fine_extra) <- sprintf("lineage_%s", seq_len(ncol(fine_extra)))
      trajectory_write_pseudotime_09(fine_csv, names(fine_fit$pseudotime_vector), fine_fit$pseudotime_vector, label = fine_fit$labels, extra = fine_extra)
    } else {
      write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), fine_csv, row.names = FALSE)
    }

    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "slingshot",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = curve_rds,
      figure_path = fig_png,
      n_cells = ncol(seu),
      status = "ok",
      reason = "",
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      fine_output_path = fine_csv,
      lineage_count = ncol(coarse_fit$pseudotime),
      root_label = root,
      label_var = coarse,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), out_csv, row.names = FALSE)
    write.csv(data.frame(cell_id = character(0), pseudotime = numeric(0)), fine_csv, row.names = FALSE)
    data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "slingshot",
      methods_enabled = "yes",
      input_rds = unit$input_rds[[1]],
      output_path = out_csv,
      extra_path = curve_rds,
      figure_path = fig_png,
      n_cells = suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = "failed",
      reason = conditionMessage(e),
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      fine_output_path = fine_csv,
      lineage_count = 0L,
      root_label = "",
      label_var = normalize_scalar_value(unit$coarse_label_var[[1]]),
      stringsAsFactors = FALSE
    )
  })
  index_rows[[length(index_rows) + 1L]] <- row
  dynamic_outputs[[sprintf("slingshot_%s_pseudotime", trajectory_unit_file_id_09(pair_id, split_value))]] <- build_output_entry(out_csv, "csv", module_name, "Slingshot coarse-label pseudotime", base_dir = cfg$project_root)
  dynamic_outputs[[sprintf("slingshot_%s_fine_pseudotime", trajectory_unit_file_id_09(pair_id, split_value))]] <- build_output_entry(fine_csv, "csv", module_name, "Slingshot fine-cluster pseudotime", base_dir = cfg$project_root)
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(slingshot_index_cols)
write_tsv_local(index_df, cfg$trajectory_slingshot_index_tsv)

outputs <- c(
  list(slingshot_index_tsv = build_output_entry(cfg$trajectory_slingshot_index_tsv, "tsv", module_name, "Slingshot method status by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09d_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09c = cfg$module_09c_manifest_path)
)

message("09d completed. Slingshot index: ", cfg$trajectory_slingshot_index_tsv)
