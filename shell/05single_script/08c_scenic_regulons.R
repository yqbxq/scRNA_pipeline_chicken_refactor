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
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "communication_mapping_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_08.R"))
source_utf8(file.path(.script_dir, "helpers", "scenic_utils.R"))

load_required_packages(c("Seurat", "SCENIC", "RcisTarget", "AUCell", "data.table", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_08()
module_name <- "08c_scenic_regulons"
prepare_dirs_08(cfg)
set.seed(cfg$random_seed)
patch_import_rankings_08()

for (path in c(cfg$module_08a_manifest_path, cfg$module_08b_manifest_path)) {
  if (!file.exists(path)) {
    stop(sprintf("缺少 upstream manifest: %s", path), call. = FALSE)
  }
}
for (path in c(cfg$scenic_db_500bp, cfg$scenic_db_10kb, cfg$scenic_motif_ann)) {
  if (!file.exists(path)) {
    stop(sprintf("Missing SCENIC database file: %s", path), call. = FALSE)
  }
}

export_index <- read_tsv_optional(cfg$scenic_export_index_tsv)
grn_index <- read_tsv_optional(cfg$scenic_grn_index_tsv)
if (nrow(export_index) == 0 || nrow(grn_index) == 0) {
  stop("08c 需要 08a export index 和 08b GRN index", call. = FALSE)
}
tasks <- merge(export_index, grn_index[, c("layer_id", "adjacencies_tsv"), drop = FALSE], by = "layer_id", all.x = FALSE, all.y = FALSE)
if (nrow(tasks) == 0) {
  stop("08c 没有匹配到可运行的 SCENIC layer", call. = FALSE)
}

index_rows <- list()
dynamic_outputs <- list()

motif_annotations <- importAnnotations(cfg$scenic_motif_ann)
db_files <- c("500bp" = cfg$scenic_db_500bp, "10kb" = cfg$scenic_db_10kb)

for (idx in seq_len(nrow(tasks))) {
  task <- tasks[idx, , drop = FALSE]
  layer_id <- task$layer_id[[1]]
  safe_layer <- safe_id_08(layer_id)
  message("08c SCENIC regulons layer: ", layer_id)

  input_rds <- task$object_rds[[1]]
  adj_file <- task$adjacencies_tsv[[1]]
  expr_csv <- task$expr_mat_human_csv[[1]]
  ortholog_csv <- task$ortholog_csv[[1]]
  if (!file.exists(input_rds)) stop(sprintf("Missing layer object: %s", input_rds), call. = FALSE)
  if (!file.exists(adj_file)) stop(sprintf("Missing SCENIC network output: %s", adj_file), call. = FALSE)
  if (!file.exists(expr_csv)) stop(sprintf("Missing SCENIC expression export: %s", expr_csv), call. = FALSE)

  obj <- readRDS(input_rds)
  obj <- maybe_join_layers(obj)
  expr_genes_x_cells <- read_scenic_expr_csv_08(expr_csv)

  common_cells <- intersect(colnames(obj), colnames(expr_genes_x_cells))
  if (length(common_cells) == 0) {
    stop(sprintf("No overlapping cells between layer %s object and SCENIC expression matrix", layer_id), call. = FALSE)
  }
  expr_genes_x_cells <- expr_genes_x_cells[, common_cells, drop = FALSE]
  obj <- subset(obj, cells = common_cells)

  adj <- data.table::fread(adj_file)
  if (!all(c("TF", "target", "importance") %in% colnames(adj))) {
    adj <- data.table::fread(adj_file, sep = "\t")
  }

  modules <- build_scenic_modules_from_adj_08(
    adj,
    rownames(expr_genes_x_cells),
    cfg$scenic_module_top_n,
    cfg$scenic_module_min_genes
  )
  if (length(modules) == 0) {
    stop(sprintf("No valid SCENIC modules were generated for layer %s", layer_id), call. = FALSE)
  }

  paths <- scenic_regulon_paths_08(cfg, layer_id)
  ensure_dir(paths$output_dir)

  motif_results <- lapply(names(db_files), function(db_name) {
    message("Running RcisTarget on database: ", db_name, " layer: ", layer_id)
    rankings <- importRankings(db_files[[db_name]])
    enrichment <- cisTarget(modules, rankings, motif_annotations, nCores = cfg$scenic_threads)
    enrichment$db <- db_name
    as.data.frame(enrichment)
  })
  names(motif_results) <- names(db_files)

  motif_df <- data.table::rbindlist(motif_results, fill = TRUE)
  data.table::fwrite(motif_df, paths$motif_enrichment_tsv, sep = "\t")
  saveRDS(motif_results, paths$motif_enrichment_rds)

  regulons_list <- regulons_from_motif_df_08(
    motif_df,
    rownames(expr_genes_x_cells),
    cfg$scenic_module_top_n,
    cfg$scenic_regulon_min_targets,
    cfg$scenic_nes_threshold
  )
  if (length(regulons_list) == 0) {
    stop(sprintf("No regulons remained after filtering for layer %s", layer_id), call. = FALSE)
  }

  regulon_summary <- tibble::tibble(
    TF = names(regulons_list),
    TargetCount = lengths(regulons_list)
  ) %>%
    dplyr::arrange(dplyr::desc(TargetCount), TF)
  write.csv(regulon_summary, paths$regulons_csv, row.names = FALSE)

  regulon_targets_long <- dplyr::bind_rows(lapply(names(regulons_list), function(tf_name) {
    tibble::tibble(TF = tf_name, Target = regulons_list[[tf_name]])
  }))
  write.csv(regulon_targets_long, paths$regulon_targets_long_csv, row.names = FALSE)
  write_regulons_gmt_08(regulons_list, paths$regulons_gmt)

  cells_rankings <- AUCell_buildRankings(expr_genes_x_cells, plotStats = FALSE, nCores = cfg$scenic_threads)
  auc_max_rank <- max(1L, as.integer(ceiling(nrow(expr_genes_x_cells) * cfg$scenic_auc_max_rank_fraction)))
  cells_auc <- AUCell_calcAUC(regulons_list, cells_rankings, aucMaxRank = auc_max_rank)
  auc_features_x_cells <- getAUC(cells_auc)
  auc_cells_x_regulons <- t(auc_features_x_cells)
  write.csv(as.data.frame(auc_cells_x_regulons), paths$auc_matrix_csv, row.names = TRUE)

  saveRDS(
    list(
      regulons = regulons_list,
      motif_enrichment = motif_results,
      auc_matrix = auc_features_x_cells,
      ortholog_map = if (file.exists(ortholog_csv)) read.csv(ortholog_csv, stringsAsFactors = FALSE) else NULL
    ),
    paths$scenic_core_outputs_rds
  )

  if (identical(layer_id, cfg$panorama_layer_id)) {
    legacy_files <- c(
      motif_enrichment_tsv = "motif_enrichment.tsv",
      motif_enrichment_rds = "motif_enrichment.rds",
      regulons_csv = "regulons.csv",
      regulon_targets_long_csv = "regulon_targets_long.csv",
      regulons_gmt = "regulons.gmt",
      auc_matrix_csv = "auc_matrix.csv",
      scenic_core_outputs_rds = "scenic_core_outputs.rds"
    )
    for (name in names(legacy_files)) {
      file.copy(paths[[name]], file.path(cfg$scenic_output_dir, legacy_files[[name]]), overwrite = TRUE)
    }
  }

  index_rows[[length(index_rows) + 1L]] <- data.frame(
    layer_id = layer_id,
    layer_role = task$layer_role[[1]],
    object_rds = normalize_path_07(input_rds),
    expr_mat_human_csv = normalize_path_07(expr_csv),
    adjacencies_tsv = normalize_path_07(adj_file),
    regulons_csv = normalize_path_07(paths$regulons_csv),
    regulon_targets_long_csv = normalize_path_07(paths$regulon_targets_long_csv),
    regulons_gmt = normalize_path_07(paths$regulons_gmt),
    motif_enrichment_tsv = normalize_path_07(paths$motif_enrichment_tsv),
    auc_matrix_csv = normalize_path_07(paths$auc_matrix_csv),
    scenic_core_outputs_rds = normalize_path_07(paths$scenic_core_outputs_rds),
    regulon_n = length(regulons_list),
    cell_n = ncol(expr_genes_x_cells),
    gene_n = nrow(expr_genes_x_cells),
    stringsAsFactors = FALSE
  )

  dynamic_outputs[[paste0("regulons_csv_", safe_layer)]] <- build_output_entry(paths$regulons_csv, "csv", module_name, "SCENIC regulon summary", base_dir = cfg$project_root, schema = infer_schema_from_df(regulon_summary))
  dynamic_outputs[[paste0("regulon_targets_long_csv_", safe_layer)]] <- build_output_entry(paths$regulon_targets_long_csv, "csv", module_name, "SCENIC regulon target membership", base_dir = cfg$project_root, schema = infer_schema_from_df(regulon_targets_long))
  dynamic_outputs[[paste0("regulons_gmt_", safe_layer)]] <- build_output_entry(paths$regulons_gmt, "gmt", module_name, "SCENIC regulons in GMT format", base_dir = cfg$project_root)
  dynamic_outputs[[paste0("motif_enrichment_tsv_", safe_layer)]] <- build_output_entry(paths$motif_enrichment_tsv, "tsv", module_name, "RcisTarget motif enrichment table", base_dir = cfg$project_root)
  dynamic_outputs[[paste0("auc_matrix_csv_", safe_layer)]] <- build_output_entry(paths$auc_matrix_csv, "csv", module_name, "cells by regulons AUCell matrix", base_dir = cfg$project_root)
  dynamic_outputs[[paste0("scenic_core_outputs_rds_", safe_layer)]] <- build_output_entry(paths$scenic_core_outputs_rds, "rds", module_name, "SCENIC core outputs list", base_dir = cfg$project_root)
}

index_df <- dplyr::bind_rows(index_rows)
write_tsv_local(index_df, cfg$scenic_regulons_index_tsv)

if (file.exists(cfg$module_08c_manifest_path)) {
  unlink(cfg$module_08c_manifest_path)
}
fixed_outputs <- list(
  scenic_regulons_index_tsv = build_output_entry(cfg$scenic_regulons_index_tsv, "tsv", module_name, "one row per SCENIC regulon layer", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))
)
write_manifest_local(
  manifest_path = cfg$module_08c_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_08a = cfg$module_08a_manifest_path,
    module_08b = cfg$module_08b_manifest_path,
    scenic_db_500bp = cfg$scenic_db_500bp,
    scenic_db_10kb = cfg$scenic_db_10kb,
    scenic_motif_ann = cfg$scenic_motif_ann
  ),
  version = cfg$module_version,
  depends_on = list(
    module_08a = cfg$module_08a_manifest_path,
    module_08b = cfg$module_08b_manifest_path,
    scenic_resources_manifest = cfg$scenic_resource_manifest
  )
)

message("08c completed. index: ", cfg$scenic_regulons_index_tsv)
