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
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_08.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix"))

cfg <- get_single_script_config_08()
module_name <- "08a_scenic_export"
prepare_dirs_08(cfg)
set.seed(cfg$random_seed)

resolve_ortholog_map_path_08a <- function(cfg) {
  env_path <- normalize_scalar_value(Sys.getenv("SCENIC_ORTHOLOG_MAP_FILE", ""))
  if (nzchar(env_path) && file.exists(env_path)) {
    return(env_path)
  }
  legacy_pipeline_map <- file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs_for_pipeline.csv")
  if (file.exists(legacy_pipeline_map)) {
    return(legacy_pipeline_map)
  }
  map_path <- resolve_manifest_output_optional_07(
    cfg$ortholog_manifest_path,
    c("human_best", "chicken_human_best_csv", "human_best_csv")
  )
  if (nzchar(map_path) && file.exists(map_path)) {
    return(map_path)
  }
  file.path(cfg$ortholog_cache_dir, "chicken_human_orthologs.csv")
}

load_scenic_orthologs_08a <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("缺少 SCENIC ortholog map: %s", path), call. = FALSE)
  }
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  chicken_col <- pick_first_col_07(raw, c("external_gene_name", "chicken_symbol", "gene_name"))
  human_col <- pick_first_col_07(raw, c("target_gene_name", "hsapiens_homolog_associated_gene_name", "human_symbol"))
  if (!nzchar(chicken_col) || !nzchar(human_col)) {
    stop(
      sprintf(
        "SCENIC ortholog map 缺少必需列；需要 chicken symbol 和 human symbol。文件: %s",
        path
      ),
      call. = FALSE
    )
  }
  if ("target_species" %in% colnames(raw)) {
    human_rows <- tolower(trimws(as.character(raw$target_species))) == "human"
    if (any(human_rows)) {
      raw <- raw[human_rows, , drop = FALSE]
    }
  }
  orthologs <- data.frame(
    external_gene_name = trimws(as.character(raw[[chicken_col]])),
    hsapiens_homolog_associated_gene_name = trimws(as.character(raw[[human_col]])),
    stringsAsFactors = FALSE
  )
  orthologs <- orthologs[
    nzchar(orthologs$external_gene_name) &
      !is.na(orthologs$external_gene_name) &
      nzchar(orthologs$hsapiens_homolog_associated_gene_name) &
      !is.na(orthologs$hsapiens_homolog_associated_gene_name),
    ,
    drop = FALSE
  ]
  orthologs <- orthologs[!duplicated(orthologs$external_gene_name), , drop = FALSE]
  rownames(orthologs) <- NULL
  orthologs
}

selected_regulation_layers_08a <- function(cfg) {
  layers <- communication_layer_status_07(cfg)
  requested <- unique(trimws(as.character(cfg$regulation_layers)))
  requested <- requested[nzchar(requested)]
  if (length(requested) == 0) {
    requested <- cfg$panorama_layer_id
  }
  if (!any(tolower(requested) %in% c("all", "*"))) {
    layers <- layers[layers$layer_id %in% requested, , drop = FALSE]
  }
  layers
}

ortholog_map_path <- resolve_ortholog_map_path_08a(cfg)
orthologs <- load_scenic_orthologs_08a(ortholog_map_path)
map_vec <- setNames(orthologs$hsapiens_homolog_associated_gene_name, orthologs$external_gene_name)

layers <- selected_regulation_layers_08a(cfg)
if (nrow(layers) == 0) {
  stop(
    sprintf(
      "REGULATION_LAYERS 未匹配到可用注释对象: %s",
      paste(cfg$regulation_layers, collapse = ",")
    ),
    call. = FALSE
  )
}

index_rows <- list()
mapping_rows <- list()
triage_rows <- list()
dynamic_outputs <- list()

for (idx in seq_len(nrow(layers))) {
  layer_row <- layers[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  safe_layer <- safe_id_08(layer_id)
  message("08a SCENIC export layer: ", layer_id)

  obj <- load_comm_layer_object_07(layer_row)
  counts <- get_assay_matrix(obj, type = "counts")
  gene_ids <- rownames(counts)
  mapped_idx <- which(gene_ids %in% names(map_vec))
  if (length(mapped_idx) == 0) {
    stop(sprintf("layer %s 没有任何基因命中 human ortholog map", layer_id), call. = FALSE)
  }

  counts_human <- counts[mapped_idx, , drop = FALSE]
  rownames(counts_human) <- map_vec[rownames(counts_human)]
  counts_human_agg <- rowsum(as.matrix(counts_human), group = rownames(counts_human))
  counts_export <- t(counts_human_agg)

  paths <- scenic_export_paths_08(cfg, layer_id)
  ensure_dir(paths$input_dir)
  write.csv(orthologs, paths$ortholog_csv, row.names = FALSE)
  write.csv(counts_export, paths$expr_mat_human_csv, row.names = TRUE)

  if (identical(layer_id, cfg$panorama_layer_id)) {
    write.csv(orthologs, file.path(cfg$scenic_input_dir, "chicken_to_human_orthologs.csv"), row.names = FALSE)
    write.csv(counts_export, file.path(cfg$scenic_input_dir, "expr_mat_human.csv"), row.names = TRUE)
  }

  mapping_rate <- length(unique(gene_ids[mapped_idx])) / length(unique(gene_ids))
  mapping_rows[[length(mapping_rows) + 1L]] <- data.frame(
    layer_id = layer_id,
    object_rds = normalize_path_07(layer_row$object_rds[[1]]),
    chicken_gene_n = length(unique(gene_ids)),
    mapped_chicken_gene_n = length(unique(gene_ids[mapped_idx])),
    human_gene_n = nrow(counts_human_agg),
    cell_n = ncol(counts),
    mapping_rate = mapping_rate,
    ortholog_map = normalize_path_07(ortholog_map_path),
    stringsAsFactors = FALSE
  )

  if (is.finite(mapping_rate) && mapping_rate < cfg$scenic_ortholog_min_coverage) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = layer_id,
      severity = "warning",
      signal_id = "scenic_ortholog_coverage_low",
      evidence = sprintf("SCENIC human ortholog coverage %.3f is below %.3f", mapping_rate, cfg$scenic_ortholog_min_coverage),
      recommended_action = "确认对象基因名是否为 chicken gene symbol；必要时重建 00 ortholog cache。",
      manual_review_required = "yes"
    )
  }

  index_rows[[length(index_rows) + 1L]] <- data.frame(
    layer_id = layer_id,
    layer_role = layer_row$layer_role[[1]],
    object_rds = normalize_path_07(layer_row$object_rds[[1]]),
    expr_mat_human_csv = normalize_path_07(paths$expr_mat_human_csv),
    ortholog_csv = normalize_path_07(paths$ortholog_csv),
    cell_n = ncol(counts),
    chicken_gene_n = length(unique(gene_ids)),
    mapped_chicken_gene_n = length(unique(gene_ids[mapped_idx])),
    human_gene_n = nrow(counts_human_agg),
    stringsAsFactors = FALSE
  )
  dynamic_outputs[[paste0("expr_mat_human_csv_", safe_layer)]] <- build_output_entry(paths$expr_mat_human_csv, "csv", module_name, "cells by human-symbol genes SCENIC expression matrix", base_dir = cfg$project_root)
  dynamic_outputs[[paste0("ortholog_csv_", safe_layer)]] <- build_output_entry(paths$ortholog_csv, "csv", module_name, "chicken to human symbol map used for SCENIC export", base_dir = cfg$project_root)
}

index_df <- dplyr::bind_rows(index_rows)
mapping_df <- dplyr::bind_rows(mapping_rows)
triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)

write_tsv_local(index_df, cfg$scenic_export_index_tsv)
write_tsv_local(mapping_df, cfg$scenic_export_mapping_summary_tsv)
write_tsv_local(triage_df, cfg$scenic_export_triage_tsv)

if (file.exists(cfg$module_08a_manifest_path)) {
  unlink(cfg$module_08a_manifest_path)
}
fixed_outputs <- list(
  scenic_export_index_tsv = build_output_entry(cfg$scenic_export_index_tsv, "tsv", module_name, "one row per SCENIC-exported layer", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df)),
  mapping_summary_tsv = build_output_entry(cfg$scenic_export_mapping_summary_tsv, "tsv", module_name, "SCENIC ortholog mapping summary by layer", base_dir = cfg$project_root, schema = infer_schema_from_df(mapping_df)),
  triage_tsv = build_output_entry(cfg$scenic_export_triage_tsv, "tsv", module_name, "SCENIC export triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
)
write_manifest_local(
  manifest_path = cfg$module_08a_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    ortholog_map = ortholog_map_path,
    regulation_layers = cfg$regulation_layers,
    module_00 = cfg$ortholog_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_00 = cfg$ortholog_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("08a completed. index: ", cfg$scenic_export_index_tsv)
