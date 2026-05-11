.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

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
source_utf8(file.path(.script_dir, "helpers", "decoupler_utils.R"))

load_required_packages(c(
  "Seurat", "dplyr", "tidyr", "tibble", "data.table", "ggplot2",
  "pheatmap", "jsonlite", "Matrix", "decoupleR", "dorothea", "progeny"
))

cfg <- get_single_script_config_08()
module_name <- "08e_decoupler"
prepare_dirs_08(cfg)
set.seed(cfg$random_seed)

filter_network_for_expression_08e <- function(network, expr_genes, min_targets) {
  if (nrow(network) == 0) {
    return(network)
  }
  out <- network[network$target %in% expr_genes, , drop = FALSE]
  if (nrow(out) == 0) {
    return(out)
  }
  target_n <- stats::aggregate(
    target ~ source,
    data = unique(out[, c("source", "target"), drop = FALSE]),
    FUN = length
  )
  keep <- target_n$source[target_n$target >= min_targets]
  out[out$source %in% keep, , drop = FALSE]
}

ortholog_map_path <- resolve_decoupler_ortholog_map_path_08(cfg)
ortholog_map <- load_human_to_chicken_map_08(ortholog_map_path)

resource_names <- c("dorothea", "progeny")
resources <- list()
mapping_rows <- list()
fingerprint_rows <- list()
for (resource in resource_names) {
  raw_info <- load_or_build_decoupler_raw_resource_08(cfg, resource)
  mapped <- map_decoupler_network_to_chicken_08(
    raw_info$data,
    ortholog_map,
    cfg$decoupler_min_targets
  )
  cache_paths <- decoupler_cache_paths_08(cfg, resource)
  mapped_tsv <- write_resource_cache_tsv_08(mapped$network, cache_paths$mapped_tsv)
  resources[[resource]] <- list(
    raw = raw_info$data,
    raw_source = raw_info$source,
    raw_path = raw_info$path,
    mapped_network = mapped$network,
    mapped_tsv = mapped_tsv,
    mapping_summary = mapped$mapping_summary
  )
  mapping_rows[[length(mapping_rows) + 1L]] <- mapped$mapping_summary
  fingerprint_rows[[length(fingerprint_rows) + 1L]] <- resource_fingerprint_row_08(
    resource = resource,
    source = raw_info$source,
    raw_network = raw_info$data,
    mapped_network = mapped$network,
    mapped_tsv = mapped_tsv
  )
}
mapping_summary_df <- dplyr::bind_rows(mapping_rows)
fingerprint_df <- dplyr::bind_rows(fingerprint_rows)
write_tsv_local(mapping_summary_df, cfg$decoupler_network_mapping_summary_tsv)
write_tsv_local(fingerprint_df, cfg$decoupler_resource_fingerprint_tsv)

layers <- selected_regulation_layers_08(cfg)
if (nrow(layers) == 0) {
  stop(sprintf("REGULATION_LAYERS 未匹配到可用注释对象: %s", paste(cfg$regulation_layers, collapse = ",")), call. = FALSE)
}

index_rows <- list()
all_tf_activity <- list()
all_pathway_activity <- list()
dynamic_outputs <- list()

for (idx in seq_len(nrow(layers))) {
  layer_row <- layers[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  safe_layer <- safe_id_08(layer_id)
  message("08e decoupleR layer: ", layer_id)

  obj <- load_comm_layer_object_07(layer_row)
  mat_info <- get_decoupler_expression_matrix_08(obj)
  group_col <- resolve_decoupler_group_col_08(obj, cfg, layer_id)
  groups <- obj@meta.data[colnames(mat_info$matrix), group_col, drop = TRUE]
  group_mat <- switch(
    tolower(cfg$decoupler_activity_level),
    group_average = group_average_expression_08(mat_info$matrix, groups),
    cell_level = as.matrix(mat_info$matrix),
    stop(sprintf("Unsupported DECOUPLER_ACTIVITY_LEVEL: %s", cfg$decoupler_activity_level), call. = FALSE)
  )

  dorothea_network <- filter_network_for_expression_08e(
    resources$dorothea$mapped_network,
    rownames(group_mat),
    cfg$decoupler_min_targets
  )
  progeny_network <- filter_network_for_expression_08e(
    resources$progeny$mapped_network,
    rownames(group_mat),
    cfg$decoupler_min_targets
  )

  tf_result <- run_decoupler_method_08(
    group_mat,
    dorothea_network,
    cfg$decoupler_tf_method,
    cfg$decoupler_min_targets
  )
  tf_activity <- standardize_decoupler_result_08(
    tf_result,
    layer_id = layer_id,
    resource = "dorothea",
    method = cfg$decoupler_tf_method,
    network = dorothea_network
  )

  pathway_result <- run_decoupler_method_08(
    group_mat,
    progeny_network,
    cfg$decoupler_pathway_method,
    cfg$decoupler_min_targets
  )
  pathway_activity <- standardize_decoupler_result_08(
    pathway_result,
    layer_id = layer_id,
    resource = "progeny",
    method = cfg$decoupler_pathway_method,
    network = progeny_network
  )

  paths <- decoupler_paths_08(cfg, layer_id)
  ensure_dir(paths$table_dir)
  ensure_dir(paths$figure_dir)
  ensure_dir(paths$checkpoint_dir)
  write_tsv_local(tf_activity, paths$tf_activity_tsv)
  write_tsv_local(pathway_activity, paths$pathway_activity_tsv)
  write_tsv_local(mapping_summary_df, paths$network_mapping_summary_tsv)
  write_tsv_local(fingerprint_df, paths$resource_fingerprint_tsv)
  write_decoupler_heatmap_08(
    tf_activity,
    paths$tf_activity_heatmap_png,
    title = sprintf("decoupleR TF activity: %s", layer_id),
    top_n = cfg$decoupler_top_tf_n
  )
  write_decoupler_heatmap_08(
    pathway_activity,
    paths$pathway_activity_heatmap_png,
    title = sprintf("decoupleR pathway activity: %s", layer_id),
    top_n = 25L
  )
  saveRDS(
    list(
      layer_id = layer_id,
      object_rds = normalize_path_07(layer_row$object_rds[[1]]),
      activity_level = cfg$decoupler_activity_level,
      group_col = group_col,
      expression_source = mat_info$source,
      group_expression = group_mat,
      tf_activity = tf_activity,
      pathway_activity = pathway_activity,
      dorothea_network = dorothea_network,
      progeny_network = progeny_network,
      network_mapping_summary = mapping_summary_df,
      resource_fingerprint = fingerprint_df
    ),
    paths$decoupler_object_rds
  )

  index_rows[[length(index_rows) + 1L]] <- data.frame(
    layer_id = layer_id,
    layer_role = layer_row$layer_role[[1]],
    object_rds = normalize_path_07(layer_row$object_rds[[1]]),
    activity_level = cfg$decoupler_activity_level,
    group_col = group_col,
    expression_source = mat_info$source,
    group_n = ncol(group_mat),
    gene_n = nrow(group_mat),
    tf_source_n = length(unique(tf_activity$source_id[tf_activity$status == "ok"])),
    pathway_source_n = length(unique(pathway_activity$source_id[pathway_activity$status == "ok"])),
    tf_activity_tsv = normalize_path_07(paths$tf_activity_tsv),
    pathway_activity_tsv = normalize_path_07(paths$pathway_activity_tsv),
    network_mapping_summary_tsv = normalize_path_07(paths$network_mapping_summary_tsv),
    decoupler_object_rds = normalize_path_07(paths$decoupler_object_rds),
    stringsAsFactors = FALSE
  )
  all_tf_activity[[length(all_tf_activity) + 1L]] <- tf_activity
  all_pathway_activity[[length(all_pathway_activity) + 1L]] <- pathway_activity

  dynamic_outputs[[paste0("tf_activity_tsv_", safe_layer)]] <- build_output_entry(
    paths$tf_activity_tsv,
    "tsv",
    module_name,
    "decoupleR TF activity by layer group",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(tf_activity)
  )
  dynamic_outputs[[paste0("pathway_activity_tsv_", safe_layer)]] <- build_output_entry(
    paths$pathway_activity_tsv,
    "tsv",
    module_name,
    "decoupleR pathway activity by layer group",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(pathway_activity)
  )
  dynamic_outputs[[paste0("network_mapping_summary_tsv_", safe_layer)]] <- build_output_entry(
    paths$network_mapping_summary_tsv,
    "tsv",
    module_name,
    "human prior network target mapping to chicken symbols",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(mapping_summary_df)
  )
  dynamic_outputs[[paste0("tf_activity_heatmap_png_", safe_layer)]] <- build_output_entry(
    paths$tf_activity_heatmap_png,
    "png",
    module_name,
    "decoupleR TF activity heatmap",
    base_dir = cfg$project_root
  )
  dynamic_outputs[[paste0("pathway_activity_heatmap_png_", safe_layer)]] <- build_output_entry(
    paths$pathway_activity_heatmap_png,
    "png",
    module_name,
    "decoupleR pathway activity heatmap",
    base_dir = cfg$project_root
  )
  dynamic_outputs[[paste0("decoupler_object_rds_", safe_layer)]] <- build_output_entry(
    paths$decoupler_object_rds,
    "rds",
    module_name,
    "decoupleR group matrix, resources and activity outputs",
    base_dir = cfg$project_root
  )
}

index_df <- dplyr::bind_rows(index_rows)
tf_activity_df <- if (length(all_tf_activity) > 0) dplyr::bind_rows(all_tf_activity) else empty_decoupler_activity_08()
pathway_activity_df <- if (length(all_pathway_activity) > 0) dplyr::bind_rows(all_pathway_activity) else empty_decoupler_activity_08()
write_tsv_local(index_df, cfg$decoupler_index_tsv)
write_tsv_local(tf_activity_df, cfg$decoupler_tf_activity_tsv)
write_tsv_local(pathway_activity_df, cfg$decoupler_pathway_activity_tsv)

if (file.exists(cfg$module_08e_manifest_path)) {
  unlink(cfg$module_08e_manifest_path)
}
fixed_outputs <- list(
  decoupler_index_tsv = build_output_entry(
    cfg$decoupler_index_tsv,
    "tsv",
    module_name,
    "one row per decoupleR layer",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(index_df)
  ),
  tf_activity_tsv = build_output_entry(
    cfg$decoupler_tf_activity_tsv,
    "tsv",
    module_name,
    "aggregate decoupleR TF activity across layers",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(tf_activity_df)
  ),
  pathway_activity_tsv = build_output_entry(
    cfg$decoupler_pathway_activity_tsv,
    "tsv",
    module_name,
    "aggregate decoupleR pathway activity across layers",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(pathway_activity_df)
  ),
  network_mapping_summary_tsv = build_output_entry(
    cfg$decoupler_network_mapping_summary_tsv,
    "tsv",
    module_name,
    "human prior network target mapping to chicken symbols",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(mapping_summary_df)
  ),
  resource_fingerprint_tsv = build_output_entry(
    cfg$decoupler_resource_fingerprint_tsv,
    "tsv",
    module_name,
    "decoupleR resource cache and package fingerprint",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(fingerprint_df)
  )
)
write_manifest_local(
  manifest_path = cfg$module_08e_manifest_path,
  new_outputs = c(fixed_outputs, dynamic_outputs),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    ortholog_map = ortholog_map_path,
    regulation_layers = cfg$regulation_layers,
    module_00 = cfg$ortholog_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path,
    decoupler_resource_dir = cfg$decoupler_resource_dir,
    dorothea_resource_rds = resources$dorothea$raw_path,
    progeny_resource_rds = resources$progeny$raw_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_00 = cfg$ortholog_manifest_path,
    module_03d = cfg$module_03d_manifest_path,
    module_04b = cfg$module_04b_manifest_path
  )
)

message("08e completed. index: ", cfg$decoupler_index_tsv)
