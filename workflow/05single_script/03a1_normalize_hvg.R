#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_03.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03a1_normalize_hvg"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

resolve_cell_cycle_genes <- function(cfg) {
  manifest <- read_manifest_local(cfg$module_00c_manifest_path)
  cc_path <- tryCatch(
    resolve_output_local(manifest, "cc_genes"),
    error = function(e) tryCatch(resolve_output_local(manifest, "chicken_cc_genes"), error = function(e2) "")
  )
  if (!nzchar(cc_path) || !file.exists(cc_path)) {
    stop(sprintf("vars_to_regress 请求 S.Score/G2M.Score，但缺少 00c cell-cycle RDS: %s", cc_path), call. = FALSE)
  }
  readRDS(cc_path)
}

match_features_case_insensitive <- function(features, genes) {
  genes <- unique(genes[nzchar(genes)])
  exact <- intersect(genes, features)
  if (length(exact) > 0) {
    return(exact)
  }
  feature_map <- setNames(features, toupper(features))
  matched <- unname(feature_map[toupper(genes)])
  unique(matched[!is.na(matched) & nzchar(matched)])
}

ensure_cell_cycle_scores <- function(seu, cfg, vars_to_regress) {
  requested <- intersect(vars_to_regress, c("S.Score", "G2M.Score"))
  missing_scores <- setdiff(requested, colnames(seu@meta.data))
  if (length(missing_scores) == 0) {
    return(seu)
  }
  cc_genes <- resolve_cell_cycle_genes(cfg)
  s_features <- match_features_case_insensitive(rownames(seu), cc_genes$s.genes %||% character(0))
  g2m_features <- match_features_case_insensitive(rownames(seu), cc_genes$g2m.genes %||% character(0))
  if (length(s_features) == 0 || length(g2m_features) == 0) {
    stop("无法用 00c cell-cycle RDS 在对象 rownames 中匹配 S/G2M genes，不能生成 S.Score/G2M.Score。", call. = FALSE)
  }
  Seurat::CellCycleScoring(
    seu,
    s.features = s_features,
    g2m.features = g2m_features,
    set.ident = FALSE
  )
}

manifest_02b2 <- read_manifest_local(cfg$module_02b2_manifest_path)
post_qc_rds <- resolve_output_local(manifest_02b2, "post_qc_object")
if (!file.exists(post_qc_rds)) {
  stop(sprintf("缺少 02b2 post_qc_object: %s", post_qc_rds), call. = FALSE)
}

manifest_02c <- read_manifest_local(cfg$module_02c_manifest_path)
post_qc_report <- resolve_output_local(manifest_02c, "report")
if (!file.exists(post_qc_report)) {
  stop(sprintf("缺少 02c post-QC EDA report: %s", post_qc_report), call. = FALSE)
}

layer_df <- validate_layer_config(read_object_layer_config(cfg))
enabled_layers <- filter_enabled_layers(layer_df)
panorama_row <- enabled_layers[enabled_layers$layer_role == "panorama", , drop = FALSE]
panorama_spec <- layer_config_row_to_spec(panorama_row)
norm_methods <- panorama_spec$normalization_methods

output_entries <- list()
output_paths <- list()
for (norm_method in norm_methods) {
  message("03a1 normalization: ", norm_method)
  seu <- readRDS(post_qc_rds)
  seu <- maybe_join_layers(seu)
  seu <- strip_reduction_state(seu)
  seu <- ensure_cell_cycle_scores(seu, cfg, panorama_spec$vars_to_regress)
  seu <- normalize_layer(seu, norm_method, panorama_spec, cfg)
  seu@misc$panorama_layer_spec <- panorama_spec
  seu@misc$normalization_method <- norm_method
  out_rds <- file.path(cfg$panorama_normalized_dir, sprintf("%s__%s.rds", panorama_spec$layer_id, norm_method))
  saveRDS(seu, out_rds)
  output_paths[[norm_method]] <- out_rds
  output_entries[[paste0("normalized_", norm_method)]] <- build_output_entry(
    out_rds,
    "rds",
    module_name,
    sprintf("panorama normalized object using %s", norm_method),
    base_dir = cfg$project_root
  )
}

write_manifest_local(
  manifest_path = cfg$module_03a1_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    post_qc_object = post_qc_rds,
    post_qc_report = post_qc_report,
    object_layer_config = cfg$object_layer_config_file,
    cell_cycle_manifest = cfg$module_00c_manifest_path
  ),
  version = cfg$module_version,
  depends_on = list(
    module_02b2 = cfg$module_02b2_manifest_path,
    module_02c = cfg$module_02c_manifest_path
  )
)

message("03a1 完成。normalization outputs: ", paste(unlist(output_paths, use.names = FALSE), collapse = ", "))
