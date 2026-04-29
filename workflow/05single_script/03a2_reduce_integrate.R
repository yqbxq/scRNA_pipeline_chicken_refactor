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
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "triage_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "integration_dispatch.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03a2_reduce_integrate"
prepare_dirs_03(cfg)
set.seed(cfg$random_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03a1 <- read_manifest_local(cfg$module_03a1_manifest_path)

candidate_rows <- list()
triage_rows <- list()
output_entries <- list()

append_triage <- function(severity, signal_id, evidence) {
  triage_rows[[length(triage_rows) + 1]] <<- make_triage_row(
    sample_id = panorama_spec$layer_id,
    severity = severity,
    signal_id = signal_id,
    evidence = evidence
  )
}

for (norm_method in panorama_spec$normalization_methods) {
  norm_key <- paste0("normalized_", norm_method)
  norm_rds <- resolve_output_local(manifest_03a1, norm_key)
  message("03a2 PCA/UMAP for normalization: ", norm_method)
  seu <- readRDS(norm_rds)
  seu <- maybe_join_layers(seu)
  assay_name <- if (identical(norm_method, "sct") && "SCT" %in% Assays(seu)) "SCT" else "RNA"
  seu <- reduce_pca_umap(
    seu,
    panorama_spec,
    assay = assay_name,
    reduction_key_prefix = sprintf("PCA%s_", toupper(norm_method)),
    umap_name = sprintf("umap_rna_%s", norm_method)
  )

  for (integration_mode in panorama_spec$integration_mode) {
    message("03a2 integration candidate: ", norm_method, " / ", integration_mode)
    result <- dispatch_integration(
      seu,
      integration_mode,
      group_var = "orig.ident",
      layer_spec = panorama_spec,
      normalization_method = norm_method
    )
    seu_int <- result$object
    dims <- usable_reduction_dims(seu_int, panorama_spec$pca_dims, result$reduction_name)
    if (length(dims) >= 2) {
      seu_int <- RunUMAP(
        seu_int,
        reduction = result$reduction_name,
        dims = dims,
        reduction.name = result$umap_name,
        reduction.key = paste0(gsub("[^A-Za-z0-9]", "", toupper(result$umap_name)), "_"),
        verbose = FALSE
      )
    }
    seu_int@misc$panorama_layer_spec <- panorama_spec
    seu_int@misc$normalization_method <- norm_method
    seu_int@misc$integration_mode <- integration_mode
    seu_int@misc$selected_reduction_candidate <- result$reduction_name
    seu_int@misc$selected_umap_candidate <- result$umap_name

    out_rds <- file.path(cfg$panorama_reduction_dir, sprintf("%s__%s__%s.rds", panorama_spec$layer_id, norm_method, integration_mode))
    saveRDS(seu_int, out_rds)
    candidate_key <- paste0("candidate_", norm_method, "_", integration_mode)
    output_entries[[candidate_key]] <- build_output_entry(
      out_rds,
      "rds",
      module_name,
      sprintf("panorama candidate normalization=%s integration=%s", norm_method, integration_mode),
      base_dir = cfg$project_root
    )

    downgrade_reason <- normalize_scalar_value(result$diagnostics$downgrade_reason)
    if (nzchar(downgrade_reason)) {
      append_triage(
        "medium",
        "integration_method_unavailable",
        sprintf("normalization=%s; integration=%s; downgrade_reason=%s", norm_method, integration_mode, downgrade_reason)
      )
    }
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      layer_id = panorama_spec$layer_id,
      normalization = norm_method,
      integration = integration_mode,
      candidate_key = candidate_key,
      reduction_name = result$reduction_name,
      umap_name = result$umap_name,
      runtime_sec = as.numeric(result$diagnostics$runtime_sec %||% 0),
      downgrade_reason = downgrade_reason,
      out_rds = normalizePath(out_rds, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  }
}

candidate_index <- dplyr::bind_rows(candidate_rows)
candidate_index_tsv <- file.path(cfg$integration_table_dir_layer, "reduction_candidates.tsv")
write_tsv_local(candidate_index, candidate_index_tsv)

triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
triage_tsv <- file.path(cfg$integration_table_dir_layer, "integration_method_triage.tsv")
write_tsv_local(triage_df, triage_tsv)

output_entries$candidate_index_tsv <- build_output_entry(
  candidate_index_tsv,
  "tsv",
  module_name,
  "one row per normalization/integration candidate",
  base_dir = cfg$project_root,
  schema = infer_schema_from_df(candidate_index)
)
output_entries$integration_method_triage_tsv <- build_output_entry(
  triage_tsv,
  "tsv",
  module_name,
  "one row per integration dispatch triage signal",
  base_dir = cfg$project_root,
  schema = infer_schema_from_df(triage_df)
)

write_manifest_local(
  manifest_path = cfg$module_03a2_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03a1_manifest = cfg$module_03a1_manifest_path,
    object_layer_config = cfg$object_layer_config_file
  ),
  version = cfg$module_version,
  depends_on = list(module_03a1 = cfg$module_03a1_manifest_path)
)

message("03a2 完成。候选索引: ", candidate_index_tsv)
