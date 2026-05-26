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
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "reduction_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_03()
module_name <- "03a3_compute_umap"
prepare_dirs_03(cfg)

umap_seed <- as.integer(Sys.getenv("UMAP_SEED", as.character(cfg$random_seed)))
if (!is.finite(umap_seed)) {
  umap_seed <- cfg$random_seed
}
set.seed(umap_seed)

panorama_spec <- panorama_layer_spec(cfg)
manifest_03a2 <- read_manifest_local(cfg$module_03a2_manifest_path)
candidate_index_tsv <- resolve_output_local(manifest_03a2, "candidate_index_tsv")
candidate_index <- read_tsv_local(candidate_index_tsv)
if (nrow(candidate_index) == 0) {
  stop("03a2 candidate index is empty; cannot compute UMAP", call. = FALSE)
}

candidate_rows <- list()
output_entries <- list()

for (i in seq_len(nrow(candidate_index))) {
  row <- candidate_index[i, , drop = FALSE]
  candidate_key <- normalize_scalar_value(row$candidate_key[[1]], sprintf("candidate_%03d", i))
  candidate_rds <- normalize_scalar_value(row$out_rds[[1]])
  reduction_name <- normalize_scalar_value(row$reduction_name[[1]], "pca")
  umap_name <- if ("umap_name" %in% colnames(row)) {
    normalize_scalar_value(row$umap_name[[1]], "")
  } else {
    ""
  }
  if (!nzchar(umap_name)) {
    umap_name <- sprintf("umap_%s_%s", row$integration[[1]], row$normalization[[1]])
  }

  message("03a3 UMAP candidate: ", candidate_key, " / ", reduction_name, " -> ", umap_name)
  seu <- readRDS(candidate_rds)
  seu <- maybe_join_layers(seu)
  seu <- run_umap_only(seu, panorama_spec, reduction_name = reduction_name, umap_name = umap_name, seed = umap_seed)
  seu@misc$selected_umap_candidate <- umap_name

  out_rds <- file.path(cfg$panorama_umap_dir, sprintf("%s__umap.rds", candidate_key))
  saveRDS(seu, out_rds)
  output_key <- paste0(candidate_key, "_umap")
  output_entries[[output_key]] <- build_output_entry(
    out_rds,
    "rds",
    module_name,
    sprintf("panorama candidate with UMAP candidate_key=%s", candidate_key),
    base_dir = cfg$project_root
  )

  candidate_rows[[length(candidate_rows) + 1L]] <- data.frame(
    layer_id = normalize_scalar_value(row$layer_id[[1]], panorama_spec$layer_id),
    normalization = normalize_scalar_value(row$normalization[[1]], ""),
    integration = normalize_scalar_value(row$integration[[1]], ""),
    candidate_key = output_key,
    source_candidate_key = candidate_key,
    reduction_name = reduction_name,
    umap_name = umap_name,
    umap_n_neighbors = as.integer(Sys.getenv("UMAP_N_NEIGHBORS", "30")),
    umap_min_dist = as.numeric(Sys.getenv("UMAP_MIN_DIST", "0.3")),
    umap_spread = as.numeric(Sys.getenv("UMAP_SPREAD", "1.0")),
    umap_metric = Sys.getenv("UMAP_METRIC", "cosine"),
    umap_seed = umap_seed,
    out_rds = normalizePath(out_rds, winslash = "/", mustWork = FALSE),
    source_rds = normalizePath(candidate_rds, winslash = "/", mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

umap_index <- dplyr::bind_rows(candidate_rows)
umap_index_tsv <- file.path(cfg$integration_table_dir_layer, "umap_candidates.tsv")
write_tsv_local(umap_index, umap_index_tsv)
output_entries$candidate_umap_index_tsv <- build_output_entry(
  umap_index_tsv,
  "tsv",
  module_name,
  "one row per normalization/integration UMAP candidate",
  base_dir = cfg$project_root,
  schema = infer_schema_from_df(umap_index)
)

write_manifest_local(
  manifest_path = cfg$module_03a3_manifest_path,
  new_outputs = output_entries,
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    module_03a2_manifest = cfg$module_03a2_manifest_path,
    candidate_index_tsv = candidate_index_tsv
  ),
  version = cfg$module_03a3_version,
  depends_on = list(module_03a2 = cfg$module_03a2_manifest_path)
)

message("03a3 完成。UMAP 候选索引: ", umap_index_tsv)
