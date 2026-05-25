#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_enrichment_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_06e_niche"
prepare_dirs_spatial(cfg)

niche_manifest_cols_st <- function() {
  c(
    "status", "reason", "input_rds", "neighborhood_manifest", "deconv_manifest",
    "deconv_method", "recommended_method", "score_tsv", "label_tsv",
    "spatial_panorama_niched_rds", "spot_n", "niche_n", "k", "k_neighbors"
  )
}

resolve_panorama_input_st <- function(cfg) {
  candidates <- c(
    st06_manifest_output(cfg$module_04c_subcluster_eda_manifest_path, c("spatial_panorama_subannotated", "subannotated_object", "panorama_rds"), cfg),
    st06_manifest_output(cfg$module_03a_region_annotation_eda_manifest_path, c("spatial_panorama_annotated", "annotated_object", "panorama_rds"), cfg),
    cfg$spatial_panorama_subannotated_rds,
    cfg$spatial_panorama_annotated_rds,
    cfg$spatial_panorama_clustered_rds
  )
  hit <- candidates[nzchar(candidates) & file.exists(candidates)][1]
  if (is.na(hit)) "" else hit
}

resolve_neighborhood_manifest_st <- function(cfg) {
  path <- st06_manifest_output(cfg$module_06d_neighborhood_manifest_path, "neighborhood_manifest", cfg)
  if (!nzchar(path)) {
    fallback <- file.path(cfg$spatial_neighborhood_table_dir, "neighborhood_manifest.tsv")
    if (file.exists(fallback)) {
      path <- fallback
    }
  }
  path
}

read_recommended_method_st <- function(cfg) {
  if (file.exists(cfg$spatial_recommended_deconv_method_file)) {
    value <- trimws(readLines(cfg$spatial_recommended_deconv_method_file, warn = FALSE))
    value <- value[nzchar(value)]
    if (length(value) > 0) {
      return(value[[1]])
    }
  }
  cfg$spatial_deconv_primary
}

collect_deconv_manifests_st <- function(cfg) {
  specs <- list(
    rctd = list(path = cfg$module_07a_deconvolution_rctd_manifest_path, key = "rctd_manifest", fallback = file.path(cfg$spatial_rctd_table_dir, "rctd_manifest.tsv")),
    transfer = list(path = cfg$module_07b_deconvolution_transfer_manifest_path, key = "transfer_manifest", fallback = file.path(cfg$spatial_transfer_table_dir, "transfer_manifest.tsv")),
    card = list(path = cfg$module_07c_deconvolution_card_manifest_path, key = "card_manifest", fallback = file.path(cfg$spatial_card_table_dir, "card_manifest.tsv")),
    cell2location = list(path = cfg$module_07d_deconvolution_cell2location_manifest_path, key = "cell2location_manifest", fallback = file.path(cfg$spatial_c2l_table_dir, "cell2location_manifest.tsv"))
  )
  rows <- list()
  for (method in names(specs)) {
    manifest_tsv <- st06_manifest_output(specs[[method]]$path, specs[[method]]$key, cfg)
    if (!nzchar(manifest_tsv) && file.exists(specs[[method]]$fallback)) {
      manifest_tsv <- specs[[method]]$fallback
    }
    df <- st06_read_tsv(manifest_tsv)
    if (nrow(df) == 0) {
      next
    }
    df$deconv_method <- method
    df$method_manifest_tsv <- manifest_tsv
    rows[[length(rows) + 1L]] <- df
  }
  if (length(rows) == 0) {
    return(st06_empty_df(c("status", "reason", "deconv_method", "method_manifest_tsv", "proportion_wide_tsv")))
  }
  cols <- unique(unlist(lapply(rows, colnames), use.names = FALSE))
  rows <- lapply(rows, function(df) {
    for (col in setdiff(cols, colnames(df))) {
      df[[col]] <- ""
    }
    df[, cols, drop = FALSE]
  })
  do.call(rbind, rows)
}

write_niche_skip_st <- function(status, reason, input_rds = "", neighborhood_manifest = "", deconv_manifest = "", deconv_method = "", recommended_method = "") {
  data.frame(
    status = status,
    reason = reason,
    input_rds = input_rds,
    neighborhood_manifest = neighborhood_manifest,
    deconv_manifest = deconv_manifest,
    deconv_method = deconv_method,
    recommended_method = recommended_method,
    score_tsv = "",
    label_tsv = "",
    spatial_panorama_niched_rds = "",
    spot_n = 0L,
    niche_n = 0L,
    k = cfg$spatial_niche_k,
    k_neighbors = cfg$spatial_niche_k_neighbors,
    stringsAsFactors = FALSE
  )
}

run_niche_kmeans_st <- function(obj, deconv_row, cfg) {
  score_path <- st06_abs_path(st06_scalar(deconv_row$proportion_wide_tsv, st06_scalar(deconv_row$proportion_tsv)), cfg)
  scores <- st06_read_tsv(score_path)
  spot_col <- st06_pick_col(scores, c("spot_id", "barcode", "cell", "cell_id", "sample_id"))
  if (!nzchar(score_path) || !file.exists(score_path) || nrow(scores) == 0 || !nzchar(spot_col)) {
    return(list(status = "skipped_no_deconv_scores", reason = "selected deconvolution row has no readable spot-level proportion table"))
  }
  feature_cols <- setdiff(colnames(scores), c(spot_col, "section", "sample_id", "deconv_id", "tool", "method"))
  numeric_cols <- feature_cols[vapply(scores[feature_cols], function(x) any(is.finite(suppressWarnings(as.numeric(x)))), logical(1))]
  if (length(numeric_cols) < 2) {
    return(list(status = "skipped_too_few_features", reason = "need at least two numeric deconvolution proportion features"))
  }
  spot_ids <- as.character(scores[[spot_col]])
  keep <- spot_ids %in% rownames(obj@meta.data)
  scores <- scores[keep, , drop = FALSE]
  spot_ids <- spot_ids[keep]
  if (nrow(scores) < cfg$spatial_niche_k) {
    return(list(status = "skipped_too_few_spots", reason = sprintf("%s matched spots; minimum k is %s", nrow(scores), cfg$spatial_niche_k)))
  }
  mat <- as.matrix(data.frame(lapply(scores[numeric_cols], function(x) suppressWarnings(as.numeric(x)))))
  rownames(mat) <- spot_ids
  mat[!is.finite(mat)] <- 0
  k <- min(cfg$spatial_niche_k, nrow(unique(mat)), nrow(mat))
  if (k < 2) {
    return(list(status = "skipped_too_few_features", reason = "fewer than two unique deconvolution profiles"))
  }
  set.seed(cfg$random_seed)
  km <- stats::kmeans(scale(mat), centers = k, nstart = 20)
  niche_label <- paste0("niche_", sprintf("%02d", km$cluster))

  score_tsv <- file.path(cfg$spatial_niche_table_dir, "niche_scores.tsv")
  label_tsv <- file.path(cfg$spatial_niche_table_dir, "niche_labels.tsv")
  score_df <- data.frame(spot_id = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
  label_df <- data.frame(
    spot_id = rownames(mat),
    spatial_niche = niche_label,
    spatial_niche_cluster = km$cluster,
    spatial_niche_method = "deconv_proportion_kmeans",
    deconv_method = st06_scalar(deconv_row$deconv_method),
    stringsAsFactors = FALSE
  )
  st06_write_tsv(score_df, score_tsv)
  st06_write_tsv(label_df, label_tsv)

  meta <- obj@meta.data
  for (col in c("spatial_niche", "spatial_niche_cluster", "spatial_niche_method", "spatial_niche_deconv_method")) {
    if (!col %in% colnames(meta)) {
      meta[[col]] <- NA
    }
  }
  idx <- match(label_df$spot_id, rownames(meta))
  meta$spatial_niche[idx] <- label_df$spatial_niche
  meta$spatial_niche_cluster[idx] <- label_df$spatial_niche_cluster
  meta$spatial_niche_method[idx] <- label_df$spatial_niche_method
  meta$spatial_niche_deconv_method[idx] <- label_df$deconv_method
  obj@meta.data <- meta
  saveRDS(obj, cfg$spatial_panorama_niched_rds)

  list(
    status = "ok",
    reason = "",
    score_tsv = score_tsv,
    label_tsv = label_tsv,
    spatial_panorama_niched_rds = cfg$spatial_panorama_niched_rds,
    spot_n = nrow(label_df),
    niche_n = length(unique(label_df$spatial_niche))
  )
}

panorama_input <- resolve_panorama_input_st(cfg)
neighborhood_manifest <- resolve_neighborhood_manifest_st(cfg)
neighborhood <- st06_read_tsv(neighborhood_manifest)
recommended_method <- read_recommended_method_st(cfg)
deconv <- collect_deconv_manifests_st(cfg)
ok_deconv <- deconv[deconv$status == "ok", , drop = FALSE]
if (nrow(ok_deconv) > 0 && "deconv_method" %in% colnames(ok_deconv) && recommended_method %in% ok_deconv$deconv_method) {
  ok_deconv <- ok_deconv[ok_deconv$deconv_method == recommended_method, , drop = FALSE]
}

if (!nzchar(panorama_input)) {
  manifest_df <- write_niche_skip_st("skipped_no_upstream_manifest", "no annotated/subannotated spatial panorama RDS was available", neighborhood_manifest = neighborhood_manifest, recommended_method = recommended_method)
} else if (nrow(neighborhood) == 0) {
  manifest_df <- write_niche_skip_st("skipped_no_neighborhood_manifest", "06d neighborhood manifest is unavailable", panorama_input, neighborhood_manifest, recommended_method = recommended_method)
} else if (!identical(st06_scalar(neighborhood$status), "ok")) {
  manifest_df <- write_niche_skip_st("skipped_no_neighborhood_manifest", sprintf("06d neighborhood status is %s", st06_scalar(neighborhood$status)), panorama_input, neighborhood_manifest, recommended_method = recommended_method)
} else if (nrow(ok_deconv) == 0) {
  manifest_df <- write_niche_skip_st("skipped_no_deconv_manifest", "no 07 deconvolution method manifest has status ok; P-06-E remains post-07 gated", panorama_input, neighborhood_manifest, recommended_method = recommended_method)
} else if (!requireNamespace("Seurat", quietly = TRUE)) {
  manifest_df <- write_niche_skip_st("skipped_no_packages", "Seurat is not installed in the R spatial environment", panorama_input, neighborhood_manifest, st06_scalar(ok_deconv$method_manifest_tsv), st06_scalar(ok_deconv$deconv_method), recommended_method)
} else {
  obj <- readRDS(panorama_input)
  run <- run_niche_kmeans_st(obj, ok_deconv[1, , drop = FALSE], cfg)
  if (!identical(run$status, "ok")) {
    manifest_df <- write_niche_skip_st(run$status, run$reason, panorama_input, neighborhood_manifest, st06_scalar(ok_deconv$method_manifest_tsv), st06_scalar(ok_deconv$deconv_method), recommended_method)
  } else {
    manifest_df <- data.frame(
      status = "ok",
      reason = "",
      input_rds = panorama_input,
      neighborhood_manifest = neighborhood_manifest,
      deconv_manifest = st06_scalar(ok_deconv$method_manifest_tsv),
      deconv_method = st06_scalar(ok_deconv$deconv_method),
      recommended_method = recommended_method,
      score_tsv = run$score_tsv,
      label_tsv = run$label_tsv,
      spatial_panorama_niched_rds = run$spatial_panorama_niched_rds,
      spot_n = run$spot_n,
      niche_n = run$niche_n,
      k = cfg$spatial_niche_k,
      k_neighbors = cfg$spatial_niche_k_neighbors,
      stringsAsFactors = FALSE
    )
  }
}

manifest_tsv <- file.path(cfg$spatial_niche_table_dir, "niche_manifest.tsv")
report_path <- file.path(cfg$spatial_niche_figure_dir, "niche_report.md")
st06_write_tsv(manifest_df[, niche_manifest_cols_st(), drop = FALSE], manifest_tsv)
report_lines <- c(
  "# Spatial Niche Derivation",
  "",
  sprintf("- status: `%s`", manifest_df$status[[1]]),
  sprintf("- recommended_method: `%s`", manifest_df$recommended_method[[1]]),
  sprintf("- deconv_method: `%s`", manifest_df$deconv_method[[1]]),
  "",
  "## Manifest",
  render_markdown_table_local(manifest_df[, niche_manifest_cols_st(), drop = FALSE])
)
write_markdown_local(report_lines, report_path)

st06_write_manifest_local(
  manifest_path = cfg$module_06e_niche_manifest_path,
  new_outputs = list(
    niche_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "single-row spatial niche status manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial niche derivation report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    spatial_06d_neighborhood = cfg$module_06d_neighborhood_manifest_path,
    spatial_07e_deconvolution_compare = cfg$module_07e_deconvolution_compare_manifest_path
  ),
  version = cfg$module_06_version,
  depends_on = list(
    spatial_06d_neighborhood = cfg$module_06d_neighborhood_manifest_path,
    spatial_07e_deconvolution_compare = cfg$module_07e_deconvolution_compare_manifest_path
  )
)

message(sprintf("spatial niche derivation complete: %s", manifest_df$status[[1]]))
