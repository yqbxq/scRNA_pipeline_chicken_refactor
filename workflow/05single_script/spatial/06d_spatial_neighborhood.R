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
module_name <- "spatial_06d_neighborhood"
prepare_dirs_spatial(cfg)

empty_neighborhood_manifest <- function() {
  st06_empty_df(c(
    "status", "reason", "input_rds", "spot_tsv", "label_column", "spot_n",
    "label_n", "radius", "permutations", "coord_type", "interaction_matrix_tsv",
    "nhood_enrichment_zscore_tsv", "co_occurrence_tsv", "sidecar_summary_tsv"
  ))
}

resolve_panorama_input <- function(cfg) {
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

write_skip_manifest <- function(status, reason, input_rds = "") {
  spot_tsv <- file.path(cfg$spatial_neighborhood_table_dir, "neighborhood_spots.tsv")
  st06_write_tsv(st06_empty_df(c("spot_id", "x", "y", "label", "section_id")), spot_tsv)
  data.frame(
    status = status,
    reason = reason,
    input_rds = input_rds,
    spot_tsv = spot_tsv,
    label_column = "",
    spot_n = 0L,
    label_n = 0L,
    radius = cfg$spatial_neighborhood_radius,
    permutations = cfg$spatial_neighborhood_perms,
    coord_type = cfg$spatial_neighborhood_coord_type,
    interaction_matrix_tsv = "",
    nhood_enrichment_zscore_tsv = "",
    co_occurrence_tsv = "",
    sidecar_summary_tsv = "",
    stringsAsFactors = FALSE
  )
}

panorama_input <- resolve_panorama_input(cfg)
if (!nzchar(panorama_input)) {
  manifest_df <- write_skip_manifest("skipped_no_upstream_manifest", "no annotated/subannotated spatial panorama RDS was available")
} else if (!requireNamespace("Seurat", quietly = TRUE)) {
  manifest_df <- write_skip_manifest("skipped_no_packages", "Seurat is not installed in the R spatial environment", panorama_input)
} else {
  obj <- readRDS(panorama_input)
  meta <- obj@meta.data
  label_col <- c("sub_region", "region", "seurat_clusters")[c("sub_region", "region", "seurat_clusters") %in% colnames(meta)][1]
  x_col <- c("x", "X", "imagecol", "pxl_col_in_fullres", "array_col", "col")[c("x", "X", "imagecol", "pxl_col_in_fullres", "array_col", "col") %in% colnames(meta)][1]
  y_col <- c("y", "Y", "imagerow", "pxl_row_in_fullres", "array_row", "row")[c("y", "Y", "imagerow", "pxl_row_in_fullres", "array_row", "row") %in% colnames(meta)][1]
  if (is.na(label_col) || is.na(x_col) || is.na(y_col)) {
    manifest_df <- write_skip_manifest("failed", "metadata must contain label column region/sub_region and x/y coordinate columns", panorama_input)
  } else {
    spots <- data.frame(
      spot_id = rownames(meta),
      x = suppressWarnings(as.numeric(meta[[x_col]])),
      y = suppressWarnings(as.numeric(meta[[y_col]])),
      label = as.character(meta[[label_col]]),
      section_id = if ("section_id" %in% colnames(meta)) as.character(meta$section_id) else "",
      stringsAsFactors = FALSE
    )
    spots <- spots[is.finite(spots$x) & is.finite(spots$y) & nzchar(spots$label), , drop = FALSE]
    spot_tsv <- file.path(cfg$spatial_neighborhood_table_dir, "neighborhood_spots.tsv")
    st06_write_tsv(spots, spot_tsv)
    if (nrow(spots) < 3 || length(unique(spots$label)) < 2) {
      manifest_df <- write_skip_manifest("skipped_too_few_spots", "need at least 3 spots and 2 region labels", panorama_input)
      manifest_df$spot_tsv <- spot_tsv
      manifest_df$spot_n <- nrow(spots)
      manifest_df$label_n <- length(unique(spots$label))
      manifest_df$label_column <- label_col
    } else {
      out_prefix <- file.path(cfg$spatial_neighborhood_table_dir, "neighborhood")
      cmd <- c(
        cfg$spatial_neighborhood_py,
        "--spots", spot_tsv,
        "--out-prefix", out_prefix,
        "--radius", as.character(cfg$spatial_neighborhood_radius),
        "--perms", as.character(cfg$spatial_neighborhood_perms),
        "--seed", as.character(cfg$random_seed)
      )
      run <- tryCatch(
        system2(cfg$py_spatial_bin, args = cmd, stdout = TRUE, stderr = TRUE),
        error = function(e) structure(conditionMessage(e), status = 127)
      )
      code <- attr(run, "status")
      if (is.null(code)) {
        code <- 0L
      }
      reason <- paste(as.character(run), collapse = " ")
      status <- if (identical(as.integer(code), 0L)) "ok" else if (grepl("skipped_no_packages", reason)) "skipped_no_packages" else if (grepl("skipped_no_python_env", reason)) "skipped_no_python_env" else if (grepl("skipped_too_few_spots", reason)) "skipped_too_few_spots" else "failed_sidecar"
      manifest_df <- data.frame(
        status = status,
        reason = if (identical(status, "ok")) "" else reason,
        input_rds = panorama_input,
        spot_tsv = spot_tsv,
        label_column = label_col,
        spot_n = nrow(spots),
        label_n = length(unique(spots$label)),
        radius = cfg$spatial_neighborhood_radius,
        permutations = cfg$spatial_neighborhood_perms,
        coord_type = cfg$spatial_neighborhood_coord_type,
        interaction_matrix_tsv = paste0(out_prefix, "_interaction_matrix.tsv"),
        nhood_enrichment_zscore_tsv = paste0(out_prefix, "_nhood_enrichment_zscore.tsv"),
        co_occurrence_tsv = paste0(out_prefix, "_co_occurrence.tsv"),
        sidecar_summary_tsv = paste0(out_prefix, "_summary.tsv"),
        stringsAsFactors = FALSE
      )
    }
  }
}

manifest_tsv <- file.path(cfg$spatial_neighborhood_table_dir, "neighborhood_manifest.tsv")
report_path <- file.path(cfg$spatial_neighborhood_figure_dir, "neighborhood_report.md")
st06_write_tsv(manifest_df, manifest_tsv)
report_lines <- c(
  "# Spatial Neighborhood",
  "",
  sprintf("- status: `%s`", manifest_df$status[[1]]),
  sprintf("- label_column: `%s`", manifest_df$label_column[[1]]),
  sprintf("- radius: `%s`", cfg$spatial_neighborhood_radius),
  sprintf("- permutations: `%s`", cfg$spatial_neighborhood_perms),
  "",
  "## Manifest",
  render_markdown_table_local(manifest_df)
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_06d_neighborhood_manifest_path,
  new_outputs = list(
    neighborhood_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "single-row spatial neighborhood status manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial neighborhood report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(panorama = panorama_input, sidecar = cfg$spatial_neighborhood_py),
  version = cfg$module_06_version,
  depends_on = list(spatial_06c_region_enrichment_eda = cfg$module_06c_region_enrichment_eda_manifest_path)
)

message(sprintf("spatial neighborhood complete: %s", manifest_df$status[[1]]))
