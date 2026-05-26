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

load_required_packages(c("Seurat", "Matrix", "jsonlite"))

cfg <- get_spatial_script_config()
module_name <- "spatial_02a_integration_eda"
prepare_dirs_spatial(cfg)
set.seed(cfg$random_seed)

parse_cli <- function(args) {
  out <- list(
    selected_method_tsv = cfg$selected_method_tsv,
    layer_file = cfg$spatial_object_layer_file,
    marker_panel_dir = file.path(cfg$config_dir, "marker_panels"),
    modes = cfg$integration_run_modes
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--selected-method-tsv") && i < length(args)) {
      out$selected_method_tsv <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--layer-file") && i < length(args)) {
      out$layer_file <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--marker-panel-dir") && i < length(args)) {
      out$marker_panel_dir <- args[[i + 1L]]
      i <- i + 2L
    } else if (identical(arg, "--modes") && i < length(args)) {
      out$modes <- args[[i + 1L]]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
layer_settings <- spatial_layer_settings(args$layer_file, "panorama_st")
hvg_n <- as.integer(spatial_numeric_or(layer_settings$hvg_nfeatures %||% layer_settings$hvg_n, cfg$default_hvg_n))
dims <- spatial_parse_dims(layer_settings$pca_dims %||% "1:30", fallback = 1:30)
npcs <- max(dims)
modes <- spatial_tokenize(args$modes)
if (length(modes) == 0 || any(tolower(modes) == "all")) {
  modes <- c("none", "harmony", "cca")
}
modes <- unique(tolower(modes))
if (!"none" %in% modes) {
  modes <- c("none", modes)
}

obj_list <- load_post_norm_objects(cfg, args$selected_method_tsv)
panorama <- merge_spatial_panorama(obj_list)
hvg <- compute_panorama_hvg(panorama, method = "selected", hvg_n = hvg_n)

panel <- tryCatch(read_spatial_region_panel(args$marker_panel_dir, layer_id = "panorama_st"), error = function(e) NULL)
results <- list()
bio_tables <- list()

started <- proc.time()[["elapsed"]]
none_result <- tryCatch(run_integration_none(panorama, hvg, npcs = npcs), error = function(e) e)
if (inherits(none_result, "error")) {
  stop(sprintf("required none integration failed: %s", conditionMessage(none_result)), call. = FALSE)
}
panorama <- none_result
lisi_none <- compute_lisi(panorama, "pca_none", group_col = "section_id")
bio_none <- tryCatch(compute_biological_consistency(panorama, panel, "pca_none"), error = function(e) data.frame())
bio_tables[["none"]] <- bio_none
results[["none"]] <- list(
  status = "ok",
  reduction = "pca_none",
  umap = "umap_none",
  mean_lisi = lisi_none$mean_lisi,
  lisi_implementation = lisi_none$implementation,
  biological_consistency = if (nrow(bio_none) == 0) NA_real_ else mean(bio_none$marker_aligned_cluster_frac, na.rm = TRUE),
  message = lisi_none$message,
  timing_sec = round(proc.time()[["elapsed"]] - started, 3)
)

if ("harmony" %in% modes) {
  started <- proc.time()[["elapsed"]]
  harmony_obj <- tryCatch(run_integration_harmony(panorama, hvg, npcs = npcs, group_by = "section_id"), error = function(e) e)
  if (inherits(harmony_obj, "error")) {
    results[["harmony"]] <- list(status = if (!requireNamespace("harmony", quietly = TRUE)) "failed_dependency" else "failed", reduction = "harmony", umap = "umap_harmony", message = conditionMessage(harmony_obj), timing_sec = round(proc.time()[["elapsed"]] - started, 3))
  } else {
    panorama <- harmony_obj
    lisi_harmony <- compute_lisi(panorama, "harmony", group_col = "section_id")
    bio_harmony <- tryCatch(compute_biological_consistency(panorama, panel, "harmony"), error = function(e) data.frame())
    bio_tables[["harmony"]] <- bio_harmony
    results[["harmony"]] <- list(status = "ok", reduction = "harmony", umap = "umap_harmony", mean_lisi = lisi_harmony$mean_lisi, lisi_implementation = lisi_harmony$implementation, biological_consistency = if (nrow(bio_harmony) == 0) NA_real_ else mean(bio_harmony$marker_aligned_cluster_frac, na.rm = TRUE), message = lisi_harmony$message, timing_sec = round(proc.time()[["elapsed"]] - started, 3))
  }
}

if ("cca" %in% modes) {
  started <- proc.time()[["elapsed"]]
  cca_obj <- tryCatch(run_integration_cca(panorama, hvg, npcs = npcs), error = function(e) e)
  if (inherits(cca_obj, "error")) {
    results[["cca"]] <- list(status = "failed_dependency", reduction = "cca_integrated", umap = "umap_cca", message = conditionMessage(cca_obj), timing_sec = round(proc.time()[["elapsed"]] - started, 3))
  } else {
    panorama <- cca_obj
    lisi_cca <- compute_lisi(panorama, "cca_integrated", group_col = "section_id")
    bio_cca <- tryCatch(compute_biological_consistency(panorama, panel, "cca_integrated"), error = function(e) data.frame())
    bio_tables[["cca"]] <- bio_cca
    results[["cca"]] <- list(status = "ok", reduction = "cca_integrated", umap = "umap_cca", mean_lisi = lisi_cca$mean_lisi, lisi_implementation = lisi_cca$implementation, biological_consistency = if (nrow(bio_cca) == 0) NA_real_ else mean(bio_cca$marker_aligned_cluster_frac, na.rm = TRUE), message = lisi_cca$message, timing_sec = round(proc.time()[["elapsed"]] - started, 3))
  }
}

panorama@misc$integration <- modifyList(
  panorama@misc$integration %||% list(),
  list(
    run_modes = modes,
    recommended_mode = cfg$default_integration_mode,
    selected_mode = select_integration_mode(args$layer_file, default = cfg$default_integration_mode),
    selected_method_tsv = args$selected_method_tsv
  )
)
saveRDS(panorama, cfg$spatial_panorama_integrated_rds)

lisi_summary <- summarize_integration(results)
bio_df <- if (length(bio_tables) == 0) data.frame(
  mode = character(),
  region = character(),
  n_top_spots = integer(),
  marker_aligned_cluster_frac = numeric(),
  integration_mode = character(),
  stringsAsFactors = FALSE
) else do.call(rbind, lapply(names(bio_tables), function(mode) {
  df <- bio_tables[[mode]]
  if (ncol(df) == 0) {
    df <- data.frame(mode = character(), region = character(), n_top_spots = integer(), marker_aligned_cluster_frac = numeric(), stringsAsFactors = FALSE)
  }
  df$integration_mode <- mode
  df
}))
recommended <- data.frame(
  integration_mode_default = cfg$default_integration_mode,
  recommended_mode = cfg$default_integration_mode,
  final_choice = select_integration_mode(args$layer_file, default = cfg$default_integration_mode),
  decided_by = if (nzchar(Sys.getenv("SPATIAL_INTEGRATION_MODE", unset = "")) && !identical(Sys.getenv("SPATIAL_INTEGRATION_MODE", unset = ""), cfg$default_integration_mode)) "env_override" else "spatial_object_layers",
  decision_note = "Default remains none because stage and section may be confounded at low N; review diagnostics before changing spatial_object_layers.tsv.",
  stringsAsFactors = FALSE
)

spatial_write_tsv(lisi_summary, file.path(cfg$spatial_integration_compare_dir, "lisi_summary.tsv"))
spatial_write_tsv(bio_df, file.path(cfg$spatial_integration_compare_dir, "biological_consistency.tsv"))
spatial_write_tsv(recommended, file.path(cfg$spatial_integration_compare_dir, "recommended_mode.tsv"))

report_lines <- c(
  "# Spatial Integration EDA",
  "",
  sprintf("- Panorama object: `%s`", relative_path_local(cfg$spatial_panorama_integrated_rds, cfg$project_root)),
  sprintf("- Modes requested: `%s`", paste(modes, collapse = ", ")),
  sprintf("- Recommended mode: `%s`", cfg$default_integration_mode),
  "",
  "UMAP embeddings are computed by `spatial_02a2_compute_umap` after this reduction/integration EDA stage.",
  "",
  "## Mode Summary",
  render_markdown_table_local(lisi_summary),
  "",
  "## Biological Consistency",
  render_markdown_table_local(bio_df),
  "",
  "## Recommendation",
  render_markdown_table_local(recommended)
)
write_markdown_local(report_lines, file.path(cfg$spatial_integration_compare_dir, "report.md"))

write_manifest_local(
  manifest_path = cfg$module_02a_integration_manifest_path,
  new_outputs = list(
    panorama_integrated = build_output_entry(cfg$spatial_panorama_integrated_rds, "rds", module_name, "merged spatial panorama with integration reductions", base_dir = cfg$project_root),
    lisi_summary = build_output_entry(file.path(cfg$spatial_integration_compare_dir, "lisi_summary.tsv"), "tsv", module_name, "one row per integration mode LISI summary", base_dir = cfg$project_root, schema = infer_schema_from_df(lisi_summary)),
    biological_consistency = build_output_entry(file.path(cfg$spatial_integration_compare_dir, "biological_consistency.tsv"), "tsv", module_name, "marker-module biological consistency by mode", base_dir = cfg$project_root, schema = infer_schema_from_df(bio_df)),
    recommended_mode = build_output_entry(file.path(cfg$spatial_integration_compare_dir, "recommended_mode.tsv"), "tsv", module_name, "selected integration mode decision", base_dir = cfg$project_root, schema = infer_schema_from_df(recommended)),
    report = build_output_entry(file.path(cfg$spatial_integration_compare_dir, "report.md"), "md", module_name, "spatial integration EDA report", base_dir = cfg$project_root),
    figures_dir = build_output_entry(cfg$spatial_integration_figure_dir, "directory", module_name, "integration diagnostic figures", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(selected_method_tsv = args$selected_method_tsv, spatial_object_layer_file = args$layer_file),
  version = cfg$module_03_version,
  depends_on = list(spatial_02_normalize = cfg$module_02_norm_manifest_path)
)

message(sprintf("spatial integration EDA complete: %d sections, %d modes", length(obj_list), length(results)))
