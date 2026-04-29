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
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite", "Matrix"))

cfg <- get_single_script_config_05()
module_name <- "05b_pseudobulk_de"
prepare_dirs_05(cfg)
set.seed(cfg$random_seed)

require_formal_pkg_05 <- function(pkgs, label) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("%s missing required R packages: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
}

run_formal_pbds_05 <- function(seu, vars) {
  require_formal_pkg_05(c("SingleCellExperiment", "S4Vectors", "limma", "edgeR", "muscat"), "formal pseudobulk DE")

  meta_df <- standardize_design_metadata_df_05(seu@meta.data)
  meta_df$group_id <- trimws(as.character(meta_df[[vars$group_var]]))
  meta_df$sample_id <- trimws(as.character(meta_df$sample_id))
  meta_df$cluster_id <- trimws(as.character(seu$cluster_id))
  meta_df$batch <- if (vars$batch_var %in% colnames(meta_df)) trimws(as.character(meta_df[[vars$batch_var]])) else "default"
  meta_df$batch[!nzchar(meta_df$batch)] <- "default"

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = get_assay_matrix(seu, assay = "RNA", type = "counts")),
    colData = S4Vectors::DataFrame(meta_df)
  )
  sce <- muscat::prepSCE(sce, kid = "cluster_id", sid = "sample_id", gid = "group_id", drop = FALSE)
  pb <- muscat::aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

  sample_info <- unique(meta_df[, c("sample_id", "group_id", "batch"), drop = FALSE])
  sample_info <- sample_info[match(colnames(pb), sample_info$sample_id), , drop = FALSE]
  sample_info$group_id <- factor(sample_info$group_id, levels = c(vars$ident_1, vars$ident_2))
  rownames(sample_info) <- sample_info$sample_id

  design_formula <- if (length(unique(sample_info$batch)) >= 2) "~ 0 + group_id + batch" else "~ 0 + group_id"
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_info)
  rownames(design) <- rownames(sample_info)
  contrast <- limma::makeContrasts(
    contrasts = sprintf("group_id%s-group_id%s", make.names(vars$ident_1), make.names(vars$ident_2)),
    levels = design
  )

  res <- muscat::pbDS(
    pb = pb,
    method = "limma-trend",
    design = design,
    contrast = contrast,
    min_cells = 10,
    filter = "both",
    verbose = FALSE
  )
  tidy_res <- muscat::resDS(sce, res, bind = "row")
  tidy_res$comparison_id <- vars$comparison_id
  tidy_res
}

layer_status_df <- deg_layer_status(cfg)
comparison_df <- read_deg_comparison_sheet(cfg)
marker_manifest <- file.path(cfg$marker_discovery_table_dir, "marker_discovery_manifest.tsv")

report_lines <- c(
  "# 05b Pseudobulk DE",
  "",
  sprintf("- comparison_sheet: `%s`", cfg$comparison_sheet),
  sprintf("- min_biological_replicates_default: `%s`", cfg$min_biological_replicates),
  ""
)
manifest_rows <- list()

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("05b pseudobulk layer: ", layer_id)
  obj <- load_layer_for_deg(cfg, layer_row)

  layer_paths <- pseudobulk_paths_05(cfg, layer_id)
  ensure_dir(layer_paths$table_dir)
  aggregation <- aggregate_cluster_sample_counts_05(obj, cluster_var = "cluster_id")
  saveRDS(aggregation$counts, layer_paths$aggregation_rds)
  write_tsv_local(aggregation$metadata, layer_paths$aggregation_tsv)
  write_tsv_local(
    data.frame(
      layer_id = layer_id,
      aggregation_rds = normalizePath(layer_paths$aggregation_rds, winslash = "/", mustWork = FALSE),
      aggregation_metadata_tsv = normalizePath(layer_paths$aggregation_tsv, winslash = "/", mustWork = FALSE),
      pseudobulk_columns = ncol(aggregation$counts),
      stringsAsFactors = FALSE
    ),
    layer_paths$aggregation_manifest_tsv
  )

  report_lines <- c(
    report_lines,
    sprintf("## Layer `%s`", layer_id),
    sprintf("- aggregation_counts: `%s`", layer_paths$aggregation_rds),
    sprintf("- aggregation_metadata: `%s`", layer_paths$aggregation_tsv)
  )

  for (cmp_idx in seq_len(nrow(comparison_df))) {
    comparison_row <- comparison_df[cmp_idx, , drop = FALSE]
    if (!comparison_applies_to_layer(comparison_row, layer_id)) {
      next
    }

    vars <- resolve_comparison_vars_05(obj, comparison_row)
    cmp_paths <- pseudobulk_paths_05(cfg, layer_id, vars$comparison_id)
    subset_result <- subset_cells_for_comparison(obj, vars)
    marker_fallback <- marker_discovery_paths_05(cfg, layer_id, vars$comparison_id)$exploratory_tsv

    if (!identical(subset_result$status, "ok")) {
      write_tsv_local(empty_gate_summary_05(), cmp_paths$gate_summary_tsv)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          inference_status = "skipped",
          reason = subset_result$reason,
          warning_banner = "",
          exploratory_marker_discovery = marker_fallback,
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = vars$comparison_id,
        inference_status = "skipped",
        aggregation_rds = normalizePath(layer_paths$aggregation_rds, winslash = "/", mustWork = FALSE),
        ds_results_tsv = "",
        gate_summary_tsv = normalizePath(cmp_paths$gate_summary_tsv, winslash = "/", mustWork = FALSE),
        status_tsv = normalizePath(cmp_paths$status_tsv, winslash = "/", mustWork = FALSE),
        exploratory_marker_discovery = marker_fallback,
        stringsAsFactors = FALSE
      )
      report_lines <- c(report_lines, sprintf("### `%s`", vars$comparison_id), sprintf("- status: skipped; reason: %s", subset_result$reason))
      next
    }

    obj_sub <- subset_result$object
    gate <- replicate_gate_summary_05(obj_sub@meta.data, vars)
    write_tsv_local(gate$summary, cmp_paths$gate_summary_tsv)

    if (!gate$pass) {
      status <- if (vars$force_exploratory) "exploratory_forced" else "exploratory_only"
      warning_banner <- annotation_warning_banner_05(gate$summary, vars, forced = vars$force_exploratory)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          inference_status = status,
          reason = gate$reason,
          warning_banner = warning_banner,
          exploratory_marker_discovery = marker_fallback,
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = vars$comparison_id,
        inference_status = status,
        aggregation_rds = normalizePath(layer_paths$aggregation_rds, winslash = "/", mustWork = FALSE),
        ds_results_tsv = "",
        gate_summary_tsv = normalizePath(cmp_paths$gate_summary_tsv, winslash = "/", mustWork = FALSE),
        status_tsv = normalizePath(cmp_paths$status_tsv, winslash = "/", mustWork = FALSE),
        exploratory_marker_discovery = marker_fallback,
        stringsAsFactors = FALSE
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", vars$comparison_id),
        warning_banner,
        sprintf("- status: %s", status),
        sprintf("- reason: %s", gate$reason),
        sprintf("- exploratory_fallback: `%s`", marker_fallback)
      )
      next
    }

    formal_res <- run_formal_pbds_05(obj_sub, vars)
    write_tsv_local(formal_res, cmp_paths$ds_results_tsv)
    saveRDS(formal_res, cmp_paths$result_rds)
    write_tsv_local(
      data.frame(
        comparison_id = vars$comparison_id,
        layer_id = layer_id,
        inference_status = "formal",
        reason = gate$reason,
        warning_banner = "",
        exploratory_marker_discovery = marker_fallback,
        stringsAsFactors = FALSE
      ),
      cmp_paths$status_tsv
    )
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      comparison_id = vars$comparison_id,
      inference_status = "formal",
      aggregation_rds = normalizePath(layer_paths$aggregation_rds, winslash = "/", mustWork = FALSE),
      ds_results_tsv = normalizePath(cmp_paths$ds_results_tsv, winslash = "/", mustWork = FALSE),
      gate_summary_tsv = normalizePath(cmp_paths$gate_summary_tsv, winslash = "/", mustWork = FALSE),
      status_tsv = normalizePath(cmp_paths$status_tsv, winslash = "/", mustWork = FALSE),
      exploratory_marker_discovery = marker_fallback,
      stringsAsFactors = FALSE
    )
    report_lines <- c(report_lines, sprintf("### `%s`", vars$comparison_id), "- status: formal", sprintf("- formal_results: `%s`", cmp_paths$ds_results_tsv))
  }
}

manifest_df <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "inference_status", "aggregation_rds",
  "ds_results_tsv", "gate_summary_tsv", "status_tsv", "exploratory_marker_discovery"
))
manifest_tsv <- file.path(cfg$pseudobulk_table_dir, "pseudobulk_manifest.tsv")
write_tsv_local(manifest_df, manifest_tsv)
report_lines <- c(report_lines, "", sprintf("- manifest: `%s`", manifest_tsv))
report_path <- file.path(cfg$pseudobulk_report_dir, "report.md")
write_markdown_local(report_lines, report_path)

if (file.exists(cfg$module_05b_manifest_path)) {
  unlink(cfg$module_05b_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_05b_manifest_path,
  new_outputs = list(
    pseudobulk_manifest_tsv = build_output_entry(manifest_tsv, "tsv", module_name, "one row per layer/comparison pseudobulk DE status", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "pseudobulk DE report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(layer_status_tsv = cfg$layer_status_file, comparison_sheet = cfg$comparison_sheet, marker_manifest_tsv = marker_manifest),
  version = cfg$module_version,
  depends_on = list(module_05a = cfg$module_05a_manifest_path)
)

message("05b completed. manifest: ", manifest_tsv)
