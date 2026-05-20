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
module_name <- "spatial_06c_region_enrichment_eda"
prepare_dirs_spatial(cfg)

go_manifest_path <- st06_manifest_output(cfg$module_06a_region_go_manifest_path, "region_go_manifest", cfg)
kegg_manifest_path <- st06_manifest_output(cfg$module_06b_region_kegg_manifest_path, "region_kegg_manifest", cfg)
go_manifest <- st06_read_tsv(go_manifest_path)
kegg_manifest <- st06_read_tsv(kegg_manifest_path)
all_manifest <- rbind(
  if (nrow(go_manifest) > 0) go_manifest else st06_empty_df(st06_enrichment_manifest_cols()),
  if (nrow(kegg_manifest) > 0) kegg_manifest else st06_empty_df(st06_enrichment_manifest_cols())
)

if (nrow(all_manifest) == 0) {
  summary_df <- data.frame(status = "skipped_no_upstream_manifest", reason = "06a/06b manifests are unavailable", stringsAsFactors = FALSE)
  matrix_df <- st06_empty_df(c("region_value", "comparison_id", "layer_id", "pathway_id", "pathway_label", "score"))
  eda_status <- "skipped_no_upstream_manifest"
  eda_reason <- "06a/06b manifests are unavailable"
} else {
  all_manifest$significant_term_n_num <- suppressWarnings(as.numeric(all_manifest$significant_term_n))
  all_manifest$significant_term_n_num[!is.finite(all_manifest$significant_term_n_num)] <- 0
  keys <- unique(all_manifest[, intersect(c("region_value", "comparison_id", "layer_id", "analysis_type", "ontology", "status"), colnames(all_manifest)), drop = FALSE])
  summary_rows <- lapply(seq_len(nrow(keys)), function(idx) {
    key <- keys[idx, , drop = FALSE]
    hit <- all_manifest
    for (col in colnames(key)) {
      hit <- hit[as.character(hit[[col]]) == as.character(key[[col]][[1]]), , drop = FALSE]
    }
    data.frame(key, manifest_row_n = nrow(hit), significant_term_n = sum(hit$significant_term_n_num, na.rm = TRUE), stringsAsFactors = FALSE)
  })
  summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else st06_empty_df(c("status", "manifest_row_n", "significant_term_n"))

  result_rows <- list()
  for (idx in seq_len(nrow(all_manifest))) {
    item <- all_manifest[idx, , drop = FALSE]
    sig_path <- st06_abs_path(item$enrichment_significant_tsv, cfg)
    sig <- st06_read_tsv(sig_path)
    if (nrow(sig) == 0) {
      next
    }
    id_col <- st06_pick_col(sig, c("ID", "id", "term_id", "native"))
    label_col <- st06_pick_col(sig, c("Description", "description", "term_name", "name"))
    score_col <- st06_pick_col(sig, c("qvalue", "p.adjust", "p_value", "pvalue"))
    if (!nzchar(id_col)) {
      next
    }
    score <- if (nzchar(score_col)) -log10(pmax(suppressWarnings(as.numeric(sig[[score_col]])), 1e-300)) else rep(NA_real_, nrow(sig))
    result_rows[[length(result_rows) + 1L]] <- data.frame(
      region_value = st06_scalar(item$region_value),
      comparison_id = st06_scalar(item$comparison_id),
      layer_id = st06_scalar(item$layer_id),
      pathway_id = as.character(sig[[id_col]]),
      pathway_label = if (nzchar(label_col)) as.character(sig[[label_col]]) else as.character(sig[[id_col]]),
      score = score,
      analysis_type = st06_scalar(item$analysis_type),
      ontology = st06_scalar(item$ontology),
      stringsAsFactors = FALSE
    )
  }
  matrix_df <- if (length(result_rows) > 0) do.call(rbind, result_rows) else st06_empty_df(c("region_value", "comparison_id", "layer_id", "pathway_id", "pathway_label", "score", "analysis_type", "ontology"))
  eda_status <- "ok"
  eda_reason <- ""
}

summary_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_enrichment_summary.tsv")
matrix_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_pathway_matrix.tsv")
manifest_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_enrichment_eda_manifest.tsv")
report_path <- file.path(cfg$spatial_enrichment_report_dir, "region_enrichment_eda.md")
st06_write_tsv(summary_df, summary_tsv)
st06_write_tsv(matrix_df, matrix_tsv)
manifest_df <- data.frame(
  status = eda_status,
  reason = eda_reason,
  go_manifest = go_manifest_path,
  kegg_manifest = kegg_manifest_path,
  summary_tsv = summary_tsv,
  matrix_tsv = matrix_tsv,
  report = report_path,
  stringsAsFactors = FALSE
)
st06_write_tsv(manifest_df, manifest_tsv)

report_lines <- c(
  "# Spatial Region Enrichment EDA",
  "",
  sprintf("- status: `%s`", eda_status),
  sprintf("- summary: `%s`", relative_path_local(summary_tsv, cfg$project_root)),
  sprintf("- matrix: `%s`", relative_path_local(matrix_tsv, cfg$project_root)),
  "",
  "## GO/KEGG Status Summary",
  render_markdown_table_local(st06_status_summary(all_manifest)),
  "",
  "## Region Summary Preview",
  render_markdown_table_local(utils::head(summary_df, 50))
)
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_06c_region_enrichment_eda_manifest_path,
  new_outputs = list(
    region_enrichment_eda_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "single-row EDA status manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    region_enrichment_summary = build_output_entry(summary_tsv, "tsv", module_name, "region-level GO/KEGG term count summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    region_pathway_matrix = build_output_entry(matrix_tsv, "tsv", module_name, "long region-pathway score matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(matrix_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial region enrichment EDA report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_06a_region_go = cfg$module_06a_region_go_manifest_path, spatial_06b_region_kegg = cfg$module_06b_region_kegg_manifest_path),
  version = cfg$module_06_version,
  depends_on = list(spatial_06a_region_go = cfg$module_06a_region_go_manifest_path, spatial_06b_region_kegg = cfg$module_06b_region_kegg_manifest_path)
)

message(sprintf("spatial region enrichment EDA complete: %s", eda_status))
