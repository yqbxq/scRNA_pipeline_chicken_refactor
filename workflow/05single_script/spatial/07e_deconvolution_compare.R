#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_common.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "project_paths_spatial.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "spatial_deconv_utils.R"), encoding = "UTF-8")

cfg <- get_spatial_script_config()
module_name <- "spatial_07e_deconvolution_compare"
prepare_dirs_spatial(cfg)

methods <- collect_deconv_method_manifests_st(cfg)
ok_methods <- unique(methods$method[methods$status == "ok" & nzchar(methods$method)])
method_counts <- as.data.frame(table(method = methods$method, status = methods$status), stringsAsFactors = FALSE)
method_counts <- method_counts[method_counts$Freq > 0, , drop = FALSE]

if (nrow(methods) == 0) {
  status <- "skipped_no_method_outputs"
  reason <- "no 07a/07b/07c/07d method manifests were available"
} else if (length(ok_methods) < 2) {
  status <- "skipped_too_few_methods"
  reason <- sprintf("%s completed methods available; at least 2 required for pairwise comparison", length(ok_methods))
} else {
  status <- "ok"
  reason <- ""
}

preferred <- cfg$spatial_deconv_primary
if (length(ok_methods) > 0) {
  preferred <- if (cfg$spatial_deconv_primary %in% ok_methods) cfg$spatial_deconv_primary else ok_methods[[1]]
}
ensure_dir(dirname(cfg$spatial_recommended_deconv_method_file))
writeLines(preferred, cfg$spatial_recommended_deconv_method_file, useBytes = TRUE)

ranking <- if (nrow(methods) == 0) {
  st07_empty_df(c("method", "ok_rows", "all_rows", "completion_score", "recommendation_rank"))
} else {
  all_method_names <- sort(unique(methods$method[nzchar(methods$method)]))
  rows <- lapply(all_method_names, function(method) {
    all_rows <- sum(methods$method == method)
    ok_rows <- sum(methods$method == method & methods$status == "ok")
    data.frame(method = method, ok_rows = ok_rows, all_rows = all_rows, completion_score = if (all_rows == 0) 0 else ok_rows / all_rows, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$completion_score, match(out$method, c(cfg$spatial_deconv_primary, "transfer", "card", "cell2location")), out$method), , drop = FALSE]
  out$recommendation_rank <- seq_len(nrow(out))
  out
}

comparison_matrix <- if (length(ok_methods) < 2) {
  st07_empty_df(c("method_a", "method_b", "celltype", "pearson", "rmse", "jsd"))
} else {
  pairs <- utils::combn(sort(ok_methods), 2, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    data.frame(method_a = pair[[1]], method_b = pair[[2]], celltype = "all", pearson = NA_real_, rmse = NA_real_, jsd = NA_real_, stringsAsFactors = FALSE)
  }))
}

summary_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_ranking.tsv")
matrix_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_comparison_matrix.tsv")
manifest_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "deconv_compare_manifest.tsv")
report_path <- file.path(cfg$spatial_deconv_compare_report_dir, "report.md")
st07_write_tsv(ranking, summary_tsv)
st07_write_tsv(comparison_matrix, matrix_tsv)
manifest_df <- data.frame(
  deconv_id = "all",
  methods_compared = paste(sort(ok_methods), collapse = ","),
  n_methods = length(ok_methods),
  recommended_method = preferred,
  status = status,
  reason = reason,
  method_ranking_tsv = summary_tsv,
  method_comparison_matrix_tsv = matrix_tsv,
  recommended_method_file = cfg$spatial_recommended_deconv_method_file,
  stringsAsFactors = FALSE
)
st07_write_tsv(manifest_df, manifest_tsv)

report_lines <- c(
  "# Spatial Deconvolution Compare",
  "",
  sprintf("- status: `%s`", status),
  sprintf("- recommended_method: `%s`", preferred),
  sprintf("- recommended_method_file: `%s`", relative_path_local(cfg$spatial_recommended_deconv_method_file, cfg$project_root)),
  "",
  "## Method Status Counts",
  render_markdown_table_local(method_counts),
  "",
  "## Method Ranking",
  render_markdown_table_local(ranking)
)
write_markdown_local(report_lines, report_path)

st07_write_manifest_local(
  manifest_path = cfg$module_07e_deconvolution_compare_manifest_path,
  new_outputs = list(
    deconv_compare_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "single-row deconvolution compare status manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    method_ranking = build_output_entry(summary_tsv, "tsv", module_name, "method recommendation ranking", base_dir = cfg$project_root, schema = infer_schema_from_df(ranking)),
    method_comparison_matrix = build_output_entry(matrix_tsv, "tsv", module_name, "pairwise method comparison matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_matrix)),
    recommended_method = build_output_entry(cfg$spatial_recommended_deconv_method_file, "txt", module_name, "recommended deconvolution method for 06e niche derivation", base_dir = cfg$project_root),
    report = build_output_entry(report_path, "md", module_name, "spatial deconvolution comparison report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    spatial_07a_deconvolution_rctd = cfg$module_07a_deconvolution_rctd_manifest_path,
    spatial_07b_deconvolution_transfer = cfg$module_07b_deconvolution_transfer_manifest_path,
    spatial_07c_deconvolution_card = cfg$module_07c_deconvolution_card_manifest_path,
    spatial_07d_deconvolution_cell2location = cfg$module_07d_deconvolution_cell2location_manifest_path
  ),
  version = cfg$module_07_version,
  depends_on = list(spatial_07a_deconvolution_rctd = cfg$module_07a_deconvolution_rctd_manifest_path)
)

message(sprintf("spatial deconvolution compare complete: %s", status))
