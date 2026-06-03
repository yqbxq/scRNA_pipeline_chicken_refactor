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
if (!nzchar(go_manifest_path) && file.exists(file.path(cfg$spatial_region_go_table_dir, "region_go_manifest.tsv"))) {
  go_manifest_path <- file.path(cfg$spatial_region_go_table_dir, "region_go_manifest.tsv")
}
kegg_manifest_path <- st06_manifest_output(cfg$module_06b_region_kegg_manifest_path, "region_kegg_manifest", cfg)
if (!nzchar(kegg_manifest_path) && file.exists(file.path(cfg$spatial_region_kegg_table_dir, "region_kegg_manifest.tsv"))) {
  kegg_manifest_path <- file.path(cfg$spatial_region_kegg_table_dir, "region_kegg_manifest.tsv")
}
go_manifest <- st06_read_tsv(go_manifest_path)
kegg_manifest <- st06_read_tsv(kegg_manifest_path)
all_manifest <- rbind(
  if (nrow(go_manifest) > 0) go_manifest else st06_empty_df(st06_enrichment_manifest_cols()),
  if (nrow(kegg_manifest) > 0) kegg_manifest else st06_empty_df(st06_enrichment_manifest_cols())
)

read_condition_values_st <- function(cfg) {
  samples <- st06_read_tsv(cfg$sample_sheet)
  condition_col <- st06_pick_col(samples, c("condition", "stage", "group_id", "group", "sample_group"))
  if (!nzchar(condition_col)) {
    return(character(0))
  }
  values <- unique(trimws(as.character(samples[[condition_col]])))
  values[nzchar(values)]
}

infer_condition_value_st <- function(item, conditions) {
  direct_col <- st06_pick_col(item, c("condition_value", "condition", "stage", "group_id", "group"))
  if (nzchar(direct_col)) {
    direct_value <- st06_scalar(item[[direct_col]])
    if (nzchar(direct_value)) {
      return(direct_value)
    }
  }
  if (length(conditions) == 0) {
    return("all")
  }
  haystack <- tolower(paste(
    st06_scalar(item$comparison_id),
    st06_scalar(item$region_value),
    st06_scalar(item$layer_id),
    sep = " "
  ))
  hits <- conditions[vapply(tolower(conditions), function(cond) grepl(cond, haystack, fixed = TRUE), logical(1))]
  if (length(hits) > 0) hits[[1]] else "unassigned"
}

cross_stage_empty_st <- function() {
  st06_empty_df(c(
    "group_by", "layer_id", "analysis_type", "ontology", "pathway_id", "pathway_label",
    "condition_a", "condition_b", "score_a", "score_b", "delta_log10_padj",
    "abs_delta_log10_padj"
  ))
}

build_cross_stage_drift_st <- function(matrix_df, top_k) {
  required <- c("group_by", "layer_id", "analysis_type", "ontology", "pathway_id", "pathway_label", "condition_value", "score")
  if (nrow(matrix_df) == 0 || any(!required %in% colnames(matrix_df))) {
    return(cross_stage_empty_st())
  }
  df <- matrix_df
  df$score_num <- suppressWarnings(as.numeric(df$score))
  df <- df[is.finite(df$score_num) & nzchar(df$condition_value) & !df$condition_value %in% c("all", "unassigned"), , drop = FALSE]
  if (nrow(df) == 0) {
    return(cross_stage_empty_st())
  }
  key_cols <- c("group_by", "layer_id", "analysis_type", "ontology", "pathway_id")
  keys <- unique(df[, key_cols, drop = FALSE])
  rows <- list()
  for (idx in seq_len(nrow(keys))) {
    key <- keys[idx, , drop = FALSE]
    hit <- df
    for (col in key_cols) {
      hit <- hit[as.character(hit[[col]]) == as.character(key[[col]][[1]]), , drop = FALSE]
    }
    if (nrow(hit) == 0) {
      next
    }
    condition_scores <- stats::aggregate(score_num ~ condition_value, data = hit, FUN = max, na.rm = TRUE)
    if (nrow(condition_scores) < 2) {
      next
    }
    condition_scores <- condition_scores[order(condition_scores$condition_value), , drop = FALSE]
    pairs <- utils::combn(seq_len(nrow(condition_scores)), 2)
    label <- st06_scalar(hit$pathway_label[1], st06_scalar(hit$pathway_id[1]))
    for (pair_idx in seq_len(ncol(pairs))) {
      a <- condition_scores[pairs[1, pair_idx], , drop = FALSE]
      b <- condition_scores[pairs[2, pair_idx], , drop = FALSE]
      delta <- b$score_num[[1]] - a$score_num[[1]]
      rows[[length(rows) + 1L]] <- data.frame(
        group_by = key$group_by,
        layer_id = key$layer_id,
        analysis_type = key$analysis_type,
        ontology = key$ontology,
        pathway_id = key$pathway_id,
        pathway_label = label,
        condition_a = a$condition_value[[1]],
        condition_b = b$condition_value[[1]],
        score_a = a$score_num[[1]],
        score_b = b$score_num[[1]],
        delta_log10_padj = delta,
        abs_delta_log10_padj = abs(delta),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(cross_stage_empty_st())
  }
  out <- do.call(rbind, rows)
  out <- out[order(-out$abs_delta_log10_padj, out$pathway_id), , drop = FALSE]
  out[seq_len(min(nrow(out), top_k)), colnames(cross_stage_empty_st()), drop = FALSE]
}

if (nrow(all_manifest) == 0) {
  summary_df <- data.frame(status = "skipped_no_upstream_manifest", reason = "06a/06b manifests are unavailable", stringsAsFactors = FALSE)
  matrix_df <- st06_empty_df(c("region_value", "comparison_id", "group_by", "layer_id", "condition_value", "pathway_id", "pathway_label", "score", "analysis_type", "ontology"))
  cross_stage_df <- cross_stage_empty_st()
  eda_status <- "skipped_no_upstream_manifest"
  eda_reason <- "06a/06b manifests are unavailable"
} else {
  condition_values <- read_condition_values_st(cfg)
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
      group_by = st06_scalar(item$group_by),
      layer_id = st06_scalar(item$layer_id),
      condition_value = infer_condition_value_st(item, condition_values),
      pathway_id = as.character(sig[[id_col]]),
      pathway_label = if (nzchar(label_col)) as.character(sig[[label_col]]) else as.character(sig[[id_col]]),
      score = score,
      analysis_type = st06_scalar(item$analysis_type),
      ontology = st06_scalar(item$ontology),
      stringsAsFactors = FALSE
    )
  }
  matrix_df <- if (length(result_rows) > 0) do.call(rbind, result_rows) else st06_empty_df(c("region_value", "comparison_id", "group_by", "layer_id", "condition_value", "pathway_id", "pathway_label", "score", "analysis_type", "ontology"))
  cross_stage_df <- build_cross_stage_drift_st(matrix_df, cfg$spatial_enrich_drift_top_k)
  eda_status <- "ok"
  eda_reason <- ""
}

summary_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_enrichment_summary.tsv")
matrix_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_pathway_matrix.tsv")
cross_stage_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "cross_stage_top_pathways.tsv")
manifest_tsv <- file.path(cfg$spatial_region_enrichment_eda_table_dir, "region_enrichment_eda_manifest.tsv")
report_path <- file.path(cfg$spatial_enrichment_report_dir, "region_enrichment_eda.md")
svg_enrichment_handoff_path <- Sys.getenv("SVG_ENRICHMENT_HANDOFF_TSV", unset = file.path(cfg$spatial_table_dir, "09_svg", "09c_consensus", "svg_enrichment_handoff.tsv"))
svg_enrichment_handoff <- st06_read_tsv(svg_enrichment_handoff_path)
st06_write_tsv(summary_df, summary_tsv)
st06_write_tsv(matrix_df, matrix_tsv)
st06_write_tsv(cross_stage_df, cross_stage_tsv)
manifest_df <- data.frame(
  status = eda_status,
  reason = eda_reason,
  go_manifest = go_manifest_path,
  kegg_manifest = kegg_manifest_path,
  summary_tsv = summary_tsv,
  matrix_tsv = matrix_tsv,
  cross_stage_top_pathways_tsv = cross_stage_tsv,
  svg_enrichment_handoff_tsv = ifelse(file.exists(svg_enrichment_handoff_path), svg_enrichment_handoff_path, ""),
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
  sprintf("- cross_stage_top_pathways: `%s`", relative_path_local(cross_stage_tsv, cfg$project_root)),
  "",
  "## GO/KEGG Status Summary",
  render_markdown_table_local(st06_status_summary(all_manifest)),
  "",
  "## Region Summary Preview",
  render_markdown_table_local(utils::head(summary_df, 50)),
  "",
  "## Cross-stage Pathway Drift",
  render_markdown_table_local(utils::head(cross_stage_df, 50)),
  "",
  "## ST09 SVG Enrichment Handoff",
  if (file.exists(svg_enrichment_handoff_path)) render_markdown_table_local(utils::head(svg_enrichment_handoff, 50)) else "No ST09 SVG enrichment handoff was available."
)
write_markdown_local(report_lines, report_path)

st06_write_manifest_local(
  manifest_path = cfg$module_06c_region_enrichment_eda_manifest_path,
  new_outputs = list(
    region_enrichment_eda_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "single-row EDA status manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    region_enrichment_summary = build_output_entry(summary_tsv, "tsv", module_name, "region-level GO/KEGG term count summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    region_pathway_matrix = build_output_entry(matrix_tsv, "tsv", module_name, "long region-pathway score matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(matrix_df)),
    cross_stage_top_pathways = build_output_entry(cross_stage_tsv, "tsv", module_name, "top pathway score drift across inferred stages", base_dir = cfg$project_root, schema = infer_schema_from_df(cross_stage_df)),
    svg_enrichment_handoff = build_output_entry(svg_enrichment_handoff_path, "tsv", module_name, "optional ST09 SVG gene set handoff", base_dir = cfg$project_root, schema = infer_schema_from_df(svg_enrichment_handoff)),
    report = build_output_entry(report_path, "md", module_name, "spatial region enrichment EDA report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_06a_region_go = cfg$module_06a_region_go_manifest_path, spatial_06b_region_kegg = cfg$module_06b_region_kegg_manifest_path, svg_enrichment_handoff = svg_enrichment_handoff_path),
  version = cfg$module_06_version,
  depends_on = list(spatial_06a_region_go = cfg$module_06a_region_go_manifest_path, spatial_06b_region_kegg = cfg$module_06b_region_kegg_manifest_path)
)

message(sprintf("spatial region enrichment EDA complete: %s", eda_status))
