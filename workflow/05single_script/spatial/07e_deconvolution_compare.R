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

read_deconv_proportions_07e <- function(row, cfg) {
  path <- st07_abs_path(row$proportion_tsv[[1]] %||% "", cfg)
  df <- st07_read_tsv(path)
  required <- c("spot_id", "cell_type", "proportion")
  if (nrow(df) == 0 || !all(required %in% colnames(df))) {
    return(st07_empty_df(c("method", "deconv_id", "section", required)))
  }
  df$proportion <- suppressWarnings(as.numeric(df$proportion))
  df <- df[is.finite(df$proportion), required, drop = FALSE]
  df$method <- row$method[[1]]
  df$deconv_id <- row$deconv_id[[1]] %||% "default_deconv"
  df$section <- row$section[[1]] %||% "all"
  df[, c("method", "deconv_id", "section", required), drop = FALSE]
}

safe_pearson_07e <- function(a, b) {
  if (length(a) < 2 || stats::sd(a) == 0 || stats::sd(b) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(stats::cor(a, b, method = "pearson", use = "pairwise.complete.obs")))
}

safe_jsd_07e <- function(a, b) {
  a[!is.finite(a) | a < 0] <- 0
  b[!is.finite(b) | b < 0] <- 0
  if (sum(a) <= 0 || sum(b) <= 0) {
    return(NA_real_)
  }
  p <- a / sum(a)
  q <- b / sum(b)
  m <- 0.5 * (p + q)
  kl <- function(x, y) sum(ifelse(x > 0 & y > 0, x * log(x / y), 0))
  sqrt(0.5 * kl(p, m) + 0.5 * kl(q, m))
}

compare_method_pair_07e <- function(props, method_a, method_b, deconv_id, section) {
  a <- props[props$method == method_a & props$deconv_id == deconv_id & props$section == section, , drop = FALSE]
  b <- props[props$method == method_b & props$deconv_id == deconv_id & props$section == section, , drop = FALSE]
  merged <- merge(a, b, by = c("spot_id", "cell_type"), suffixes = c("_a", "_b"))
  if (nrow(merged) == 0) {
    return(st07_empty_df(c("deconv_id", "section", "method_a", "method_b", "celltype", "n_common_spots", "pearson", "rmse", "jsd")))
  }
  celltypes <- sort(unique(merged$cell_type))
  rows <- lapply(celltypes, function(ct) {
    hit <- merged[merged$cell_type == ct, , drop = FALSE]
    data.frame(
      deconv_id = deconv_id,
      section = section,
      method_a = method_a,
      method_b = method_b,
      celltype = ct,
      n_common_spots = length(unique(hit$spot_id)),
      pearson = safe_pearson_07e(hit$proportion_a, hit$proportion_b),
      rmse = sqrt(mean((hit$proportion_a - hit$proportion_b)^2)),
      jsd = safe_jsd_07e(hit$proportion_a, hit$proportion_b),
      stringsAsFactors = FALSE
    )
  })
  all_row <- data.frame(
    deconv_id = deconv_id,
    section = section,
    method_a = method_a,
    method_b = method_b,
    celltype = "all",
    n_common_spots = length(unique(merged$spot_id)),
    pearson = safe_pearson_07e(merged$proportion_a, merged$proportion_b),
    rmse = sqrt(mean((merged$proportion_a - merged$proportion_b)^2)),
    jsd = safe_jsd_07e(merged$proportion_a, merged$proportion_b),
    stringsAsFactors = FALSE
  )
  do.call(rbind, c(list(all_row), rows))
}

build_comparison_matrix_07e <- function(methods, cfg) {
  ok_rows <- methods[methods$status == "ok" & nzchar(methods$method), , drop = FALSE]
  if (nrow(ok_rows) == 0) {
    return(st07_empty_df(c("deconv_id", "section", "method_a", "method_b", "celltype", "n_common_spots", "pearson", "rmse", "jsd")))
  }
  prop_rows <- lapply(seq_len(nrow(ok_rows)), function(i) read_deconv_proportions_07e(ok_rows[i, , drop = FALSE], cfg))
  props <- do.call(rbind, prop_rows)
  if (nrow(props) == 0) {
    return(st07_empty_df(c("deconv_id", "section", "method_a", "method_b", "celltype", "n_common_spots", "pearson", "rmse", "jsd")))
  }
  groups <- unique(props[, c("deconv_id", "section"), drop = FALSE])
  rows <- list()
  for (idx in seq_len(nrow(groups))) {
    group <- groups[idx, , drop = FALSE]
    hit <- props[props$deconv_id == group$deconv_id[[1]] & props$section == group$section[[1]], , drop = FALSE]
    methods_here <- sort(unique(hit$method))
    if (length(methods_here) < 2) {
      next
    }
    for (pair in utils::combn(methods_here, 2, simplify = FALSE)) {
      rows[[length(rows) + 1L]] <- compare_method_pair_07e(props, pair[[1]], pair[[2]], group$deconv_id[[1]], group$section[[1]])
    }
  }
  if (length(rows) == 0) {
    return(st07_empty_df(c("deconv_id", "section", "method_a", "method_b", "celltype", "n_common_spots", "pearson", "rmse", "jsd")))
  }
  do.call(rbind, rows)
}

build_method_ranking_07e <- function(methods, comparison_matrix, cfg) {
  if (nrow(methods) == 0) {
    return(st07_empty_df(c("method", "ok_rows", "all_rows", "completion_score", "consensus_pearson", "mean_rmse", "mean_jsd", "primary_bonus", "recommendation_score", "recommendation_rank")))
  }
  all_method_names <- sort(unique(methods$method[nzchar(methods$method)]))
  rows <- lapply(all_method_names, function(method) {
    all_rows <- sum(methods$method == method)
    ok_rows <- sum(methods$method == method & methods$status == "ok")
    hit <- comparison_matrix[comparison_matrix$celltype == "all" & (comparison_matrix$method_a == method | comparison_matrix$method_b == method), , drop = FALSE]
    consensus_pearson <- if (nrow(hit) == 0 || all(is.na(hit$pearson))) NA_real_ else mean(hit$pearson, na.rm = TRUE)
    mean_rmse <- if (nrow(hit) == 0 || all(is.na(hit$rmse))) NA_real_ else mean(hit$rmse, na.rm = TRUE)
    mean_jsd <- if (nrow(hit) == 0 || all(is.na(hit$jsd))) NA_real_ else mean(hit$jsd, na.rm = TRUE)
    completion_score <- if (all_rows == 0) 0 else ok_rows / all_rows
    primary_bonus <- if (identical(method, cfg$spatial_deconv_primary)) 1 else 0
    consensus_score <- if (is.na(consensus_pearson)) 0 else max(0, min(1, (consensus_pearson + 1) / 2))
    data.frame(
      method = method,
      ok_rows = ok_rows,
      all_rows = all_rows,
      completion_score = completion_score,
      consensus_pearson = consensus_pearson,
      mean_rmse = mean_rmse,
      mean_jsd = mean_jsd,
      primary_bonus = primary_bonus,
      recommendation_score = 0.6 * consensus_score + 0.3 * completion_score + 0.1 * primary_bonus,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-out$recommendation_score, -out$completion_score, match(out$method, c(cfg$spatial_deconv_primary, "rctd", "cell2location", "card", "transfer")), out$method), , drop = FALSE]
  out$recommendation_rank <- seq_len(nrow(out))
  out
}

methods <- collect_deconv_method_manifests_st(cfg)
ok_methods <- unique(methods$method[methods$status == "ok" & nzchar(methods$method)])
method_counts <- as.data.frame(table(method = methods$method, status = methods$status), stringsAsFactors = FALSE)
method_counts <- method_counts[method_counts$Freq > 0, , drop = FALSE]
comparison_matrix <- build_comparison_matrix_07e(methods, cfg)

if (nrow(methods) == 0) {
  status <- "skipped_no_method_outputs"
  reason <- "no 07a/07b/07c/07d method manifests were available"
} else if (length(unique(ok_methods)) < 2 || nrow(comparison_matrix) == 0) {
  status <- "skipped_too_few_methods"
  reason <- sprintf("%s completed methods available; at least 2 required for pairwise comparison", length(ok_methods))
} else {
  status <- "ok"
  reason <- ""
}

preferred <- cfg$spatial_deconv_primary
ranking <- build_method_ranking_07e(methods, comparison_matrix, cfg)
if (nrow(ranking) > 0) {
  ok_ranking <- ranking[ranking$ok_rows > 0, , drop = FALSE]
  if (nrow(ok_ranking) > 0) {
    preferred <- ok_ranking$method[[1]]
  }
}
ensure_dir(dirname(cfg$spatial_recommended_deconv_method_file))
writeLines(preferred, cfg$spatial_recommended_deconv_method_file, useBytes = TRUE)

summary_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_ranking.tsv")
method_summary_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_summary.tsv")
matrix_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_comparison_matrix.tsv")
method_pairwise_metrics_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "method_pairwise_metrics.tsv")
celltype_consistency_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "celltype_method_consistency.tsv")
spot_consistency_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "spot_method_consistency.tsv")
recommended_by_celltype_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "recommended_method_by_celltype.tsv")
deconv_evidence_tier_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "deconv_evidence_tier.tsv")
manifest_tsv <- file.path(cfg$spatial_deconv_compare_table_dir, "deconv_compare_manifest.tsv")
report_path <- file.path(cfg$spatial_deconv_compare_report_dir, "report.md")
st07_write_tsv(ranking, summary_tsv)
st07_write_tsv(ranking, method_summary_tsv)
st07_write_tsv(comparison_matrix, matrix_tsv)
st07_write_tsv(comparison_matrix, method_pairwise_metrics_tsv)
celltype_consistency <- comparison_matrix[comparison_matrix$celltype != "all", , drop = FALSE]
spot_consistency <- st07_empty_df(c("spot_id", "method_a", "method_b", "dominant_celltype_agreement", "entropy_delta", "status", "reason"))
recommended_by_celltype <- if (nrow(celltype_consistency) == 0) {
  st07_empty_df(c("celltype", "recommended_method", "mean_pearson", "mean_rmse", "status"))
} else {
  do.call(rbind, lapply(sort(unique(celltype_consistency$celltype)), function(ct) {
    hit <- celltype_consistency[celltype_consistency$celltype == ct, , drop = FALSE]
    data.frame(
      celltype = ct,
      recommended_method = preferred,
      mean_pearson = if (all(is.na(hit$pearson))) NA_real_ else mean(hit$pearson, na.rm = TRUE),
      mean_rmse = if (all(is.na(hit$rmse))) NA_real_ else mean(hit$rmse, na.rm = TRUE),
      status = if (all(is.na(hit$pearson))) "skipped_no_metric" else "ok",
      stringsAsFactors = FALSE
    )
  }))
}
deconv_evidence_tier <- data.frame(
  method = ranking$method %||% character(),
  completion_score = ranking$completion_score %||% numeric(),
  consensus_pearson = ranking$consensus_pearson %||% numeric(),
  deconv_evidence_tier = ifelse((ranking$completion_score %||% 0) >= 1 & !is.na(ranking$consensus_pearson %||% NA_real_) & ranking$consensus_pearson >= 0.8, "primary", ifelse((ranking$completion_score %||% 0) > 0, "supporting", "blocked")),
  stringsAsFactors = FALSE
)
st07_write_tsv(celltype_consistency, celltype_consistency_tsv)
st07_write_tsv(spot_consistency, spot_consistency_tsv)
st07_write_tsv(recommended_by_celltype, recommended_by_celltype_tsv)
st07_write_tsv(deconv_evidence_tier, deconv_evidence_tier_tsv)
manifest_df <- data.frame(
  deconv_id = "all",
  methods_compared = paste(sort(ok_methods), collapse = ","),
  n_methods = length(ok_methods),
  recommended_method = preferred,
  status = status,
  reason = reason,
  method_ranking_tsv = summary_tsv,
  method_summary_tsv = method_summary_tsv,
  method_comparison_matrix_tsv = matrix_tsv,
  method_pairwise_metrics_tsv = method_pairwise_metrics_tsv,
  celltype_method_consistency_tsv = celltype_consistency_tsv,
  spot_method_consistency_tsv = spot_consistency_tsv,
  recommended_method_by_celltype_tsv = recommended_by_celltype_tsv,
  deconv_evidence_tier_tsv = deconv_evidence_tier_tsv,
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
    method_summary = build_output_entry(method_summary_tsv, "tsv", module_name, "method completion and recommendation summary", base_dir = cfg$project_root, schema = infer_schema_from_df(ranking)),
    method_comparison_matrix = build_output_entry(matrix_tsv, "tsv", module_name, "pairwise method comparison matrix", base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_matrix)),
    method_pairwise_metrics = build_output_entry(method_pairwise_metrics_tsv, "tsv", module_name, "pairwise method Pearson/RMSE/JSD metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(comparison_matrix)),
    celltype_method_consistency = build_output_entry(celltype_consistency_tsv, "tsv", module_name, "celltype-specific method consistency metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(celltype_consistency)),
    spot_method_consistency = build_output_entry(spot_consistency_tsv, "tsv", module_name, "spot-level method consistency placeholder until dominant agreement is materialized", base_dir = cfg$project_root, schema = infer_schema_from_df(spot_consistency)),
    recommended_method_by_celltype = build_output_entry(recommended_by_celltype_tsv, "tsv", module_name, "celltype-level recommended deconvolution method", base_dir = cfg$project_root, schema = infer_schema_from_df(recommended_by_celltype)),
    deconv_evidence_tier = build_output_entry(deconv_evidence_tier_tsv, "tsv", module_name, "deconvolution method evidence tier derived from completion and consensus", base_dir = cfg$project_root, schema = infer_schema_from_df(deconv_evidence_tier)),
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
