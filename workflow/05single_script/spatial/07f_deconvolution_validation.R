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
module_name <- "spatial_07f_deconvolution_validation"
prepare_dirs_spatial(cfg)

read_method_props_07f <- function(row, cfg) {
  path <- st07_abs_path(row$proportion_tsv[[1]] %||% "", cfg)
  df <- st07_read_tsv(path)
  required <- c("spot_id", "cell_type", "proportion")
  if (nrow(df) == 0 || !all(required %in% colnames(df))) {
    return(st07_empty_df(c("method", "deconv_id", "section", required)))
  }
  df <- df[, required, drop = FALSE]
  df$proportion <- suppressWarnings(as.numeric(df$proportion))
  df <- df[is.finite(df$proportion), , drop = FALSE]
  df$method <- row$method[[1]]
  df$deconv_id <- row$deconv_id[[1]] %||% "default_deconv"
  df$section <- row$section[[1]] %||% "all"
  df[, c("method", "deconv_id", "section", required), drop = FALSE]
}

safe_pearson_07f <- function(a, b) {
  if (length(a) < 2 || stats::sd(a) == 0 || stats::sd(b) == 0) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(stats::cor(a, b, method = "pearson", use = "pairwise.complete.obs")))
}

normalize_truth_07f <- function(df) {
  required <- c("spot_id", "cell_type", "true_proportion")
  if (nrow(df) == 0 || !all(required %in% colnames(df))) {
    return(st07_empty_df(required))
  }
  df <- df[, required, drop = FALSE]
  df$true_proportion <- suppressWarnings(as.numeric(df$true_proportion))
  df <- df[is.finite(df$true_proportion), , drop = FALSE]
  split_rows <- split(seq_len(nrow(df)), df$spot_id)
  df$true_proportion <- ave(df$true_proportion, df$spot_id, FUN = function(x) {
    x[x < 0] <- 0
    total <- sum(x)
    if (total > 0) x / total else x
  })
  df[unlist(split_rows, use.names = FALSE), , drop = FALSE]
}

generate_dirichlet_truth_07f <- function(props, cfg) {
  spot_ids <- sort(unique(props$spot_id))
  celltypes <- sort(unique(props$cell_type))
  if (length(spot_ids) == 0 || length(celltypes) < 2) {
    return(st07_empty_df(c("spot_id", "cell_type", "true_proportion")))
  }
  set.seed(as.integer(cfg$random_seed))
  n_spots <- min(length(spot_ids), cfg$spatial_validation_n_spots)
  spot_ids <- sample(spot_ids, n_spots)
  alpha <- max(as.numeric(cfg$spatial_validation_dirichlet_alpha), 1e-6)
  mat <- matrix(stats::rgamma(n_spots * length(celltypes), shape = alpha, rate = 1), nrow = n_spots, ncol = length(celltypes))
  mat <- sweep(mat, 1, rowSums(mat), "/")
  rows <- lapply(seq_len(nrow(mat)), function(i) {
    data.frame(spot_id = spot_ids[[i]], cell_type = celltypes, true_proportion = as.numeric(mat[i, ]), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

compare_to_truth_07f <- function(props, truth) {
  methods <- sort(unique(props$method))
  rows <- lapply(methods, function(method) {
    pred <- props[props$method == method, , drop = FALSE]
    merged <- merge(truth, pred, by = c("spot_id", "cell_type"))
    if (nrow(merged) == 0) {
      return(data.frame(method = method, n_spots = 0L, n_celltypes = 0L, mean_rmse = NA_real_, mean_pearson = NA_real_, dominant_cell_type_accuracy = NA_real_, stringsAsFactors = FALSE))
    }
    celltype_rows <- lapply(sort(unique(merged$cell_type)), function(ct) {
      hit <- merged[merged$cell_type == ct, , drop = FALSE]
      data.frame(
        cell_type = ct,
        rmse = sqrt(mean((hit$proportion - hit$true_proportion)^2)),
        pearson = safe_pearson_07f(hit$proportion, hit$true_proportion),
        stringsAsFactors = FALSE
      )
    })
    celltype_df <- do.call(rbind, celltype_rows)
    pred_wide <- reshape(merged[, c("spot_id", "cell_type", "proportion")], idvar = "spot_id", timevar = "cell_type", direction = "wide")
    truth_wide <- reshape(merged[, c("spot_id", "cell_type", "true_proportion")], idvar = "spot_id", timevar = "cell_type", direction = "wide")
    common_spots <- intersect(pred_wide$spot_id, truth_wide$spot_id)
    pred_wide <- pred_wide[match(common_spots, pred_wide$spot_id), , drop = FALSE]
    truth_wide <- truth_wide[match(common_spots, truth_wide$spot_id), , drop = FALSE]
    pred_cols <- grep("^proportion\\.", colnames(pred_wide), value = TRUE)
    truth_cols <- sub("^proportion\\.", "true_proportion.", pred_cols)
    truth_cols <- truth_cols[truth_cols %in% colnames(truth_wide)]
    pred_cols <- sub("^true_proportion\\.", "proportion.", truth_cols)
    dominant_acc <- NA_real_
    if (length(pred_cols) > 0 && nrow(pred_wide) > 0) {
      pred_dom <- sub("^proportion\\.", "", pred_cols[max.col(as.matrix(pred_wide[, pred_cols, drop = FALSE]), ties.method = "first")])
      truth_dom <- sub("^true_proportion\\.", "", truth_cols[max.col(as.matrix(truth_wide[, truth_cols, drop = FALSE]), ties.method = "first")])
      dominant_acc <- mean(pred_dom == truth_dom)
    }
    data.frame(
      method = method,
      n_spots = length(unique(merged$spot_id)),
      n_celltypes = length(unique(merged$cell_type)),
      mean_rmse = mean(celltype_df$rmse, na.rm = TRUE),
      mean_pearson = mean(celltype_df$pearson, na.rm = TRUE),
      dominant_cell_type_accuracy = dominant_acc,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

methods <- collect_deconv_method_manifests_st(cfg)
ok_rows <- methods[methods$status == "ok" & nzchar(methods$method), , drop = FALSE]
prop_rows <- lapply(seq_len(nrow(ok_rows)), function(i) read_method_props_07f(ok_rows[i, , drop = FALSE], cfg))
props <- if (length(prop_rows) == 0) st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion")) else do.call(rbind, prop_rows)

truth_input <- Sys.getenv("SPATIAL_VALIDATION_TRUTH_TSV", unset = "")
truth_source <- "dirichlet_from_deconv_spots"
if (nzchar(truth_input) && file.exists(truth_input)) {
  truth <- normalize_truth_07f(st07_read_tsv(truth_input))
  truth_source <- truth_input
} else {
  truth <- generate_dirichlet_truth_07f(props, cfg)
}

if (nrow(props) == 0) {
  status <- "skipped_no_method_outputs"
  reason <- "no ok 07a/07b/07c/07d method proportion outputs were available"
  summary_df <- data.frame(method = character(), n_spots = integer(), n_celltypes = integer(), mean_rmse = numeric(), mean_pearson = numeric(), dominant_cell_type_accuracy = numeric(), stringsAsFactors = FALSE)
} else if (nrow(truth) == 0) {
  status <- "failed_synthesis"
  reason <- "synthetic truth generation produced no rows"
  summary_df <- data.frame(method = character(), n_spots = integer(), n_celltypes = integer(), mean_rmse = numeric(), mean_pearson = numeric(), dominant_cell_type_accuracy = numeric(), stringsAsFactors = FALSE)
} else {
  summary_df <- compare_to_truth_07f(props, truth)
  status <- if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_rmse))) "ok" else "failed_validation"
  reason <- if (identical(status, "ok")) "" else "no method output overlapped synthetic truth"
}

validation_id <- "deconv_validation_default"
out_dir <- file.path(cfg$spatial_deconv_validation_table_dir, validation_id)
ensure_dir(out_dir)
truth_tsv <- file.path(out_dir, "synthetic_truth.tsv")
summary_tsv <- file.path(out_dir, "method_summary.tsv")
st07_write_tsv(truth, truth_tsv)
st07_write_tsv(summary_df, summary_tsv)
manifest_tsv <- file.path(cfg$spatial_deconv_validation_table_dir, "validation_manifest.tsv")
report_path <- file.path(cfg$spatial_deconv_validation_report_dir, "report.md")
manifest_df <- data.frame(
  validation_id = validation_id,
  synthetic_n_spots = cfg$spatial_validation_n_spots,
  truth_source = truth_source,
  scdesign3_available = requireNamespace("scDesign3", quietly = TRUE),
  methods_validated = paste(summary_df$method, collapse = ","),
  mean_rmse_best_method = if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_rmse))) min(summary_df$mean_rmse, na.rm = TRUE) else NA_real_,
  status = status,
  reason = reason,
  synthetic_truth_tsv = truth_tsv,
  method_summary_tsv = summary_tsv,
  stringsAsFactors = FALSE
)
st07_write_tsv(manifest_df, manifest_tsv)

report_lines <- c(
  "# Spatial Deconvolution Validation",
  "",
  sprintf("- status: `%s`", status),
  sprintf("- synthetic_n_spots: `%s`", cfg$spatial_validation_n_spots),
  sprintf("- dirichlet_alpha: `%s`", cfg$spatial_validation_dirichlet_alpha),
  sprintf("- truth_source: `%s`", truth_source),
  sprintf("- scDesign3_available: `%s`", requireNamespace("scDesign3", quietly = TRUE)),
  "",
  "## Manifest",
  render_markdown_table_local(manifest_df),
  "",
  "## Method Summary",
  render_markdown_table_local(summary_df)
)
write_markdown_local(report_lines, report_path)

st07_write_manifest_local(
  manifest_path = cfg$module_07f_deconvolution_validation_manifest_path,
  new_outputs = list(
    validation_manifest = build_output_entry(manifest_tsv, "tsv", module_name, "deconvolution validation manifest", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    synthetic_truth = build_output_entry(truth_tsv, "tsv", module_name, "synthetic spot truth proportions", base_dir = cfg$project_root),
    method_summary = build_output_entry(summary_tsv, "tsv", module_name, "per-method validation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    report = build_output_entry(report_path, "md", module_name, "spatial deconvolution validation report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file),
  version = cfg$module_07_version,
  depends_on = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file)
)

message(sprintf("spatial deconvolution validation complete: %s", status))
