#!/usr/bin/env Rscript

.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
})

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/helpers/report_utils.R"), encoding = "UTF-8")
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

truth_wide_07f <- function(truth) {
  if (nrow(truth) == 0) {
    return(st07_empty_df(c("spot_id")))
  }
  reshape(truth, idvar = "spot_id", timevar = "cell_type", direction = "wide")
}

prediction_wide_07f <- function(props) {
  if (nrow(props) == 0) {
    return(st07_empty_df(c("method", "spot_id")))
  }
  reshape(props[, c("method", "spot_id", "cell_type", "proportion"), drop = FALSE], idvar = c("method", "spot_id"), timevar = "cell_type", direction = "wide")
}

celltype_metrics_07f <- function(props, truth) {
  merged <- merge(truth, props, by = c("spot_id", "cell_type"))
  if (nrow(merged) == 0) {
    return(data.frame(method = character(), cell_type = character(), n_spots = integer(), rmse = numeric(), pearson = numeric(), stringsAsFactors = FALSE))
  }
  rows <- lapply(split(merged, list(merged$method, merged$cell_type), drop = TRUE), function(hit) {
    data.frame(
      method = hit$method[[1]],
      cell_type = hit$cell_type[[1]],
      n_spots = length(unique(hit$spot_id)),
      rmse = sqrt(mean((hit$proportion - hit$true_proportion)^2)),
      pearson = safe_pearson_07f(hit$proportion, hit$true_proportion),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

spot_metrics_07f <- function(props, truth) {
  merged <- merge(truth, props, by = c("spot_id", "cell_type"))
  if (nrow(merged) == 0) {
    return(data.frame(method = character(), spot_id = character(), rmse = numeric(), pearson = numeric(), dominant_match = logical(), stringsAsFactors = FALSE))
  }
  rows <- lapply(split(merged, list(merged$method, merged$spot_id), drop = TRUE), function(hit) {
    data.frame(
      method = hit$method[[1]],
      spot_id = hit$spot_id[[1]],
      rmse = sqrt(mean((hit$proportion - hit$true_proportion)^2)),
      pearson = safe_pearson_07f(hit$proportion, hit$true_proportion),
      dominant_match = hit$cell_type[which.max(hit$proportion)] == hit$cell_type[which.max(hit$true_proportion)],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

gate_status_from_validation_07f <- function(status, validation_mode, summary_df, cfg) {
  if (validation_mode == "dirichlet_only") {
    return(list(gate_status = "WARN", interpretation_allowed = "exploratory", reason = "dirichlet_only validation is smoke/fallback evidence and cannot unlock I11 PASS."))
  }
  if (!identical(status, "ok")) {
    return(list(gate_status = "FAIL", interpretation_allowed = "no", reason = status))
  }
  ok_methods <- nrow(summary_df[is.finite(summary_df$mean_rmse), , drop = FALSE])
  best_rmse <- if (ok_methods > 0) min(summary_df$mean_rmse, na.rm = TRUE) else Inf
  best_cor <- if (ok_methods > 0) max(summary_df$mean_pearson, na.rm = TRUE) else -Inf
  if (ok_methods >= 2 && best_rmse <= cfg$spatial_validation_rmse_pass && best_cor >= cfg$spatial_validation_cor_pass) {
    return(list(gate_status = "PASS", interpretation_allowed = "yes", reason = "synthetic deconvolution benchmark passed thresholds."))
  }
  if (ok_methods >= 1 && best_rmse <= cfg$spatial_validation_rmse_warn && best_cor >= cfg$spatial_validation_cor_warn) {
    return(list(gate_status = "WARN", interpretation_allowed = "exploratory", reason = "synthetic deconvolution benchmark met warning thresholds."))
  }
  list(gate_status = "FAIL", interpretation_allowed = "no", reason = "synthetic deconvolution benchmark did not meet thresholds.")
}

methods <- collect_deconv_method_manifests_st(cfg)
ok_rows <- methods[methods$status == "ok" & nzchar(methods$method), , drop = FALSE]
prop_rows <- lapply(seq_len(nrow(ok_rows)), function(i) read_method_props_07f(ok_rows[i, , drop = FALSE], cfg))
props <- if (length(prop_rows) == 0) st07_empty_df(c("method", "deconv_id", "section", "spot_id", "cell_type", "proportion")) else do.call(rbind, prop_rows)

truth_input <- Sys.getenv("SPATIAL_VALIDATION_TRUTH_TSV", unset = "")
validation_mode <- tolower(cfg$spatial_validation_mode)
if (nzchar(truth_input) && file.exists(truth_input)) {
  validation_mode <- "external_truth"
}
truth_source <- if (identical(validation_mode, "external_truth")) truth_input else "dirichlet_from_deconv_spots"
if (identical(validation_mode, "scdesign3") && !requireNamespace("scDesign3", quietly = TRUE)) {
  truth <- st07_empty_df(c("spot_id", "cell_type", "true_proportion"))
  truth_source <- "scdesign3_unavailable"
} else if (identical(validation_mode, "external_truth") && nzchar(truth_input) && file.exists(truth_input)) {
  truth <- normalize_truth_07f(st07_read_tsv(truth_input))
} else {
  validation_mode <- "dirichlet_only"
  truth <- generate_dirichlet_truth_07f(props, cfg)
}

if (identical(cfg$spatial_validation_mode, "scdesign3") && !requireNamespace("scDesign3", quietly = TRUE)) {
  status <- "skipped_no_packages"
  reason <- "SPATIAL_VALIDATION_MODE=scdesign3 requested but scDesign3 is not installed in this R environment"
  summary_df <- data.frame(method = character(), n_spots = integer(), n_celltypes = integer(), mean_rmse = numeric(), mean_pearson = numeric(), dominant_cell_type_accuracy = numeric(), stringsAsFactors = FALSE)
} else if (nrow(props) == 0) {
  status <- "skipped_no_method_outputs"
  reason <- "no ok 07a/07b/07c/07d method proportion outputs were available"
  summary_df <- data.frame(method = character(), n_spots = integer(), n_celltypes = integer(), mean_rmse = numeric(), mean_pearson = numeric(), dominant_cell_type_accuracy = numeric(), stringsAsFactors = FALSE)
} else if (nrow(truth) == 0) {
  status <- "failed_synthesis"
  reason <- "synthetic truth generation produced no rows"
  summary_df <- data.frame(method = character(), n_spots = integer(), n_celltypes = integer(), mean_rmse = numeric(), mean_pearson = numeric(), dominant_cell_type_accuracy = numeric(), stringsAsFactors = FALSE)
} else {
  summary_df <- compare_to_truth_07f(props, truth)
  status <- if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_rmse))) if (identical(validation_mode, "dirichlet_only")) "ok_smoke" else "ok" else "failed_validation"
  reason <- if (status %in% c("ok", "ok_smoke")) "" else "no method output overlapped synthetic truth"
}

validation_id <- "deconv_validation_default"
out_dir <- file.path(cfg$spatial_deconv_validation_table_dir, validation_id)
ensure_dir(out_dir)
truth_tsv <- file.path(out_dir, "synthetic_truth.tsv")
truth_wide_tsv <- file.path(out_dir, "synthetic_truth_wide.tsv")
generation_summary_tsv <- file.path(out_dir, "synthetic_generation_summary.tsv")
prediction_long_tsv <- file.path(out_dir, "method_prediction_long.tsv")
prediction_wide_tsv <- file.path(out_dir, "method_prediction_wide.tsv")
summary_tsv <- file.path(out_dir, "method_summary.tsv")
celltype_metrics_tsv <- file.path(out_dir, "celltype_metrics.tsv")
spot_metrics_tsv <- file.path(out_dir, "spot_metrics.tsv")
recommended_consistency_tsv <- file.path(out_dir, "recommended_method_consistency.tsv")
question_gate_tsv <- file.path(cfg$spatial_deconv_validation_table_dir, "spatial_question_gate_status.tsv")
st07_write_tsv(truth, truth_tsv)
st07_write_tsv(truth_wide_07f(truth), truth_wide_tsv)
generation_summary <- data.frame(
  validation_mode = validation_mode,
  truth_source = truth_source,
  synthetic_cell_generation = if (identical(validation_mode, "scdesign3") && requireNamespace("scDesign3", quietly = TRUE)) "pending_runtime" else "not_used",
  synthetic_spot_generation = if (nrow(truth) > 0) "ok_truth_table" else "not_generated",
  truth_spot_n = length(unique(truth$spot_id %||% character())),
  truth_celltype_n = length(unique(truth$cell_type %||% character())),
  stringsAsFactors = FALSE
)
st07_write_tsv(generation_summary, generation_summary_tsv)
st07_write_tsv(props, prediction_long_tsv)
st07_write_tsv(prediction_wide_07f(props), prediction_wide_tsv)
st07_write_tsv(summary_df, summary_tsv)
celltype_metrics <- celltype_metrics_07f(props, truth)
spot_metrics <- spot_metrics_07f(props, truth)
st07_write_tsv(celltype_metrics, celltype_metrics_tsv)
st07_write_tsv(spot_metrics, spot_metrics_tsv)
recommended_consistency <- data.frame(
  recommended_method = if (nrow(summary_df) > 0) summary_df$method[order(summary_df$mean_rmse, -summary_df$mean_pearson)][[1]] else "",
  validation_status = status,
  validation_mode = validation_mode,
  best_rmse = if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_rmse))) min(summary_df$mean_rmse, na.rm = TRUE) else NA_real_,
  best_pearson = if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_pearson))) max(summary_df$mean_pearson, na.rm = TRUE) else NA_real_,
  stringsAsFactors = FALSE
)
st07_write_tsv(recommended_consistency, recommended_consistency_tsv)
gate <- gate_status_from_validation_07f(status, validation_mode, summary_df, cfg)
question_gates <- do.call(rbind, list(
  data.frame(question_id = "I08_deconv_panorama", module = module_name, gate_status = gate$gate_status, interpretation_allowed = gate$interpretation_allowed, reason = gate$reason, stringsAsFactors = FALSE),
  data.frame(question_id = "I09_deconv_GC_subtype", module = module_name, gate_status = gate$gate_status, interpretation_allowed = gate$interpretation_allowed, reason = gate$reason, stringsAsFactors = FALSE),
  data.frame(question_id = "I10_deconv_TC_subtype", module = module_name, gate_status = "PLANNED", interpretation_allowed = "no", reason = "TC subtype deconvolution remains planned until TC subtype labels are stable.", stringsAsFactors = FALSE),
  data.frame(question_id = "I11_deconv_validation", module = module_name, gate_status = gate$gate_status, interpretation_allowed = gate$interpretation_allowed, reason = gate$reason, stringsAsFactors = FALSE),
  data.frame(question_id = "I12_ST_region_compo", module = module_name, gate_status = if (gate$gate_status %in% c("PASS", "WARN")) "WARN" else "FAIL", interpretation_allowed = if (gate$gate_status %in% c("PASS", "WARN")) "exploratory" else "no", reason = "ST region composition depends on 07e recommendation and 07f validation status.", stringsAsFactors = FALSE)
))
st07_write_tsv(question_gates, question_gate_tsv)
manifest_tsv <- file.path(cfg$spatial_deconv_validation_table_dir, "validation_manifest.tsv")
report_path <- file.path(cfg$spatial_deconv_validation_report_dir, "report.md")
manifest_df <- data.frame(
  validation_id = validation_id,
  validation_mode = validation_mode,
  synthetic_n_spots = cfg$spatial_validation_n_spots,
  truth_source = truth_source,
  scdesign3_available = requireNamespace("scDesign3", quietly = TRUE),
  methods_validated = paste(summary_df$method, collapse = ","),
  mean_rmse_best_method = if (nrow(summary_df) > 0 && any(is.finite(summary_df$mean_rmse))) min(summary_df$mean_rmse, na.rm = TRUE) else NA_real_,
  status = status,
  reason = reason,
  synthetic_truth_tsv = truth_tsv,
  synthetic_truth_wide_tsv = truth_wide_tsv,
  synthetic_generation_summary_tsv = generation_summary_tsv,
  method_summary_tsv = summary_tsv,
  celltype_metrics_tsv = celltype_metrics_tsv,
  spot_metrics_tsv = spot_metrics_tsv,
  spatial_question_gate_status_tsv = question_gate_tsv,
  stringsAsFactors = FALSE
)
st07_write_tsv(manifest_df, manifest_tsv)

report_lines <- c(
  "# Spatial Deconvolution Validation",
  "",
  sprintf("- status: `%s`", status),
  sprintf("- validation_mode: `%s`", validation_mode),
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
    synthetic_generation_summary = build_output_entry(generation_summary_tsv, "tsv", module_name, "synthetic generation status and provenance", base_dir = cfg$project_root, schema = infer_schema_from_df(generation_summary)),
    synthetic_truth_wide = build_output_entry(truth_wide_tsv, "tsv", module_name, "wide synthetic spot truth proportions", base_dir = cfg$project_root),
    method_prediction_long = build_output_entry(prediction_long_tsv, "tsv", module_name, "method prediction long table used for validation", base_dir = cfg$project_root),
    method_prediction_wide = build_output_entry(prediction_wide_tsv, "tsv", module_name, "method prediction wide table used for validation", base_dir = cfg$project_root),
    method_summary = build_output_entry(summary_tsv, "tsv", module_name, "per-method validation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    celltype_metrics = build_output_entry(celltype_metrics_tsv, "tsv", module_name, "per-method per-celltype validation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(celltype_metrics)),
    spot_metrics = build_output_entry(spot_metrics_tsv, "tsv", module_name, "per-method per-spot validation metrics", base_dir = cfg$project_root, schema = infer_schema_from_df(spot_metrics)),
    recommended_method_consistency = build_output_entry(recommended_consistency_tsv, "tsv", module_name, "validation-derived method recommendation consistency", base_dir = cfg$project_root, schema = infer_schema_from_df(recommended_consistency)),
    spatial_question_gate_status = build_output_entry(question_gate_tsv, "tsv", module_name, "I08/I09/I10/I11/I12 deconvolution question gates", base_dir = cfg$project_root, schema = infer_schema_from_df(question_gates)),
    report = build_output_entry(report_path, "md", module_name, "spatial deconvolution validation report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file),
  version = cfg$module_07_version,
  depends_on = list(spatial_reference_inventory = cfg$spatial_reference_inventory_file)
)

message(sprintf("spatial deconvolution validation complete: %s", status))
