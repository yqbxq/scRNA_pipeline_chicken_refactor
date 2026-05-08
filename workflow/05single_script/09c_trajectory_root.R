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
source_utf8(file.path(.script_dir, "helpers", "project_paths_06.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_07.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_09.R"))
source_utf8(file.path(.script_dir, "helpers", "trajectory_utils.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite", "Matrix", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09c_trajectory_root"
prepare_dirs_09(cfg)

root_inference_cols <- c(
  "pair_id", "split_value", "method", "recommended_root", "score",
  "agree_with_prior", "status", "reason"
)
root_index_cols <- c(
  "pair_id", "split_value", "selected_root", "selected_source",
  "prior_root", "cytotrace_root", "velocity_root", "n_methods_ok",
  "status", "reason", "root_inference_tsv"
)
split_agreement_cols <- c("pair_id", "split_values", "selected_roots", "status", "reason")

first_or_empty_09c <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) "" else x[[1]]
}

cluster_labels_09c <- function(seu, label_var) {
  labels <- as.character(seu@meta.data[[label_var]])
  labels[is.na(labels) | !nzchar(labels)] <- "unknown"
  labels
}

prior_root_row_09c <- function(pair_id, split_value, prior_root, labels) {
  prior_root <- normalize_scalar_value(prior_root, "auto")
  if (prior_root %in% c("auto", "*")) {
    return(data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "prior",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "skipped_auto",
      reason = "root_group is auto or wildcard",
      stringsAsFactors = FALSE
    ))
  }
  present <- prior_root %in% unique(labels)
  data.frame(
    pair_id = pair_id,
    split_value = split_value,
    method = "prior",
    recommended_root = if (present) prior_root else "",
    score = if (present) 1 else 0,
    agree_with_prior = if (present) "yes" else "no",
    status = if (present) "ok" else "missing_prior_label",
    reason = if (present) "" else sprintf("root_group=%s absent from label_var", prior_root),
    stringsAsFactors = FALSE
  )
}

cytotrace_scores_09c <- function(seu) {
  counts <- trajectory_counts_matrix_09(seu)
  if (requireNamespace("CytoTRACE", quietly = TRUE)) {
    result <- tryCatch(
      CytoTRACE::CytoTRACE(as.matrix(counts)),
      error = function(e) e
    )
    if (!inherits(result, "error")) {
      score <- result$CytoTRACE %||% result$cytotrace %||% result$CytoTRACE_Score
      if (!is.null(score)) {
        score <- as.numeric(score)
        names(score) <- colnames(counts)
        return(list(score = score, status = "ok", reason = "CytoTRACE package"))
      }
    }
  }

  nfeature <- tryCatch(Matrix::colSums(counts > 0), error = function(e) rep(NA_real_, ncol(counts)))
  ncount <- tryCatch(Matrix::colSums(counts), error = function(e) rep(NA_real_, ncol(counts)))
  proxy <- trajectory_scale01_09(log1p(as.numeric(nfeature))) + trajectory_scale01_09(log1p(as.numeric(ncount)))
  names(proxy) <- colnames(counts)
  list(score = proxy, status = "proxy", reason = "CytoTRACE package unavailable; used nFeature/nCount proxy")
}

cytotrace_root_row_09c <- function(seu, pair_id, split_value, label_var, prior_root) {
  labels <- cluster_labels_09c(seu, label_var)
  score <- cytotrace_scores_09c(seu)
  df <- data.frame(label = labels, score = as.numeric(score$score[colnames(seu)]), stringsAsFactors = FALSE)
  df <- df[is.finite(df$score), , drop = FALSE]
  if (nrow(df) == 0) {
    return(data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "cytotrace",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "skipped_no_scores",
      reason = score$reason,
      stringsAsFactors = FALSE
    ))
  }
  med <- stats::aggregate(score ~ label, data = df, FUN = median)
  med <- med[order(med$score, decreasing = TRUE), , drop = FALSE]
  root <- med$label[[1]]
  data.frame(
    pair_id = pair_id,
    split_value = split_value,
    method = if (identical(score$status, "ok")) "cytotrace" else "cytotrace_proxy",
    recommended_root = root,
    score = med$score[[1]],
    agree_with_prior = if (normalize_scalar_value(prior_root) %in% c("", "auto", "*")) "not_applicable" else if (identical(root, prior_root)) "yes" else "no",
    status = score$status,
    reason = score$reason,
    stringsAsFactors = FALSE
  )
}

velocity_root_row_09c <- function(cfg, pair_id, split_value, prior_root) {
  unit_id <- trajectory_unit_file_id_09(pair_id, split_value)
  candidates <- c(
    file.path(cfg$velocity_output_dir, sprintf("velocity_root_terminal_%s.tsv", unit_id)),
    file.path(cfg$velocity_output_dir, sprintf("velocity_root_terminal_%s__%s.tsv", safe_id_09(pair_id), safe_id_09(split_value))),
    file.path(cfg$velocity_output_dir, sprintf("velocity_root_terminal_%s.tsv", safe_id_09(pair_id))),
    Sys.glob(file.path(cfg$velocity_output_dir, sprintf("*%s*root*terminal*.tsv", unit_id)))
  )
  candidates <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (length(candidates) == 0) {
    return(data.frame(
      pair_id = pair_id,
      split_value = display_scalar_value(split_value, "pooled"),
      method = "velocity",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "skipped_no_velocity",
      reason = "10h velocity root/terminal output not found",
      stringsAsFactors = FALSE
    ))
  }
  vt <- read_tsv_optional(candidates[[1]])
  if (nrow(vt) == 0) {
    return(data.frame(
      pair_id = pair_id,
      split_value = display_scalar_value(split_value, "pooled"),
      method = "velocity",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "skipped_empty_velocity",
      reason = candidates[[1]],
      stringsAsFactors = FALSE
    ))
  }
  cluster_col <- intersect(c("cluster", "label", "cell_type", "cell_subtype", "state"), colnames(vt))
  score_col <- intersect(c("initial_score", "root_score", "initial_probability", "initial_states"), colnames(vt))
  if (length(cluster_col) == 0 || length(score_col) == 0) {
    return(data.frame(
      pair_id = pair_id,
      split_value = display_scalar_value(split_value, "pooled"),
      method = "velocity",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "skipped_bad_velocity_schema",
      reason = sprintf("missing cluster/score columns in %s", candidates[[1]]),
      stringsAsFactors = FALSE
    ))
  }
  vt[[score_col[[1]]]] <- suppressWarnings(as.numeric(vt[[score_col[[1]]]]))
  vt <- vt[order(vt[[score_col[[1]]]], decreasing = TRUE), , drop = FALSE]
  root <- as.character(vt[[cluster_col[[1]]]][[1]])
  data.frame(
    pair_id = pair_id,
    split_value = display_scalar_value(split_value, "pooled"),
    method = "velocity",
    recommended_root = root,
    score = vt[[score_col[[1]]]][[1]],
    agree_with_prior = if (normalize_scalar_value(prior_root) %in% c("", "auto", "*")) "not_applicable" else if (identical(root, prior_root)) "yes" else "no",
    status = "ok",
    reason = candidates[[1]],
    stringsAsFactors = FALSE
  )
}

select_root_09c <- function(inference, prior_root) {
  ok <- inference[inference$status %in% c("ok", "proxy") & nzchar(inference$recommended_root), , drop = FALSE]
  if (nrow(ok) == 0) {
    return(list(root = "", source = "", status = "failed_no_root", reason = "no root method returned a label"))
  }
  prior_root <- normalize_scalar_value(prior_root, "auto")
  if (!prior_root %in% c("auto", "*")) {
    prior_hit <- ok[ok$method == "prior" & ok$recommended_root == prior_root, , drop = FALSE]
    if (nrow(prior_hit) > 0) {
      return(list(root = prior_root, source = "prior", status = "ok", reason = "prior label present"))
    }
  }
  tab <- sort(table(ok$recommended_root), decreasing = TRUE)
  root <- names(tab)[[1]]
  source <- paste(ok$method[ok$recommended_root == root], collapse = ",")
  list(root = root, source = source, status = "ok", reason = "method vote")
}

units <- trajectory_execution_units_09(cfg)
root_rows <- list()
index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  input_rds <- unit$input_rds[[1]]
  out_path <- trajectory_root_inference_path_09(cfg, pair_id, if (identical(split_value, "pooled")) "" else split_value)
  ensure_dir(dirname(out_path))

  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    seu <- readRDS(input_rds)
    label_var <- normalize_scalar_value(unit$coarse_label_var[[1]])
    labels <- cluster_labels_09c(seu, label_var)
    prior_root <- normalize_scalar_value(unit$root_group[[1]], "auto")
    inference <- dplyr::bind_rows(
      prior_root_row_09c(pair_id, split_value, prior_root, labels),
      cytotrace_root_row_09c(seu, pair_id, split_value, label_var, prior_root),
      velocity_root_row_09c(cfg, pair_id, split_value, prior_root)
    )
    selected <- select_root_09c(inference, prior_root)
    write_tsv_local(inference, out_path)
    list(inference = inference, selected = selected, n_cells = ncol(seu), status = selected$status, reason = selected$reason)
  }, error = function(e) {
    inference <- data.frame(
      pair_id = pair_id,
      split_value = split_value,
      method = "root_inference",
      recommended_root = "",
      score = NA_real_,
      agree_with_prior = "not_applicable",
      status = "failed",
      reason = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    write_tsv_local(inference, out_path)
    list(inference = inference, selected = list(root = "", source = "", status = "failed", reason = conditionMessage(e)), n_cells = 0L, status = "failed", reason = conditionMessage(e))
  })

  root_rows[[length(root_rows) + 1L]] <- result$inference
  index_rows[[length(index_rows) + 1L]] <- data.frame(
    pair_id = pair_id,
    split_value = split_value,
    selected_root = result$selected$root,
    selected_source = result$selected$source,
    prior_root = normalize_scalar_value(unit$root_group[[1]], "auto"),
    cytotrace_root = first_or_empty_09c(result$inference$recommended_root[grepl("^cytotrace", result$inference$method)]),
    velocity_root = first_or_empty_09c(result$inference$recommended_root[result$inference$method == "velocity"]),
    n_methods_ok = sum(result$inference$status %in% c("ok", "proxy")),
    status = result$status,
    reason = result$reason,
    root_inference_tsv = out_path,
    runtime_s = round(proc.time()[["elapsed"]] - started, 3),
    stringsAsFactors = FALSE
  )
  dynamic_outputs[[sprintf("root_inference_%s", trajectory_unit_file_id_09(pair_id, split_value))]] <- build_output_entry(
    out_path,
    "tsv",
    module_name,
    "root inference votes for one trajectory pair/split",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(result$inference)
  )
}

root_index <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(root_index_cols)
write_tsv_local(root_index, cfg$trajectory_root_index_tsv)

split_rows <- list()
if (nrow(root_index) > 0) {
  for (pair_id in unique(root_index$pair_id)) {
    sub <- root_index[root_index$pair_id == pair_id, , drop = FALSE]
    split_values <- unique(sub$split_value)
    roots <- unique(sub$selected_root[nzchar(sub$selected_root)])
    if (length(split_values) > 1) {
      split_rows[[length(split_rows) + 1L]] <- data.frame(
        pair_id = pair_id,
        split_values = paste(split_values, collapse = ","),
        selected_roots = paste(sprintf("%s=%s", sub$split_value, sub$selected_root), collapse = ","),
        status = if (length(roots) <= 1) "ok" else "split_mode_root_disagreement",
        reason = if (length(roots) <= 1) "" else "selected roots differ across split values",
        stringsAsFactors = FALSE
      )
    }
  }
}
split_agreement <- if (length(split_rows) > 0) dplyr::bind_rows(split_rows) else trajectory_empty_df_09(split_agreement_cols)
write_tsv_local(split_agreement, cfg$trajectory_root_split_agreement_tsv)

outputs <- c(
  list(
    trajectory_root_index = build_output_entry(cfg$trajectory_root_index_tsv, "tsv", module_name, "selected root by trajectory pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(root_index)),
    split_root_agreement = build_output_entry(cfg$trajectory_root_split_agreement_tsv, "tsv", module_name, "cross-split selected root agreement", base_dir = cfg$project_root, schema = infer_schema_from_df(split_agreement))
  ),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09c_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, trajectory_pairs = cfg$trajectory_pairs_sheet, velocity_output_dir = cfg$velocity_output_dir),
  depends_on = list(module_09a = cfg$module_09a_manifest_path)
)

message("09c completed. root index: ", cfg$trajectory_root_index_tsv)
