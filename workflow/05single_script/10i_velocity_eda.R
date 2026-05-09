#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_10.R"), encoding = "UTF-8")

load_required_packages(c("dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_10()
module_name <- "10i_velocity_eda"
prepare_dirs_10(cfg)

empty_module_status_10i <- function() {
  data.frame(
    pair_id = character(0),
    value = character(0),
    method = character(0),
    methods_enabled = character(0),
    status = character(0),
    n_cells_used = integer(0),
    runtime_s = numeric(0),
    metric_value = character(0),
    notes = character(0),
    stringsAsFactors = FALSE
  )
}

read_manifest_output_or_fallback_10i <- function(manifest_path, key, fallback_path) {
  path <- velocity_manifest_output_optional_10(manifest_path, key)
  if (!nzchar(path)) fallback_path else path
}

read_tsv_with_cols_10i <- function(path, cols) {
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    return(velocity_empty_df_10(cols))
  }
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df
}

read_csv_optional_10i <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

safe_numeric_10i <- function(x) {
  suppressWarnings(as.numeric(x))
}

method_index_10i <- function(manifest_path, key, fallback_path, cols) {
  path <- read_manifest_output_or_fallback_10i(manifest_path, key, fallback_path)
  read_tsv_with_cols_10i(path, cols)
}

status_from_method_index_10i <- function(index_df, fallback_method, metric_fun = NULL) {
  if (nrow(index_df) == 0) {
    return(empty_module_status_10i())
  }
  rows <- lapply(seq_len(nrow(index_df)), function(i) {
    row <- index_df[i, , drop = FALSE]
    notes <- normalize_scalar_value(row$reason[[1]])
    if (!nzchar(notes) && "output_path" %in% colnames(row)) {
      notes <- normalize_scalar_value(row$output_path[[1]])
    }
    metric_value <- ""
    if (!is.null(metric_fun)) {
      metric_value <- normalize_scalar_value(metric_fun(row))
    }
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = normalize_scalar_value(row$method[[1]], fallback_method),
      methods_enabled = display_scalar_value(row$methods_enabled[[1]], "yes"),
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = suppressWarnings(as.integer(row$n_cells[[1]])),
      runtime_s = suppressWarnings(as.numeric(row$runtime_s[[1]])),
      metric_value = metric_value,
      notes = notes,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_loom_10i <- function(loom_index) {
  if (nrow(loom_index) == 0) {
    return(empty_module_status_10i())
  }
  rows <- lapply(seq_len(nrow(loom_index)), function(i) {
    row <- loom_index[i, , drop = FALSE]
    data.frame(
      pair_id = normalize_scalar_value(row$sample_id[[1]], "__PROJECT__"),
      value = "sample",
      method = "velocyto_loom",
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = NA_integer_,
      runtime_s = suppressWarnings(as.numeric(row$runtime_s[[1]])),
      metric_value = "",
      notes = normalize_scalar_value(row$reason[[1]], normalize_scalar_value(row$loom_path[[1]])),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_reference_10i <- function(reference_index) {
  if (nrow(reference_index) == 0) {
    return(empty_module_status_10i())
  }
  rows <- lapply(seq_len(nrow(reference_index)), function(i) {
    row <- reference_index[i, , drop = FALSE]
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = "velocity_reference",
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = suppressWarnings(as.integer(row$n_cells[[1]])),
      runtime_s = NA_real_,
      metric_value = "",
      notes = normalize_scalar_value(row$reason[[1]], normalize_scalar_value(row$metadata_csv[[1]])),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

metric_value_from_driver_10i <- function(row) {
  sprintf("driver_n=%s; velocity_gene_n=%s", display_scalar_value(row$driver_n[[1]], "NA"), display_scalar_value(row$velocity_gene_n[[1]], "NA"))
}

metric_value_from_consistency_index_10i <- function(row) {
  normalize_scalar_value(row$split_compare_path[[1]])
}

metric_value_from_root_terminal_10i <- function(row) {
  normalize_scalar_value(row$root_terminal_tsv[[1]])
}

metric_value_from_cellrank_10i <- function(row) {
  paste(
    c(
      sprintf("fate=%s", display_scalar_value(row$fate_csv[[1]], "NA")),
      sprintf("macrostates=%s", display_scalar_value(row$macrostates_tsv[[1]], "NA"))
    ),
    collapse = "; "
  )
}

metric_value_from_qc_10i <- function(qc_df) {
  if (nrow(qc_df) == 0) {
    return(velocity_empty_df_10(c("pair_id", "split_value", "n_cells", "unspliced_spliced_ratio", "velocity_confidence_mean", "velocity_length_mean", "status", "reason")))
  }
  cols <- c("pair_id", "split_value", "n_cells", "unspliced_spliced_ratio", "velocity_confidence_mean", "velocity_length_mean", "status", "reason")
  for (col in cols) {
    if (!col %in% colnames(qc_df)) qc_df[[col]] <- ""
  }
  qc_df[, cols, drop = FALSE]
}

summarize_cellrank_fate_10i <- function(cellrank_index) {
  rows <- list()
  ok <- cellrank_index[cellrank_index$status == "ok" & nzchar(cellrank_index$macrostates_tsv), , drop = FALSE]
  for (i in seq_len(nrow(ok))) {
    macro <- read_tsv_optional(ok$macrostates_tsv[[i]])
    if (nrow(macro) == 0) next
    state_col <- if ("max_fate_state" %in% colnames(macro)) "max_fate_state" else ""
    prob_col <- if ("max_fate_probability" %in% colnames(macro)) "max_fate_probability" else ""
    if (!nzchar(state_col)) next
    states <- as.character(macro[[state_col]])
    tab <- sort(table(states[nzchar(states)]), decreasing = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      pair_id = ok$pair_id[[i]],
      split_value = display_scalar_value(ok$split_value[[i]], "pooled"),
      cell_n = nrow(macro),
      fate_state_n = length(tab),
      top_fate_state = if (length(tab) > 0) names(tab)[[1]] else "",
      top_fate_cell_n = if (length(tab) > 0) unname(tab[[1]]) else 0L,
      max_fate_probability_median = if (nzchar(prob_col)) stats::median(safe_numeric_10i(macro[[prob_col]]), na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    return(velocity_empty_df_10(c("pair_id", "split_value", "cell_n", "fate_state_n", "top_fate_state", "top_fate_cell_n", "max_fate_probability_median")))
  }
  dplyr::bind_rows(rows)
}

build_pair_review_lines_10i <- function(pair_ids, qc_df, driver_index, consistency_summary, split_compare, cellrank_summary, root_terminal_index) {
  if (length(pair_ids) == 0) {
    return("No velocity rows available.")
  }
  lines <- character()
  for (pair_id in pair_ids) {
    lines <- c(lines, sprintf("### %s", pair_id))
    qc <- qc_df[qc_df$pair_id == pair_id, , drop = FALSE]
    if (nrow(qc) > 0) {
      qc$velocity_confidence_mean <- safe_numeric_10i(qc$velocity_confidence_mean)
      lines <- c(lines, sprintf("- scVelo QC units: %s; median confidence: %s.", nrow(qc), fmt_num(stats::median(qc$velocity_confidence_mean, na.rm = TRUE), 3)))
    } else {
      lines <- c(lines, "- scVelo QC units: no data.")
    }
    drivers <- driver_index[driver_index$pair_id == pair_id, , drop = FALSE]
    if (nrow(drivers) > 0) {
      drivers$driver_n <- suppressWarnings(as.integer(drivers$driver_n))
      lines <- c(lines, sprintf("- Driver rows: %s; max driver count: %s.", nrow(drivers), max(drivers$driver_n, na.rm = TRUE)))
    }
    fate <- cellrank_summary[cellrank_summary$pair_id == pair_id, , drop = FALSE]
    if (nrow(fate) > 0) {
      lines <- c(lines, sprintf("- CellRank fate: %s unit(s); top fate states: %s.", nrow(fate), paste(unique(fate$top_fate_state), collapse = ", ")))
    }
    cons <- consistency_summary[consistency_summary$pair_id == pair_id & consistency_summary$status == "ok", , drop = FALSE]
    if (nrow(cons) > 0) {
      lines <- c(lines, sprintf("- Velocity/pseudotime consistency metrics available: %s.", nrow(cons)))
    }
    split <- split_compare[split_compare$pair_id == pair_id, , drop = FALSE]
    if (nrow(split) > 0) {
      lines <- c(
        lines,
        "- Cross-condition velocity:",
        render_markdown_table_local(utils::head(split, 20))
      )
    }
    rt <- root_terminal_index[root_terminal_index$pair_id == pair_id, , drop = FALSE]
    if (nrow(rt) > 0) {
      lines <- c(lines, sprintf("- Root/terminal bridge rows: %s; statuses: %s.", nrow(rt), paste(unique(rt$status), collapse = ", ")))
    }
    lines <- c(lines, "")
  }
  lines
}

loom_index <- read_tsv_with_cols_10i(
  read_manifest_output_or_fallback_10i(cfg$module_10a_manifest_path, "velocity_loom_index", cfg$velocity_loom_index_tsv),
  velocity_loom_index_cols_10
)
reference_index <- velocity_read_reference_index_10(cfg, include_not_ok = TRUE)
scvelo_index <- velocity_read_scvelo_index_10(cfg, include_not_ok = TRUE)
scvelo_qc <- read_tsv_with_cols_10i(
  read_manifest_output_or_fallback_10i(cfg$module_10c_manifest_path, "velocity_qc_tsv", cfg$velocity_scvelo_qc_tsv),
  velocity_scvelo_qc_cols_10
)
steady_index <- method_index_10i(cfg$module_10d_manifest_path, "velocyto_steady_index_tsv", cfg$velocity_velocyto_steady_index_tsv, velocity_velocyto_steady_index_cols_10)
driver_index <- method_index_10i(cfg$module_10e_manifest_path, "velocity_driver_index_tsv", cfg$velocity_driver_index_tsv, velocity_driver_index_cols_10)
driver_overlap <- read_tsv_with_cols_10i(
  read_manifest_output_or_fallback_10i(cfg$module_10e_manifest_path, "velocity_driver_overlap_tsv", cfg$velocity_driver_overlap_tsv),
  velocity_driver_overlap_cols_10
)
cellrank_index <- velocity_read_cellrank_index_10(cfg, include_not_ok = TRUE)
consistency_index <- method_index_10i(cfg$module_10g_manifest_path, "velocity_consistency_index_tsv", cfg$velocity_consistency_index_tsv, velocity_consistency_index_cols_10)
consistency_summary <- read_tsv_with_cols_10i(
  read_manifest_output_or_fallback_10i(cfg$module_10g_manifest_path, "velocity_consistency_summary_tsv", cfg$velocity_consistency_summary_tsv),
  velocity_consistency_summary_cols_10
)
split_compare <- read_tsv_with_cols_10i(
  read_manifest_output_or_fallback_10i(cfg$module_10g_manifest_path, "velocity_split_compare_tsv", cfg$velocity_split_compare_tsv),
  velocity_split_compare_cols_10
)
root_terminal_index <- method_index_10i(cfg$module_10h_manifest_path, "velocity_root_terminal_index_tsv", cfg$velocity_root_terminal_index_tsv, velocity_root_terminal_index_cols_10)

module_status <- dplyr::bind_rows(
  status_from_loom_10i(loom_index),
  status_from_reference_10i(reference_index),
  status_from_method_index_10i(scvelo_index, "scvelo_dynamical"),
  status_from_method_index_10i(steady_index, "velocyto_steady_state"),
  status_from_method_index_10i(driver_index, "scvelo_drivers", metric_value_from_driver_10i),
  status_from_method_index_10i(cellrank_index, "cellrank", metric_value_from_cellrank_10i),
  status_from_method_index_10i(consistency_index, "velocity_consistency", metric_value_from_consistency_index_10i),
  status_from_method_index_10i(root_terminal_index, "velocity_root_terminal", metric_value_from_root_terminal_10i)
)
if (nrow(module_status) == 0) {
  module_status <- empty_module_status_10i()
}
for (col in velocity_module_status_cols_10) {
  if (!col %in% colnames(module_status)) module_status[[col]] <- ""
}
module_status <- module_status[, velocity_module_status_cols_10, drop = FALSE]
write_tsv_local(module_status, cfg$velocity_module_status_tsv)

triage_rows <- list()
not_ok <- module_status[
  !module_status$status %in% c("ok", "ok_existing", "skipped_disabled") &
    nzchar(module_status$status),
  ,
  drop = FALSE
]
for (i in seq_len(nrow(not_ok))) {
  row <- not_ok[i, , drop = FALSE]
  triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
    sample_id = row$pair_id[[1]],
    severity = if (grepl("^failed", row$status[[1]])) "error" else "warning",
    signal_id = "velocity_module_status_not_ok",
    evidence = sprintf("%s value=%s status=%s notes=%s", row$method[[1]], row$value[[1]], row$status[[1]], row$notes[[1]]),
    recommended_action = "Review the module-specific velocity index and decide whether to rerun, waive, or interpret the row as exploratory."
  )
}

qc_summary <- metric_value_from_qc_10i(scvelo_qc)
if (nrow(qc_summary) > 0) {
  qc_summary$velocity_confidence_mean_num <- safe_numeric_10i(qc_summary$velocity_confidence_mean)
  qc_summary$unspliced_spliced_ratio_num <- safe_numeric_10i(qc_summary$unspliced_spliced_ratio)
  low_conf <- qc_summary[is.finite(qc_summary$velocity_confidence_mean_num) & qc_summary$velocity_confidence_mean_num < cfg$cellrank_min_velocity_confidence, , drop = FALSE]
  for (i in seq_len(nrow(low_conf))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = low_conf$pair_id[[i]],
      severity = "warning",
      signal_id = "velocity_low_confidence",
      evidence = sprintf("split=%s velocity_confidence_mean=%s threshold=%s", low_conf$split_value[[i]], low_conf$velocity_confidence_mean[[i]], cfg$cellrank_min_velocity_confidence),
      recommended_action = "Treat CellRank fate and velocity directionality as exploratory for this unit unless the low confidence is expected."
    )
  }
  no_unspliced <- qc_summary[is.finite(qc_summary$unspliced_spliced_ratio_num) & qc_summary$unspliced_spliced_ratio_num <= 0, , drop = FALSE]
  for (i in seq_len(nrow(no_unspliced))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = no_unspliced$pair_id[[i]],
      severity = "warning",
      signal_id = "velocity_no_unspliced_signal",
      evidence = sprintf("split=%s unspliced_spliced_ratio=%s", no_unspliced$split_value[[i]], no_unspliced$unspliced_spliced_ratio[[i]]),
      recommended_action = "Check loom generation and GTF compatibility before interpreting RNA velocity for this unit."
    )
  }
}

if (nrow(consistency_summary) > 0) {
  consistency_summary$value_num <- safe_numeric_10i(consistency_summary$value)
  low_consistency <- consistency_summary[
    consistency_summary$status == "ok" &
      consistency_summary$metric_id %in% c("dynamical_vs_stochastic_cosine", "dynamical_vs_velocyto_length_spearman", "latent_time_vs_slingshot_spearman") &
      is.finite(consistency_summary$value_num) &
      consistency_summary$value_num < 0.3,
    ,
    drop = FALSE
  ]
  for (i in seq_len(nrow(low_consistency))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = low_consistency$pair_id[[i]],
      severity = "warning",
      signal_id = "velocity_consistency_low",
      evidence = sprintf("split=%s metric=%s value=%s", low_consistency$split_value[[i]], low_consistency$metric_id[[i]], low_consistency$value[[i]]),
      recommended_action = "Do not use velocity as directional support without checking root selection, loom quality, and trajectory method agreement."
    )
  }
}

split_not_ok <- split_compare[nzchar(split_compare$status) & split_compare$status != "ok", , drop = FALSE]
for (i in seq_len(nrow(split_not_ok))) {
  triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
    sample_id = split_not_ok$pair_id[[i]],
    severity = if (grepl("^failed", split_not_ok$status[[i]])) "error" else "warning",
    signal_id = "velocity_split_compare_not_ok",
    evidence = sprintf("comparison=%s metric=%s status=%s reason=%s", split_not_ok$comparison_id[[i]], split_not_ok$metric_id[[i]], split_not_ok$status[[i]], split_not_ok$reason[[i]]),
    recommended_action = "Inspect 10g split comparison before making cross-condition velocity claims."
  )
}

triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
for (col in colnames(empty_triage_df(include_sample = TRUE))) {
  if (!col %in% colnames(triage_df)) triage_df[[col]] <- ""
}
triage_df <- triage_df[, colnames(empty_triage_df(include_sample = TRUE)), drop = FALSE]
write_tsv_local(triage_df, cfg$velocity_triage_tsv)

status_counts <- if (nrow(module_status) > 0) {
  as.data.frame(table(module_status$method, module_status$status), stringsAsFactors = FALSE)
} else {
  data.frame(method = character(0), status = character(0), n = integer(0), stringsAsFactors = FALSE)
}
if (nrow(status_counts) > 0) colnames(status_counts) <- c("method", "status", "n")

cellrank_summary <- summarize_cellrank_fate_10i(cellrank_index)
pair_ids <- unique(c(reference_index$pair_id, scvelo_index$pair_id, consistency_summary$pair_id, split_compare$pair_id))
pair_ids <- pair_ids[nzchar(pair_ids) & pair_ids != "__PROJECT__"]
pair_review <- build_pair_review_lines_10i(pair_ids, qc_summary, driver_index, consistency_summary, split_compare, cellrank_summary, root_terminal_index)

ok_n <- sum(module_status$status %in% c("ok", "ok_existing"), na.rm = TRUE)
not_ok_n <- sum(!module_status$status %in% c("ok", "ok_existing", "skipped_disabled"), na.rm = TRUE)

report_lines <- build_report_lines_v04(
  title = "10 RNA Velocity EDA",
  header_bullets = c(
    sprintf("Module status rows: %s ok/existing, %s non-ok requiring review.", ok_n, not_ok_n),
    sprintf("Triage signals collected: %s.", nrow(triage_df)),
    "Review velocity_finalize only after loom/reference inputs, scVelo quality, CellRank fate, velocity-pseudotime consistency, and split-mode velocity comparisons are checked."
  ),
  key_files = list(
    module_status = cfg$velocity_module_status_tsv,
    triage = cfg$velocity_triage_tsv,
    report = cfg$velocity_report_md,
    scvelo_qc = cfg$velocity_scvelo_qc_tsv,
    consistency = cfg$velocity_consistency_summary_tsv,
    split_compare = cfg$velocity_split_compare_tsv
  ),
  review_focus = c(
    "Approve velocity_finalize only after failed rows are waived or rerun.",
    "Use CellRank fate probabilities only for units with adequate cell count and velocity confidence.",
    "For split rows, inspect the cross-condition velocity section before claiming syf/f5 fate or flow-field differences."
  ),
  extra_sections = list(
    "Module Status Counts" = render_markdown_table_local(status_counts),
    "Module Status Preview" = render_markdown_table_local(utils::head(module_status, 100)),
    "scVelo QC Summary" = render_markdown_table_local(utils::head(qc_summary[, setdiff(colnames(qc_summary), c("velocity_confidence_mean_num", "unspliced_spliced_ratio_num")), drop = FALSE], 80)),
    "Driver Overlap" = render_markdown_table_local(driver_overlap),
    "CellRank Fate Summary" = render_markdown_table_local(cellrank_summary),
    "Velocity Consistency Metrics" = render_markdown_table_local(utils::head(consistency_summary[, setdiff(colnames(consistency_summary), "value_num"), drop = FALSE], 100)),
    "Split-Mode Cross-Condition Velocity" = render_markdown_table_local(split_compare),
    "Root/Terminal Bridge" = render_markdown_table_local(root_terminal_index),
    "Per Row Review" = pair_review
  ),
  triage_df = triage_df
)

ensure_dir(dirname(cfg$velocity_report_md))
write_markdown_local(report_lines, cfg$velocity_report_md)

outputs <- list(
  report_md = build_output_entry(cfg$velocity_report_md, "md", module_name, "final 10 RNA velocity EDA report", base_dir = cfg$project_root),
  velocity_module_status_tsv = build_output_entry(cfg$velocity_module_status_tsv, "tsv", module_name, "10 velocity module status table", base_dir = cfg$project_root, schema = infer_schema_from_df(module_status)),
  velocity_triage_tsv = build_output_entry(cfg$velocity_triage_tsv, "tsv", module_name, "10 velocity triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
)

velocity_manifest_10(
  cfg,
  cfg$module_10i_manifest_path,
  module_name,
  outputs = outputs,
  inputs = list(
    module_10a = cfg$module_10a_manifest_path,
    module_10b = cfg$module_10b_manifest_path,
    module_10c = cfg$module_10c_manifest_path,
    module_10d = cfg$module_10d_manifest_path,
    module_10e = cfg$module_10e_manifest_path,
    module_10f = cfg$module_10f_manifest_path,
    module_10g = cfg$module_10g_manifest_path,
    module_10h = cfg$module_10h_manifest_path
  ),
  depends_on = list(
    module_10c = cfg$module_10c_manifest_path,
    module_10d = cfg$module_10d_manifest_path,
    module_10e = cfg$module_10e_manifest_path,
    module_10f = cfg$module_10f_manifest_path,
    module_10g = cfg$module_10g_manifest_path,
    module_10h = cfg$module_10h_manifest_path
  )
)

message("10i completed. velocity report: ", cfg$velocity_report_md)
