#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")

load_required_packages(c("dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09l_trajectory_eda"
prepare_dirs_09(cfg)

module_status_cols_09l <- c(
  "pair_id", "value", "method", "methods_enabled", "status",
  "n_cells_used", "runtime_s", "notes"
)

empty_module_status_09l <- function() {
  data.frame(
    pair_id = character(0),
    value = character(0),
    method = character(0),
    methods_enabled = character(0),
    status = character(0),
    n_cells_used = integer(0),
    runtime_s = numeric(0),
    notes = character(0),
    stringsAsFactors = FALSE
  )
}

read_manifest_output_or_fallback_09l <- function(manifest_path, key, fallback_path) {
  path <- trajectory_manifest_output_optional_09(manifest_path, key)
  if (!nzchar(path)) fallback_path else path
}

read_tsv_with_cols_09l <- function(path, cols) {
  df <- read_tsv_optional(path)
  if (nrow(df) == 0) {
    return(trajectory_empty_df_09(cols))
  }
  for (col in cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- ""
    }
  }
  df
}

method_index_09l <- function(manifest_path, key, fallback_path, extra_cols = character()) {
  trajectory_read_method_index_09(
    manifest_path,
    key,
    fallback_path,
    cols = unique(c(trajectory_method_index_cols_09, extra_cols))
  )
}

status_from_method_index_09l <- function(index_df, fallback_method) {
  if (nrow(index_df) == 0) {
    return(empty_module_status_09l())
  }
  rows <- lapply(seq_len(nrow(index_df)), function(i) {
    row <- index_df[i, , drop = FALSE]
    method <- normalize_scalar_value(row$method[[1]], fallback_method)
    if ("option" %in% colnames(row) && nzchar(normalize_scalar_value(row$option[[1]]))) {
      method <- paste(method, row$option[[1]], sep = ":")
    }
    notes <- normalize_scalar_value(row$reason[[1]])
    if (!nzchar(notes) && nzchar(normalize_scalar_value(row$output_path[[1]]))) {
      notes <- normalize_scalar_value(row$output_path[[1]])
    }
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = method,
      methods_enabled = display_scalar_value(row$methods_enabled[[1]], "yes"),
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = suppressWarnings(as.integer(row$n_cells[[1]])),
      runtime_s = suppressWarnings(as.numeric(row$runtime_s[[1]])),
      notes = notes,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_input_09l <- function(input_status) {
  if (nrow(input_status) == 0) {
    return(empty_module_status_09l())
  }
  rows <- lapply(seq_len(nrow(input_status)), function(i) {
    row <- input_status[i, , drop = FALSE]
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = "input",
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = suppressWarnings(as.integer(row$cell_n_after[[1]])),
      runtime_s = NA_real_,
      notes = normalize_scalar_value(row$reason[[1]]),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_root_09l <- function(root_index) {
  if (nrow(root_index) == 0) {
    return(empty_module_status_09l())
  }
  rows <- lapply(seq_len(nrow(root_index)), function(i) {
    row <- root_index[i, , drop = FALSE]
    notes <- paste(
      c(
        sprintf("selected_root=%s", display_scalar_value(row$selected_root[[1]], "NA")),
        sprintf("source=%s", display_scalar_value(row$selected_source[[1]], "NA")),
        sprintf("non_prior_vote=%s", display_scalar_value(row$non_prior_vote_root[[1]], "NA")),
        sprintf("disagreeing_methods=%s", display_scalar_value(row$disagreeing_methods[[1]], "none")),
        sprintf("velocity_status=%s", display_scalar_value(row$velocity_status[[1]], "NA")),
        normalize_scalar_value(row$reason[[1]])
      ),
      collapse = "; "
    )
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = "root",
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = NA_integer_,
      runtime_s = NA_real_,
      notes = notes,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_figures_09l <- function(figure_index) {
  if (nrow(figure_index) == 0) {
    return(empty_module_status_09l())
  }
  rows <- lapply(seq_len(nrow(figure_index)), function(i) {
    row <- figure_index[i, , drop = FALSE]
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_value[[1]], "pooled"),
      method = paste("figures", display_scalar_value(row$figure_role[[1]], "summary"), sep = ":"),
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = suppressWarnings(as.integer(row$n_cells[[1]])),
      runtime_s = NA_real_,
      notes = normalize_scalar_value(row$reason[[1]], normalize_scalar_value(row$png_path[[1]])),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

status_from_split_compare_09l <- function(split_index) {
  if (nrow(split_index) == 0) {
    return(empty_module_status_09l())
  }
  rows <- lapply(seq_len(nrow(split_index)), function(i) {
    row <- split_index[i, , drop = FALSE]
    data.frame(
      pair_id = normalize_scalar_value(row$pair_id[[1]], "__PROJECT__"),
      value = display_scalar_value(row$split_values[[1]], "NA"),
      method = "split_compare",
      methods_enabled = "yes",
      status = display_scalar_value(row$status[[1]], "not_available"),
      n_cells_used = NA_integer_,
      runtime_s = NA_real_,
      notes = normalize_scalar_value(row$reason[[1]], normalize_scalar_value(row$figure_path[[1]])),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

input_status <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09b_manifest_path, "input_status_tsv", cfg$trajectory_input_status_tsv),
  c("pair_id", "split_value", "cell_n_after", "status", "reason")
)
root_index <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09c_manifest_path, "trajectory_root_index", cfg$trajectory_root_index_tsv),
  c(
    "pair_id", "split_value", "selected_root", "selected_source", "status", "reason",
    "non_prior_vote_root", "non_prior_vote_sources", "prior_disagreement", "disagreeing_methods",
    "velocity_status"
  )
)
root_agreement <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09c_manifest_path, "split_root_agreement", cfg$trajectory_root_split_agreement_tsv),
  c("pair_id", "split_values", "selected_roots", "status", "reason")
)
slingshot_index <- method_index_09l(cfg$module_09d_manifest_path, "slingshot_index_tsv", cfg$trajectory_slingshot_index_tsv, c("lineage_count", "root_label", "label_var"))
monocle3_index <- method_index_09l(cfg$module_09e_manifest_path, "monocle3_index_tsv", cfg$trajectory_monocle3_index_tsv, c("root_cell_n", "root_label"))
monocle2_index <- method_index_09l(cfg$module_09e2_manifest_path, "monocle2_index_tsv", cfg$trajectory_monocle2_index_tsv)
paga_index <- method_index_09l(cfg$module_09f_manifest_path, "paga_dpt_index_tsv", cfg$trajectory_paga_index_tsv, c("h5ad_path"))
palantir_index <- method_index_09l(cfg$module_09g_manifest_path, "palantir_index_tsv", cfg$trajectory_palantir_index_tsv, c("fate_path"))
tradeseq_index <- method_index_09l(
  cfg$module_09h_manifest_path,
  "tradeseq_index_tsv",
  cfg$trajectory_tradeseq_index_tsv,
  c("option", "association_path", "pattern_path", "driver_path", "conditiontest_path")
)
consensus_index <- method_index_09l(
  cfg$module_09i_manifest_path,
  "consensus_index_tsv",
  cfg$trajectory_consensus_index_tsv,
  c("method_count", "methods_used", "correlation_path", "jaccard_path", "consensus_pseudotime_path", "conflict_count")
)
velocity_index <- method_index_09l(
  cfg$module_09j_manifest_path,
  "velocity_link_index_tsv",
  cfg$trajectory_velocity_link_index_tsv,
  c("velocity_h5ad", "velocity_metric", "spearman_rho", "p_value", "common_cell_n", "direction_aligned")
)
figure_index <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09k_manifest_path, "figure_index_tsv", cfg$trajectory_figure_index_tsv),
  c("pair_id", "split_value", "figure_role", "color_var", "png_path", "n_cells", "status", "reason")
)
split_compare_index <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09m_manifest_path, "split_compare_index_tsv", cfg$trajectory_split_compare_index_tsv),
  c("pair_id", "split_values", "figure_path", "status", "reason")
)

module_status <- dplyr::bind_rows(
  status_from_input_09l(input_status),
  status_from_root_09l(root_index),
  status_from_method_index_09l(slingshot_index, "slingshot"),
  status_from_method_index_09l(monocle3_index, "monocle3"),
  status_from_method_index_09l(monocle2_index, "monocle2"),
  status_from_method_index_09l(paga_index, "paga_dpt"),
  status_from_method_index_09l(palantir_index, "palantir"),
  status_from_method_index_09l(tradeseq_index, "tradeseq"),
  status_from_method_index_09l(consensus_index, "consensus"),
  status_from_method_index_09l(velocity_index, "velocity_link"),
  status_from_figures_09l(figure_index),
  status_from_split_compare_09l(split_compare_index)
)
if (nrow(module_status) == 0) {
  module_status <- empty_module_status_09l()
}
module_status <- module_status[, module_status_cols_09l, drop = FALSE]
write_tsv_local(module_status, cfg$trajectory_module_status_tsv)

triage_rows <- list()
input_triage <- read_tsv_optional(cfg$trajectory_inputs_triage_tsv)
if (nrow(input_triage) > 0) {
  triage_rows[[length(triage_rows) + 1L]] <- input_triage
}
consensus_triage <- read_tsv_optional(cfg$trajectory_consensus_triage_tsv)
if (nrow(consensus_triage) > 0) {
  triage_rows[[length(triage_rows) + 1L]] <- consensus_triage
}

if (nrow(module_status) > 0) {
  issue_status <- module_status[
    grepl("^failed", module_status$status) |
      grepl("^warning", module_status$status) |
      module_status$status %in% c("warning_low_overlap", "failed_no_root"),
    ,
    drop = FALSE
  ]
  for (i in seq_len(nrow(issue_status))) {
    row <- issue_status[i, , drop = FALSE]
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = row$pair_id[[1]],
      severity = if (grepl("^failed", row$status[[1]])) "error" else "warning",
      signal_id = "trajectory_module_status_not_ok",
      evidence = sprintf("%s value=%s status=%s notes=%s", row$method[[1]], row$value[[1]], row$status[[1]], row$notes[[1]]),
      recommended_action = "Review the module-specific index and decide whether to rerun, waive, or interpret the result as exploratory."
    )
  }
}

if (nrow(root_index) > 0 && "prior_disagreement" %in% colnames(root_index)) {
  prior_disagree <- root_index[root_index$prior_disagreement == "yes", , drop = FALSE]
  for (i in seq_len(nrow(prior_disagree))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = prior_disagree$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_prior_root_disagreement",
      evidence = sprintf(
        "split=%s prior-selected=%s non_prior_vote=%s disagreeing_methods=%s",
        prior_disagree$split_value[[i]],
        prior_disagree$selected_root[[i]],
        prior_disagree$non_prior_vote_root[[i]],
        prior_disagree$disagreeing_methods[[i]]
      ),
      recommended_action = "Prior root_group is being honored, but inspect 09c root inference before treating directionality as method-supported."
    )
  }
}

if (nrow(root_index) > 0 && "velocity_status" %in% colnames(root_index)) {
  missing_velocity <- root_index[root_index$velocity_status %in% c("skipped_no_velocity", "skipped_empty_velocity"), , drop = FALSE]
  for (i in seq_len(nrow(missing_velocity))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = missing_velocity$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_velocity_support_missing",
      evidence = sprintf("split=%s velocity_status=%s", missing_velocity$split_value[[i]], missing_velocity$velocity_status[[i]]),
      recommended_action = "Run module 10 velocity first, then rerun 09c and 09j if velocity-supported root direction is required."
    )
  }
}

if (nrow(root_agreement) > 0) {
  disagree <- root_agreement[root_agreement$status != "ok", , drop = FALSE]
  for (i in seq_len(nrow(disagree))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = disagree$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_split_mode_root_disagreement",
      evidence = sprintf("split_values=%s selected_roots=%s status=%s", disagree$split_values[[i]], disagree$selected_roots[[i]], disagree$status[[i]]),
      recommended_action = "Inspect 09c root inference tables before comparing split pseudotime direction."
    )
  }
}

option_consistency <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09m_manifest_path, "option_a_vs_b_consistency_tsv", cfg$trajectory_option_a_b_consistency_tsv),
  c("pair_id", "split_values", "option_a_driver_n", "option_b_driver_n", "overlap_n", "jaccard", "status", "reason")
)
if (nrow(option_consistency) > 0) {
  low_overlap <- option_consistency[option_consistency$status == "warning_low_overlap", , drop = FALSE]
  for (i in seq_len(nrow(low_overlap))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = low_overlap$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_option_a_b_inconsistent",
      evidence = sprintf("jaccard=%s option_a=%s option_b=%s", low_overlap$jaccard[[i]], low_overlap$option_a_driver_n[[i]], low_overlap$option_b_driver_n[[i]]),
      recommended_action = "Do not use option A split-only drivers as condition-specific claims without checking option B conditionTest."
    )
  }
}

if (nrow(velocity_index) > 0 && "direction_aligned" %in% colnames(velocity_index)) {
  mismatch <- velocity_index[velocity_index$direction_aligned == "no", , drop = FALSE]
  for (i in seq_len(nrow(mismatch))) {
    triage_rows[[length(triage_rows) + 1L]] <- make_triage_row(
      sample_id = mismatch$pair_id[[i]],
      severity = "warning",
      signal_id = "trajectory_velocity_direction_mismatch",
      evidence = sprintf("split=%s metric=%s rho=%s", mismatch$split_value[[i]], mismatch$velocity_metric[[i]], mismatch$spearman_rho[[i]]),
      recommended_action = "Review root direction and RNA velocity latent-time before using velocity as directional support."
    )
  }
}

triage_df <- if (length(triage_rows) > 0) dplyr::bind_rows(triage_rows) else empty_triage_df(include_sample = TRUE)
for (col in colnames(empty_triage_df(include_sample = TRUE))) {
  if (!col %in% colnames(triage_df)) {
    triage_df[[col]] <- ""
  }
}
triage_df <- triage_df[, colnames(empty_triage_df(include_sample = TRUE)), drop = FALSE]
write_tsv_local(triage_df, cfg$trajectory_triage_tsv)

status_counts <- if (nrow(module_status) > 0) {
  as.data.frame(table(module_status$method, module_status$status), stringsAsFactors = FALSE)
} else {
  data.frame(method = character(0), status = character(0), n = integer(0), stringsAsFactors = FALSE)
}
if (nrow(status_counts) > 0) {
  colnames(status_counts) <- c("method", "status", "n")
}

consensus_corr <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09i_manifest_path, "consensus_correlation_summary_tsv", cfg$trajectory_consensus_correlation_summary_tsv),
  c("pair_id", "split_value", "method_a", "method_b", "spearman_rho", "kendall_tau", "p_value", "n_cells")
)
tradeseq_driver_index <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09h_manifest_path, "tradeseq_driver_index_tsv", cfg$trajectory_tradeseq_driver_index_tsv),
  c("pair_id", "split_value", "option", "driver_path", "n_drivers", "status", "reason")
)
entropy_compare <- read_tsv_with_cols_09l(
  read_manifest_output_or_fallback_09l(cfg$module_09m_manifest_path, "entropy_compare_tsv", cfg$trajectory_entropy_compare_tsv),
  c("pair_id", "split_a", "split_b", "label_var", "label_value", "n_cells_a", "n_cells_b", "median_entropy_a", "median_entropy_b", "status", "reason")
)

ok_n <- sum(module_status$status == "ok", na.rm = TRUE)
not_ok_n <- sum(module_status$status != "ok", na.rm = TRUE)
triage_n <- nrow(triage_df)

report_lines <- build_report_lines_v04(
  title = "09 Trajectory EDA",
  header_bullets = c(
    sprintf("Module status rows: %s ok, %s non-ok.", ok_n, not_ok_n),
    sprintf("Triage signals collected: %s.", triage_n),
    "Final trajectory interpretation should combine input QC, root selection, method agreement, tradeSeq drivers, Palantir entropy, Monocle information, velocity consistency, and split comparison."
  ),
  key_files = list(
    module_status = cfg$trajectory_module_status_tsv,
    triage = cfg$trajectory_triage_tsv,
    report = cfg$trajectory_report_md,
    figure_index = cfg$trajectory_figure_index_tsv,
    split_compare_report = cfg$trajectory_split_compare_report_md
  ),
  review_focus = c(
    "Approve trajectory_finalize only after failed rows are waived or rerun and split-mode warnings are reviewed.",
    "Use consensus pseudotime preferentially when at least two methods agree; otherwise treat single-method pseudotime as exploratory.",
    "Do not interpret option A split-only tradeSeq drivers as condition-specific without checking option B consistency."
  ),
  extra_sections = list(
    "Module Status Counts" = render_markdown_table_local(status_counts),
    "Module Status Preview" = render_markdown_table_local(utils::head(module_status, 80)),
    "Input QC Summary" = render_markdown_table_local(utils::head(input_status, 40)),
    "Root Summary" = render_markdown_table_local(root_index),
    "Method Pseudotime Correlation" = render_markdown_table_local(utils::head(consensus_corr, 60)),
    "Velocity Consistency" = render_markdown_table_local(utils::head(velocity_index[, intersect(c("pair_id", "split_value", "velocity_metric", "spearman_rho", "common_cell_n", "direction_aligned", "status", "reason"), colnames(velocity_index)), drop = FALSE], 60)),
    "tradeSeq Driver Counts" = render_markdown_table_local(tradeseq_driver_index),
    "Palantir Entropy Compare" = render_markdown_table_local(utils::head(entropy_compare, 60)),
    "Monocle 2 Info" = render_markdown_table_local(monocle2_index),
    "Split Comparison" = render_markdown_table_local(split_compare_index)
  ),
  triage_df = triage_df
)

ensure_dir(dirname(cfg$trajectory_report_md))
write_markdown_local(report_lines, cfg$trajectory_report_md)

outputs <- list(
  report_md = build_output_entry(cfg$trajectory_report_md, "md", module_name, "final 09 trajectory EDA report", base_dir = cfg$project_root),
  trajectory_module_status_tsv = build_output_entry(cfg$trajectory_module_status_tsv, "tsv", module_name, "09 module status table", base_dir = cfg$project_root, schema = infer_schema_from_df(module_status)),
  trajectory_triage_tsv = build_output_entry(cfg$trajectory_triage_tsv, "tsv", module_name, "09 trajectory triage signals", base_dir = cfg$project_root, schema = infer_schema_from_df(triage_df))
)

trajectory_method_manifest_09(
  cfg,
  cfg$module_09l_manifest_path,
  module_name,
  outputs,
  inputs = list(
    module_09a = cfg$module_09a_manifest_path,
    module_09b = cfg$module_09b_manifest_path,
    module_09c = cfg$module_09c_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09e = cfg$module_09e_manifest_path,
    module_09e2 = cfg$module_09e2_manifest_path,
    module_09f = cfg$module_09f_manifest_path,
    module_09g = cfg$module_09g_manifest_path,
    module_09h = cfg$module_09h_manifest_path,
    module_09i = cfg$module_09i_manifest_path,
    module_09j = cfg$module_09j_manifest_path,
    module_09k = cfg$module_09k_manifest_path,
    module_09m = cfg$module_09m_manifest_path
  ),
  depends_on = list(
    module_09b = cfg$module_09b_manifest_path,
    module_09c = cfg$module_09c_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09h = cfg$module_09h_manifest_path,
    module_09i = cfg$module_09i_manifest_path,
    module_09j = cfg$module_09j_manifest_path,
    module_09k = cfg$module_09k_manifest_path,
    module_09m = cfg$module_09m_manifest_path
  )
)

message("09l completed. report: ", cfg$trajectory_report_md)
