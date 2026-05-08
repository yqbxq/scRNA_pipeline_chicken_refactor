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

load_required_packages(c("Seurat", "dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09j_trajectory_velocity_link"
prepare_dirs_09(cfg)

velocity_link_index_cols <- c(
  trajectory_method_index_cols_09,
  "velocity_h5ad", "velocity_metric", "spearman_rho", "p_value",
  "common_cell_n", "root_mean_velocity_time", "terminal_mean_velocity_time",
  "direction_aligned"
)

read_slingshot_index_09j <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09d_manifest_path,
    "slingshot_index_tsv",
    cfg$trajectory_slingshot_index_tsv,
    cols = c(trajectory_method_index_cols_09, "fine_output_path", "lineage_count", "root_label", "label_var")
  )
}

slingshot_row_for_unit_09j <- function(index_df, pair_id, split_value) {
  hit <- index_df[
    index_df$pair_id == pair_id &
      index_df$split_value == display_scalar_value(split_value, "pooled"),
    ,
    drop = FALSE
  ]
  if (nrow(hit) == 0) NULL else hit[1, , drop = FALSE]
}

velocity_h5ad_candidates_09j <- function(cfg, pair_id, split_value) {
  if (!dir.exists(cfg$velocity_output_dir)) {
    return(character(0))
  }
  files <- list.files(cfg$velocity_output_dir, pattern = "\\.h5ad$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    return(character(0))
  }
  value <- trajectory_unit_value_09(split_value)
  tokens <- unique(c(safe_id_09(pair_id), pair_id, if (nzchar(value)) c(safe_id_09(value), value) else character(0)))
  score <- vapply(files, function(path) {
    lower <- tolower(path)
    sum(vapply(tokens, function(tok) grepl(tolower(tok), lower, fixed = TRUE), logical(1)))
  }, numeric(1))
  files[order(score, decreasing = TRUE)]
}

extract_velocity_time_09j <- function(h5ad_path, out_tsv) {
  python <- trajectory_python_exe_09()
  if (!nzchar(python) || !file.exists(python)) {
    stop("python executable for velocity/anndata is unavailable", call. = FALSE)
  }
  script <- tempfile(fileext = ".py")
  writeLines(c(
    "import sys",
    "import anndata as ad",
    "import pandas as pd",
    "h5ad_path, out_tsv = sys.argv[1], sys.argv[2]",
    "adata = ad.read_h5ad(h5ad_path)",
    "obs = adata.obs.copy()",
    "cols = [c for c in ['latent_time', 'velocity_pseudotime', 'velocity_length', 'velocity_confidence'] if c in obs.columns]",
    "df = obs[cols].copy() if cols else pd.DataFrame(index=obs.index)",
    "df.insert(0, 'cell_id', obs.index.astype(str))",
    "df.to_csv(out_tsv, sep='\\t', index=False)"
  ), script, useBytes = TRUE)
  status <- system2(python, args = c(script, h5ad_path, out_tsv), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(status, "status")) && attr(status, "status") != 0) {
    stop(paste(status, collapse = "\n"), call. = FALSE)
  }
  read_tsv_optional(out_tsv)
}

direction_summary_09j <- function(seu, unit, velocity_df, metric_col) {
  label_var <- normalize_scalar_value(unit$coarse_label_var[[1]])
  root <- trajectory_selected_root_09(cfg, unit$pair_id[[1]], unit$split_value[[1]], unit$root_group[[1]])
  terminal <- normalize_scalar_value(unit$terminal_group[[1]])
  if (!label_var %in% colnames(seu@meta.data) || !metric_col %in% colnames(velocity_df)) {
    return(list(root_mean = NA_real_, terminal_mean = NA_real_, aligned = "not_evaluated"))
  }
  root_cells <- trajectory_root_cells_09(seu, label_var, root)
  terminal_cells <- trajectory_terminal_cells_09(seu, label_var, terminal)
  if (length(root_cells) == 0 || length(terminal_cells) == 0) {
    return(list(root_mean = NA_real_, terminal_mean = NA_real_, aligned = "not_evaluated"))
  }
  values <- setNames(suppressWarnings(as.numeric(velocity_df[[metric_col]])), velocity_df$cell_id)
  root_values <- values[intersect(root_cells, names(values))]
  terminal_values <- values[intersect(terminal_cells, names(values))]
  root_mean <- mean(root_values, na.rm = TRUE)
  terminal_mean <- mean(terminal_values, na.rm = TRUE)
  aligned <- if (is.finite(root_mean) && is.finite(terminal_mean)) {
    if (terminal_mean > root_mean) "yes" else "no"
  } else {
    "not_evaluated"
  }
  list(root_mean = root_mean, terminal_mean = terminal_mean, aligned = aligned)
}

run_velocity_link_09j <- function(unit, sling_row) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  link_tsv <- trajectory_method_output_path_09(cfg, "velocity_link", pair_id, split_value, "velocity_link", "tsv")
  velocity_obs_tsv <- trajectory_method_output_path_09(cfg, "velocity_link", pair_id, split_value, "velocity_obs_time", "tsv")
  fig_png <- trajectory_method_figure_path_09(cfg, "velocity_link", pair_id, split_value, "Figure_VelocityLink")
  ensure_dir(dirname(link_tsv))
  ensure_dir(dirname(fig_png))
  started <- proc.time()[["elapsed"]]

  tryCatch({
    candidates <- velocity_h5ad_candidates_09j(cfg, pair_id, split_value)
    if (length(candidates) == 0) {
      stop("skipped_no_velocity: no .h5ad files in VELOCITY_OUTPUT_DIR", call. = FALSE)
    }
    h5ad_path <- candidates[[1]]
    if (is.null(sling_row) || sling_row$status[[1]] != "ok" || !file.exists(sling_row$output_path[[1]])) {
      stop("missing successful Slingshot pseudotime", call. = FALSE)
    }
    velocity_df <- extract_velocity_time_09j(h5ad_path, velocity_obs_tsv)
    slingshot_df <- trajectory_read_pseudotime_csv_09(sling_row$output_path[[1]], "slingshot")
    metric_cols <- intersect(c("latent_time", "velocity_pseudotime"), colnames(velocity_df))
    if (length(metric_cols) == 0) {
      stop("velocity h5ad lacks latent_time and velocity_pseudotime in obs", call. = FALSE)
    }
    seu <- readRDS(unit$input_rds[[1]])
    rows <- list()
    for (metric_col in metric_cols) {
      merged <- merge(
        slingshot_df[, c("cell_id", "pseudotime"), drop = FALSE],
        velocity_df[, c("cell_id", metric_col), drop = FALSE],
        by = "cell_id",
        all = FALSE
      )
      merged$pseudotime <- suppressWarnings(as.numeric(merged$pseudotime))
      merged$velocity_time <- suppressWarnings(as.numeric(merged[[metric_col]]))
      merged <- merged[is.finite(merged$pseudotime) & is.finite(merged$velocity_time), , drop = FALSE]
      if (nrow(merged) >= 3) {
        ct <- suppressWarnings(stats::cor.test(merged$pseudotime, merged$velocity_time, method = "spearman", exact = FALSE))
        rho <- unname(ct$estimate)
        p_value <- ct$p.value
      } else {
        rho <- NA_real_
        p_value <- NA_real_
      }
      direction <- direction_summary_09j(seu, unit, velocity_df, metric_col)
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = pair_id,
        split_value = split_value,
        velocity_metric = metric_col,
        velocity_h5ad = h5ad_path,
        common_cell_n = nrow(merged),
        spearman_rho = rho,
        p_value = p_value,
        root_mean_velocity_time = direction$root_mean,
        terminal_mean_velocity_time = direction$terminal_mean,
        direction_aligned = direction$aligned,
        status = if (nrow(merged) >= 3) "ok" else "skipped_low_common_cells",
        reason = if (nrow(merged) >= 3) "" else "fewer than 3 common finite cells",
        stringsAsFactors = FALSE
      )
    }
    result_df <- dplyr::bind_rows(rows)
    write_tsv_local(result_df, link_tsv)
    plot_df <- result_df[is.finite(result_df$spearman_rho), , drop = FALSE]
    if (nrow(plot_df) == 0) {
      plot_obj <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No velocity link correlation") +
        ggplot2::theme_void()
    } else {
      plot_obj <- ggplot2::ggplot(plot_df, ggplot2::aes(x = velocity_metric, y = spearman_rho, fill = direction_aligned)) +
        ggplot2::geom_col() +
        ggplot2::theme_bw(base_size = 9) +
        ggplot2::labs(title = sprintf("Velocity link %s %s", pair_id, split_value), x = NULL, y = "Spearman rho")
    }
    save_plot_local(plot_obj, fig_png, width = 5.5, height = 4)
    best <- result_df[1, , drop = FALSE]
    data.frame(
      pair_id = pair_id, split_value = split_value, method = "velocity_link", methods_enabled = "yes",
      input_rds = unit$input_rds[[1]], output_path = link_tsv, extra_path = velocity_obs_tsv,
      figure_path = fig_png, n_cells = best$common_cell_n[[1]], status = best$status[[1]],
      reason = best$reason[[1]], runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      velocity_h5ad = best$velocity_h5ad[[1]], velocity_metric = best$velocity_metric[[1]],
      spearman_rho = best$spearman_rho[[1]], p_value = best$p_value[[1]],
      common_cell_n = best$common_cell_n[[1]], root_mean_velocity_time = best$root_mean_velocity_time[[1]],
      terminal_mean_velocity_time = best$terminal_mean_velocity_time[[1]],
      direction_aligned = best$direction_aligned[[1]], stringsAsFactors = FALSE
    )
  }, error = function(e) {
    status <- if (grepl("skipped_no_velocity", conditionMessage(e))) {
      "skipped_no_velocity"
    } else if (grepl("anndata|python executable", conditionMessage(e), ignore.case = TRUE)) {
      "skipped_h5ad_reader_unavailable"
    } else if (grepl("Slingshot", conditionMessage(e))) {
      "skipped_no_slingshot"
    } else {
      "failed"
    }
    result_df <- data.frame(
      pair_id = pair_id,
      split_value = split_value,
      velocity_metric = "",
      velocity_h5ad = "",
      common_cell_n = 0L,
      spearman_rho = NA_real_,
      p_value = NA_real_,
      root_mean_velocity_time = NA_real_,
      terminal_mean_velocity_time = NA_real_,
      direction_aligned = "not_evaluated",
      status = status,
      reason = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    write_tsv_local(result_df, link_tsv)
    data.frame(
      pair_id = pair_id, split_value = split_value, method = "velocity_link", methods_enabled = "yes",
      input_rds = unit$input_rds[[1]], output_path = link_tsv, extra_path = velocity_obs_tsv,
      figure_path = fig_png, n_cells = 0L, status = status, reason = conditionMessage(e),
      runtime_s = round(proc.time()[["elapsed"]] - started, 3),
      velocity_h5ad = "", velocity_metric = "", spearman_rho = NA_real_, p_value = NA_real_,
      common_cell_n = 0L, root_mean_velocity_time = NA_real_, terminal_mean_velocity_time = NA_real_,
      direction_aligned = "not_evaluated", stringsAsFactors = FALSE
    )
  })
}

units <- trajectory_execution_units_09(cfg)
slingshot_index <- read_slingshot_index_09j(cfg)
index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(units))) {
  unit <- units[i, , drop = FALSE]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  row <- run_velocity_link_09j(unit, slingshot_row_for_unit_09j(slingshot_index, unit$pair_id[[1]], split_value))
  index_rows[[length(index_rows) + 1L]] <- row
  dynamic_outputs[[sprintf("velocity_link_%s", trajectory_unit_file_id_09(unit$pair_id[[1]], split_value))]] <- build_output_entry(
    row$output_path[[1]],
    "tsv",
    module_name,
    "pseudotime versus RNA velocity latent-time link",
    base_dir = cfg$project_root
  )
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else trajectory_empty_df_09(velocity_link_index_cols)
write_tsv_local(index_df, cfg$trajectory_velocity_link_index_tsv)

outputs <- c(
  list(velocity_link_index_tsv = build_output_entry(cfg$trajectory_velocity_link_index_tsv, "tsv", module_name, "velocity link status by pair/split", base_dir = cfg$project_root, schema = infer_schema_from_df(index_df))),
  dynamic_outputs
)
trajectory_method_manifest_09(
  cfg,
  cfg$module_09j_manifest_path,
  module_name,
  outputs,
  inputs = list(module_09a = cfg$module_09a_manifest_path, module_09d = cfg$module_09d_manifest_path, velocity_output_dir = cfg$velocity_output_dir),
  depends_on = list(module_09a = cfg$module_09a_manifest_path, module_09d = cfg$module_09d_manifest_path)
)

message("09j completed. velocity link index: ", cfg$trajectory_velocity_link_index_tsv)
