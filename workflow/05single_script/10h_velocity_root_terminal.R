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
module_name <- "10h_velocity_root_terminal"
prepare_dirs_10(cfg)

read_delim_optional_10h <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path) || !file.exists(path) || isTRUE(file.info(path)$size == 0)) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read_tsv_optional(path)
  }
}

scale01_10h <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  keep <- is.finite(x)
  if (!any(keep)) return(out)
  rng <- range(x[keep])
  if (diff(rng) == 0) {
    out[keep] <- 0.5
  } else {
    out[keep] <- (x[keep] - rng[[1]]) / diff(rng)
  }
  out
}

shrink_probability_10h <- function(prob, cell_n, min_cells) {
  prob <- suppressWarnings(as.numeric(prob))
  cell_n <- suppressWarnings(as.numeric(cell_n))
  out <- prob
  keep <- is.finite(prob) & is.finite(cell_n)
  if (!any(keep)) {
    return(out)
  }
  weight <- pmin(1, pmax(0, cell_n[keep] / max(1, min_cells)))
  out[keep] <- 0.5 + (prob[keep] - 0.5) * weight
  out
}

empty_index_row_10h <- function(unit, out_tsv, status, reason, runtime_s = 0) {
  data.frame(
    pair_id = unit$pair_id[[1]],
    split_value = display_scalar_value(unit$split_value[[1]], "pooled"),
    method = "velocity_root_terminal",
    methods_enabled = "yes",
    input_path = unit$output_path[[1]],
    output_path = if (identical(status, "ok")) out_tsv else "",
    extra_path = normalize_scalar_value(unit$vector_tsv[[1]]),
    figure_path = "",
    n_cells = suppressWarnings(as.integer(unit$n_cells[[1]])),
    status = status,
    reason = reason,
    runtime_s = runtime_s,
    root_terminal_tsv = out_tsv,
    stringsAsFactors = FALSE
  )
}

root_terminal_table_10h <- function(unit, cellrank_index) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  vectors <- read_delim_optional_10h(unit$vector_tsv[[1]])
  if (nrow(vectors) == 0) {
    stop("missing 10c vector_tsv", call. = FALSE)
  }
  label_col <- if ("cell_type" %in% colnames(vectors) && any(nzchar(vectors$cell_type))) "cell_type" else "cluster"
  if (!label_col %in% colnames(vectors)) {
    stop("vector_tsv lacks cell_type/cluster labels", call. = FALSE)
  }
  vectors$label <- as.character(vectors[[label_col]])
  vectors <- vectors[nzchar(vectors$label), , drop = FALSE]
  if (nrow(vectors) == 0) {
    stop("no labeled cells available for root/terminal scoring", call. = FALSE)
  }
  vectors$latent01 <- if ("latent_time" %in% colnames(vectors)) scale01_10h(vectors$latent_time) else NA_real_

  fate <- data.frame(stringsAsFactors = FALSE)
  cr <- cellrank_index[cellrank_index$pair_id == pair_id & cellrank_index$split_value == split_value & cellrank_index$status == "ok", , drop = FALSE]
  if (nrow(cr) > 0) {
    fate <- read_delim_optional_10h(cr$fate_csv[[1]])
  }
  if (nrow(fate) > 0 && "cell_id" %in% colnames(fate)) {
    fate_cols <- setdiff(colnames(fate), c("pair_id", "split_value", "cell_id"))
    fate[fate_cols] <- lapply(fate[fate_cols], function(x) suppressWarnings(as.numeric(x)))
    fate$terminal_probability_cell <- if (length(fate_cols) > 0) apply(fate[, fate_cols, drop = FALSE], 1, max, na.rm = TRUE) else NA_real_
    fate$terminal_probability_cell[!is.finite(fate$terminal_probability_cell)] <- NA_real_
    vectors <- merge(vectors, fate[, c("cell_id", "terminal_probability_cell"), drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
  } else {
    vectors$terminal_probability_cell <- NA_real_
  }

  grouped <- vectors |>
    dplyr::group_by(.data$label) |>
    dplyr::summarise(
      cell_n = dplyr::n(),
      latent_time_mean = mean(.data$latent01, na.rm = TRUE),
      terminal_probability = mean(.data$terminal_probability_cell, na.rm = TRUE),
      .groups = "drop"
    )
  grouped$latent_time_mean[!is.finite(grouped$latent_time_mean)] <- NA_real_
  grouped$terminal_probability[!is.finite(grouped$terminal_probability)] <- NA_real_
  min_cells_for_score <- cfg$velocity_root_terminal_min_cells_per_cluster
  grouped$initial_probability <- shrink_probability_10h(1 - grouped$latent_time_mean, grouped$cell_n, min_cells_for_score)
  grouped$terminal_probability <- shrink_probability_10h(grouped$terminal_probability, grouped$cell_n, min_cells_for_score)
  latent_terminal_score <- shrink_probability_10h(grouped$latent_time_mean, grouped$cell_n, min_cells_for_score)
  grouped$terminal_score <- ifelse(is.finite(grouped$terminal_probability), grouped$terminal_probability, latent_terminal_score)
  grouped$initial_score <- grouped$initial_probability
  grouped$root_score <- grouped$initial_score
  grouped$status <- "ok"
  grouped$reason <- ifelse(
    grouped$cell_n < min_cells_for_score,
    sprintf("probabilities shrunk toward 0.5 because cell_n < %s", min_cells_for_score),
    ""
  )
  out <- data.frame(
    pair_id = pair_id,
    split_value = split_value,
    cluster = grouped$label,
    cell_n = grouped$cell_n,
    initial_score = grouped$initial_score,
    root_score = grouped$root_score,
    terminal_score = grouped$terminal_score,
    initial_probability = grouped$initial_probability,
    terminal_probability = grouped$terminal_probability,
    latent_time_mean = grouped$latent_time_mean,
    status = grouped$status,
    reason = grouped$reason,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$initial_score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  out
}

scvelo_index <- velocity_read_scvelo_index_10(cfg)
cellrank_index <- velocity_read_cellrank_index_10(cfg, include_not_ok = TRUE)

index_rows <- list()
dynamic_outputs <- list()

for (i in seq_len(nrow(scvelo_index))) {
  unit <- scvelo_index[i, , drop = FALSE]
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  out_tsv <- velocity_root_terminal_path_10(cfg, pair_id, split_value)
  started <- proc.time()[["elapsed"]]
  result <- tryCatch({
    rt <- root_terminal_table_10h(unit, cellrank_index)
    write_tsv_local(rt, out_tsv)
    list(ok = TRUE, index = empty_index_row_10h(unit, out_tsv, "ok", "", round(proc.time()[["elapsed"]] - started, 3)), table = rt)
  }, error = function(e) {
    rt <- velocity_empty_df_10(velocity_root_terminal_cols_10)
    write_tsv_local(rt, out_tsv)
    list(ok = FALSE, index = empty_index_row_10h(unit, out_tsv, "failed", conditionMessage(e), round(proc.time()[["elapsed"]] - started, 3)), table = rt)
  })
  index_rows[[length(index_rows) + 1L]] <- result$index
  dynamic_outputs[[sprintf("velocity_root_terminal__%s", velocity_unit_file_id_10(pair_id, split_value))]] <- build_output_entry(
    out_tsv,
    "tsv",
    module_name,
    "velocity-derived root and terminal state probabilities by label",
    base_dir = cfg$project_root,
    schema = infer_schema_from_df(result$table)
  )
}

index_df <- if (length(index_rows) > 0) dplyr::bind_rows(index_rows) else velocity_empty_df_10(velocity_root_terminal_index_cols_10)
write_tsv_local(index_df, cfg$velocity_root_terminal_index_tsv)

velocity_manifest_10(
  cfg,
  cfg$module_10h_manifest_path,
  module_name,
  outputs = c(
    list(
      velocity_root_terminal_index_tsv = build_output_entry(
        cfg$velocity_root_terminal_index_tsv,
        "tsv",
        module_name,
        "velocity root/terminal status by unit",
        base_dir = cfg$project_root,
        schema = infer_schema_from_df(index_df)
      )
    ),
    dynamic_outputs
  ),
  inputs = list(
    module_10c = cfg$module_10c_manifest_path,
    module_10f = cfg$module_10f_manifest_path
  ),
  depends_on = list(
    module_10c = cfg$module_10c_manifest_path,
    module_10f = cfg$module_10f_manifest_path
  )
)

message("10h completed. velocity root/terminal index: ", cfg$velocity_root_terminal_index_tsv)
