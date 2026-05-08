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
source_utf8(file.path(.script_dir, "helpers", "plotting_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_09.R"))
source_utf8(file.path(.script_dir, "helpers", "trajectory_utils.R"))

load_required_packages(c("Seurat", "dplyr", "jsonlite", "ggplot2"))

cfg <- get_single_script_config_09()
module_name <- "09k_trajectory_figures"
prepare_dirs_09(cfg)

figure_index_cols_09k <- c(
  "pair_id", "split_value", "figure_role", "color_var", "png_path",
  "pdf_path", "n_cells", "status", "reason"
)

consensus_index_cols_09k <- c(
  trajectory_method_index_cols_09,
  "method_count", "methods_used", "correlation_path", "jaccard_path",
  "consensus_pseudotime_path", "conflict_count"
)

slingshot_index_cols_09k <- c(
  trajectory_method_index_cols_09,
  "fine_output_path", "lineage_count", "root_label", "label_var"
)

read_consensus_index_09k <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09i_manifest_path,
    "consensus_index_tsv",
    cfg$trajectory_consensus_index_tsv,
    cols = consensus_index_cols_09k
  )
}

read_slingshot_index_09k <- function(cfg) {
  trajectory_read_method_index_09(
    cfg$module_09d_manifest_path,
    "slingshot_index_tsv",
    cfg$trajectory_slingshot_index_tsv,
    cols = slingshot_index_cols_09k
  )
}

index_row_for_unit_09k <- function(index_df, pair_id, split_value) {
  if (nrow(index_df) == 0) {
    return(NULL)
  }
  split_value <- display_scalar_value(split_value, "pooled")
  hit <- index_df[index_df$pair_id == pair_id & index_df$split_value == split_value, , drop = FALSE]
  if (nrow(hit) == 0) NULL else hit[1, , drop = FALSE]
}

read_unit_pseudotime_09k <- function(consensus_index, slingshot_index, pair_id, split_value) {
  consensus <- index_row_for_unit_09k(consensus_index, pair_id, split_value)
  if (!is.null(consensus) &&
    consensus$status[[1]] == "ok" &&
    nzchar(consensus$consensus_pseudotime_path[[1]]) &&
    file.exists(consensus$consensus_pseudotime_path[[1]])) {
    df <- read.csv(consensus$consensus_pseudotime_path[[1]], stringsAsFactors = FALSE, check.names = FALSE)
    if ("cell_id" %in% colnames(df) && "consensus_pseudotime" %in% colnames(df)) {
      return(data.frame(
        cell_id = as.character(df$cell_id),
        pseudotime = trajectory_scale01_09(df$consensus_pseudotime),
        method = "consensus",
        stringsAsFactors = FALSE
      ))
    }
  }

  sling <- index_row_for_unit_09k(slingshot_index, pair_id, split_value)
  if (!is.null(sling) &&
    sling$status[[1]] == "ok" &&
    nzchar(sling$output_path[[1]]) &&
    file.exists(sling$output_path[[1]])) {
    df <- trajectory_read_pseudotime_csv_09(sling$output_path[[1]], "slingshot")
    if (nrow(df) > 0) {
      df$pseudotime <- trajectory_scale01_09(df$pseudotime)
      return(df[, c("cell_id", "pseudotime", "method"), drop = FALSE])
    }
  }

  data.frame(cell_id = character(0), pseudotime = numeric(0), method = character(0), stringsAsFactors = FALSE)
}

candidate_label_vars_09k <- function(unit, seu) {
  candidates <- unique(c(
    "cell_subtype",
    "seurat_clusters",
    normalize_scalar_value(unit$coarse_label_var[[1]]),
    normalize_scalar_value(unit$fine_label_var[[1]])
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates[candidates %in% colnames(seu@meta.data)]
}

figure_path_09k <- function(pair_id, split_value, stem, color_var = "") {
  out_dir <- trajectory_method_figure_dir_09(cfg, "figures", pair_id, split_value)
  ensure_dir(out_dir)
  value <- display_scalar_value(split_value, "pooled")
  name <- if (nzchar(color_var)) {
    sprintf("%s_%s__%s__%s.png", stem, safe_id_09(pair_id), safe_id_09(value), safe_id_09(color_var))
  } else {
    sprintf("%s_%s__%s.png", stem, safe_id_09(pair_id), safe_id_09(value))
  }
  file.path(out_dir, name)
}

empty_plot_09k <- function(title, label) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = label, size = 4) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title)
}

build_plot_df_09k <- function(seu, unit, pseudotime_df, color_var = "") {
  emb <- trajectory_umap_09(seu)
  if (is.null(emb)) {
    return(NULL)
  }
  cells <- rownames(emb)
  meta <- seu@meta.data[cells, , drop = FALSE]
  plot_df <- data.frame(
    cell_id = cells,
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    pseudotime = NA_real_,
    stringsAsFactors = FALSE
  )
  if (nrow(pseudotime_df) > 0) {
    plot_df$pseudotime <- suppressWarnings(as.numeric(pseudotime_df$pseudotime[match(cells, pseudotime_df$cell_id)]))
  }
  if (nzchar(color_var) && color_var %in% colnames(meta)) {
    plot_df$label <- as.character(meta[[color_var]])
    plot_df$label[is.na(plot_df$label) | !nzchar(plot_df$label)] <- "unknown"
  }
  plot_df
}

plot_label_umap_09k <- function(plot_df, color_var, title) {
  if (is.null(plot_df) || nrow(plot_df) == 0 || !"label" %in% colnames(plot_df)) {
    return(empty_plot_09k(title, sprintf("No UMAP or %s labels", color_var)))
  }
  plot_df$label <- factor(plot_df$label)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = label)) +
    ggplot2::geom_point(size = 0.25, alpha = 0.78, na.rm = TRUE) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2, alpha = 1))) +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2", color = color_var)
  labels <- levels(plot_df$label)
  pal <- paper_feature_palette()
  if (length(labels) > 0 && length(labels) <= length(pal)) {
    p <- p + ggplot2::scale_color_manual(values = setNames(pal[seq_along(labels)], labels), na.value = "grey80")
  }
  p
}

plot_root_terminal_09k <- function(seu, unit, pseudotime_df, title) {
  label_var <- normalize_scalar_value(unit$coarse_label_var[[1]])
  plot_df <- build_plot_df_09k(seu, unit, pseudotime_df, label_var)
  if (is.null(plot_df) || nrow(plot_df) == 0 || !"label" %in% colnames(plot_df)) {
    return(empty_plot_09k(title, "No root/terminal label overlay"))
  }

  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  root <- trajectory_selected_root_09(cfg, pair_id, split_value, unit$root_group[[1]])
  terminal <- normalize_scalar_value(unit$terminal_group[[1]])
  plot_df$state <- "other"
  if (nzchar(root) && !root %in% c("auto", "*")) {
    plot_df$state[plot_df$label == root] <- "root"
  }
  if (nzchar(terminal) && !terminal %in% c("auto", "*")) {
    plot_df$state[plot_df$label == terminal] <- "terminal"
  }
  highlight <- plot_df[plot_df$state != "other", , drop = FALSE]
  if (nrow(highlight) == 0) {
    return(empty_plot_09k(title, "Root/terminal labels were not found"))
  }

  ggplot2::ggplot(plot_df, ggplot2::aes(x = UMAP_1, y = UMAP_2)) +
    ggplot2::geom_point(color = "grey82", size = 0.22, alpha = 0.55, na.rm = TRUE) +
    ggplot2::geom_point(
      data = highlight,
      ggplot2::aes(color = state),
      size = 0.35,
      alpha = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c(root = "#D95F02", terminal = "#1B9E77")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = title, x = "UMAP 1", y = "UMAP 2", color = NULL)
}

save_figure_row_09k <- function(plot_obj, pair_id, split_value, role, color_var, png_path, n_cells, status = "ok", reason = "") {
  pdf_path <- sub("\\.png$", ".pdf", png_path)
  save_plot_dual(plot_obj, png_path, width = 6.2, height = 4.8, pdf_path = pdf_path)
  data.frame(
    pair_id = pair_id,
    split_value = display_scalar_value(split_value, "pooled"),
    figure_role = role,
    color_var = color_var,
    png_path = png_path,
    pdf_path = pdf_path,
    n_cells = n_cells,
    status = status,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

plot_unit_figures_09k <- function(unit, consensus_index, slingshot_index) {
  pair_id <- unit$pair_id[[1]]
  split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
  started <- proc.time()[["elapsed"]]
  tryCatch({
    seu <- readRDS(unit$input_rds[[1]])
    pseudotime_df <- read_unit_pseudotime_09k(consensus_index, slingshot_index, pair_id, split_value)
    label_vars <- candidate_label_vars_09k(unit, seu)
    rows <- list()
    for (color_var in label_vars) {
      plot_df <- build_plot_df_09k(seu, unit, pseudotime_df, color_var)
      plot_obj <- plot_label_umap_09k(plot_df, color_var, sprintf("Trajectory lineages %s %s by %s", pair_id, split_value, color_var))
      png_path <- figure_path_09k(pair_id, split_value, "Figure_5F_Lineages", color_var)
      rows[[length(rows) + 1L]] <- save_figure_row_09k(plot_obj, pair_id, split_value, "lineages", color_var, png_path, ncol(seu))
    }

    root_plot <- plot_root_terminal_09k(seu, unit, pseudotime_df, sprintf("Trajectory root/terminal %s %s", pair_id, split_value))
    root_png <- figure_path_09k(pair_id, split_value, "Figure_5F_RootTerminal")
    rows[[length(rows) + 1L]] <- save_figure_row_09k(root_plot, pair_id, split_value, "root_terminal_overlay", "", root_png, ncol(seu))
    dplyr::bind_rows(rows)
  }, error = function(e) {
    png_path <- figure_path_09k(pair_id, split_value, "Figure_5F_Lineages", "failed")
    plot_obj <- empty_plot_09k(sprintf("Trajectory lineages %s %s", pair_id, split_value), conditionMessage(e))
    save_figure_row_09k(
      plot_obj,
      pair_id,
      split_value,
      "lineages",
      "",
      png_path,
      suppressWarnings(as.integer(unit$cell_n_after[[1]])),
      status = "failed",
      reason = sprintf("%s; runtime_s=%.3f", conditionMessage(e), proc.time()[["elapsed"]] - started)
    )
  })
}

plot_split_compare_09k <- function(pair_units, consensus_index, slingshot_index) {
  pair_id <- pair_units$pair_id[[1]]
  split_values <- unique(vapply(pair_units$split_value, display_scalar_value, character(1), default = "pooled"))
  split_values <- setdiff(split_values, "pooled")
  if (length(split_values) < 2L) {
    return(trajectory_empty_df_09(figure_index_cols_09k))
  }

  color_vars <- unique(c("cell_subtype", "seurat_clusters"))
  unit_var_rows <- list()
  for (i in seq_len(nrow(pair_units))) {
    unit <- pair_units[i, , drop = FALSE]
    item <- tryCatch({
      seu <- readRDS(unit$input_rds[[1]])
      list(unit = unit, seu = seu, label_vars = candidate_label_vars_09k(unit, seu))
    }, error = function(e) NULL)
    if (!is.null(item)) {
      color_vars <- unique(c(color_vars, item$label_vars))
      unit_var_rows[[length(unit_var_rows) + 1L]] <- item
    }
  }
  color_vars <- color_vars[nzchar(color_vars)]
  rows <- list()
  for (color_var in color_vars) {
    combined <- list()
    for (item in unit_var_rows) {
      unit <- item$unit
      seu <- item$seu
      split_value <- display_scalar_value(unit$split_value[[1]], "pooled")
      if (!color_var %in% colnames(seu@meta.data)) {
        next
      }
      pseudotime_df <- read_unit_pseudotime_09k(consensus_index, slingshot_index, pair_id, split_value)
      plot_df <- build_plot_df_09k(seu, unit, pseudotime_df, color_var)
      if (is.null(plot_df) || !"label" %in% colnames(plot_df)) {
        next
      }
      plot_df$split_value <- split_value
      combined[[length(combined) + 1L]] <- plot_df
    }
    if (length(combined) == 0) {
      next
    }
    plot_df <- dplyr::bind_rows(combined)
    plot_df$split_value <- factor(plot_df$split_value, levels = split_values)
    plot_obj <- plot_label_umap_09k(plot_df, color_var, sprintf("Split trajectory comparison %s by %s", pair_id, color_var)) +
      ggplot2::facet_wrap(~split_value, nrow = 1)
    png_path <- figure_path_09k(pair_id, "split_compare", "Figure_5F_Lineages", color_var)
    rows[[length(rows) + 1L]] <- save_figure_row_09k(
      plot_obj,
      pair_id,
      "split_compare",
      "split_side_by_side",
      color_var,
      png_path,
      nrow(plot_df)
    )
  }
  if (length(rows) == 0) trajectory_empty_df_09(figure_index_cols_09k) else dplyr::bind_rows(rows)
}

units <- trajectory_execution_units_09(cfg)
consensus_index <- read_consensus_index_09k(cfg)
slingshot_index <- read_slingshot_index_09k(cfg)
figure_rows <- list()

for (i in seq_len(nrow(units))) {
  figure_rows[[length(figure_rows) + 1L]] <- plot_unit_figures_09k(units[i, , drop = FALSE], consensus_index, slingshot_index)
}

if (nrow(units) > 0) {
  for (pair_id in unique(units$pair_id)) {
    pair_units <- units[units$pair_id == pair_id, , drop = FALSE]
    figure_rows[[length(figure_rows) + 1L]] <- plot_split_compare_09k(pair_units, consensus_index, slingshot_index)
  }
}

figure_index <- if (length(figure_rows) > 0) dplyr::bind_rows(figure_rows) else trajectory_empty_df_09(figure_index_cols_09k)
write_tsv_local(figure_index, cfg$trajectory_figure_index_tsv)

dynamic_outputs <- list()
if (nrow(figure_index) > 0) {
  for (i in seq_len(nrow(figure_index))) {
    key <- sprintf(
      "trajectory_figure_%s_%s_%s",
      trajectory_unit_file_id_09(figure_index$pair_id[[i]], figure_index$split_value[[i]]),
      safe_id_09(figure_index$figure_role[[i]]),
      safe_id_09(display_scalar_value(figure_index$color_var[[i]], "overlay"))
    )
    dynamic_outputs[[key]] <- build_output_entry(
      figure_index$png_path[[i]],
      "png",
      module_name,
      "trajectory summary figure",
      base_dir = cfg$project_root
    )
  }
}

outputs <- c(
  list(
    figure_index_tsv = build_output_entry(
      cfg$trajectory_figure_index_tsv,
      "tsv",
      module_name,
      "trajectory figure index by pair/split",
      base_dir = cfg$project_root,
      schema = infer_schema_from_df(figure_index)
    )
  ),
  dynamic_outputs
)

trajectory_method_manifest_09(
  cfg,
  cfg$module_09k_manifest_path,
  module_name,
  outputs,
  inputs = list(
    module_09a = cfg$module_09a_manifest_path,
    module_09c = cfg$module_09c_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09i = cfg$module_09i_manifest_path
  ),
  depends_on = list(
    module_09a = cfg$module_09a_manifest_path,
    module_09c = cfg$module_09c_manifest_path,
    module_09d = cfg$module_09d_manifest_path,
    module_09i = cfg$module_09i_manifest_path
  )
)

message("09k completed. figure index: ", cfg$trajectory_figure_index_tsv)
