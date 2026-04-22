source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "annotation_stats_helpers.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

stage_dir <- ensure_eda_stage_dir(cfg, "composition")
layer_status_df <- built_layer_status(cfg)
comparison_df <- read_comparison_sheet(cfg)

if (nrow(comparison_df) == 0) {
  stop("comparisons.tsv 中没有启用的比较设计，无法执行 composition analysis。", call. = FALSE)
}

require_formal_pkg <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "formal composition analysis 缺少依赖包: %s。请先更新 environment_r_main.yml 和 install_r_main_packages.R 对应环境。",
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

load_layer_object <- function(layer_row) {
  ann_paths <- layer_annotation_paths(cfg, layer_row$layer_id[[1]])
  obj_path <- if (file.exists(ann_paths$annotated_rds)) ann_paths$annotated_rds else layer_row$clustered_rds[[1]]
  obj <- readRDS(obj_path)
  obj <- maybe_join_layers(obj)
  obj <- standardize_design_metadata(obj)
  cluster_var <- if ("cluster_id" %in% colnames(obj@meta.data)) {
    "cluster_id"
  } else if (nzchar(layer_row$cluster_column[[1]]) && layer_row$cluster_column[[1]] %in% colnames(obj@meta.data)) {
    layer_row$cluster_column[[1]]
  } else {
    "seurat_clusters"
  }
  obj$cluster_id <- as.character(obj@meta.data[[cluster_var]])
  obj
}

plot_proportion_summary <- function(summary_df, output_png) {
  if (nrow(summary_df) == 0) {
    return(invisible(NULL))
  }
  plot_obj <- ggplot(summary_df, aes(x = sample_id, y = proportion, fill = cluster_id)) +
    geom_col() +
    theme_classic(base_size = 10) +
    labs(title = "Sample-level cluster proportions", x = NULL, y = "Proportion", fill = "Cluster")
  save_plot_dual(plot_obj, output_png, width = 8, height = 5)
}

run_formal_propeller <- function(seu, comparison_row, group_var, batch_var) {
  require_formal_pkg(c("speckle", "limma"))

  keep <- seu@meta.data[[group_var]] %in% c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]])
  meta_df <- standardize_design_metadata_df(seu@meta.data[keep, , drop = FALSE])
  meta_df$cluster_id <- as.character(seu$cluster_id[keep])
  meta_df$group_id <- trim_character(meta_df[[group_var]], "")
  if (batch_var %in% colnames(meta_df)) {
    meta_df$batch <- trim_character(meta_df[[batch_var]], "default")
  } else {
    meta_df$batch <- "default"
  }

  props <- speckle::getTransformedProps(
    clusters = meta_df$cluster_id,
    sample = meta_df$sample_id,
    transform = "logit"
  )
  sample_info <- unique(meta_df[, c("sample_id", "group_id", "batch"), drop = FALSE])
  sample_info <- sample_info[match(colnames(props$Proportions), sample_info$sample_id), , drop = FALSE]
  sample_info$group_id <- factor(sample_info$group_id, levels = c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]]))

  design_formula <- if (length(unique(sample_info$batch)) >= 2) "~ 0 + group_id + batch" else "~ 0 + group_id"
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_info)
  contrast <- limma::makeContrasts(
    contrasts = sprintf("group_id%s-group_id%s", make.names(comparison_row$ident_1[[1]]), make.names(comparison_row$ident_2[[1]])),
    levels = design
  )

  res <- speckle::propeller.ttest(
    prop.list = props,
    design = design,
    contrasts = contrast,
    robust = TRUE,
    trend = FALSE,
    sort = TRUE
  )
  res <- tibble::rownames_to_column(as.data.frame(res), "cluster_id")
  res$comparison_id <- comparison_row$comparison_id[[1]]
  res
}

report_lines <- c(
  "# Composition Analysis Report",
  "",
  sprintf("- comparison sheet: `%s`", cfg$comparison_sheet),
  sprintf("- min biological replicates per group: `%s`", cfg$min_biological_replicates),
  ""
)

manifest_rows <- list()

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  obj <- load_layer_object(layer_row)

  layer_paths <- composition_paths(cfg, layer_id)
  ensure_parent_dirs(layer_paths)

  proportion_df <- sample_level_proportion_summary(obj, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id")
  write_tsv(proportion_df, layer_paths$proportion_tsv)
  tryCatch({
    plot_proportion_summary(proportion_df, layer_paths$proportion_plot_png)
  }, error = function(e) {
    warning(sprintf("layer `%s` composition plot 失败: %s", layer_id, conditionMessage(e)), call. = FALSE)
  })

  report_lines <- c(
    report_lines,
    sprintf("## Layer `%s`", layer_id),
    sprintf("- proportion summary: `%s`", layer_paths$proportion_tsv),
    sprintf("- proportion plot: `%s`", layer_paths$proportion_plot_png)
  )

  for (cmp_idx in seq_len(nrow(comparison_df))) {
    comparison_row <- comparison_df[cmp_idx, , drop = FALSE]
    if (!comparison_applies_to_layer(comparison_row, layer_id)) {
      next
    }

    group_var <- normalize_scalar_value(comparison_row$group_var[[1]], "group_id")
    if (!group_var %in% colnames(obj@meta.data)) {
      group_var <- "group_id"
    }
    batch_var <- normalize_scalar_value(comparison_row$batch_var[[1]], "batch")
    gate <- replicate_gate_summary(
      meta_df = obj@meta.data,
      ident_1 = comparison_row$ident_1[[1]],
      ident_2 = comparison_row$ident_2[[1]],
      group_var = group_var,
      replicate_var = "biological_replicate",
      sample_var = "sample_id",
      min_reps = comparison_row$min_biological_replicates[[1]]
    )

    cmp_paths <- composition_paths(cfg, layer_id, comparison_row$comparison_id[[1]])
    ensure_parent_dirs(cmp_paths)
    write_tsv(gate$summary, cmp_paths$gate_summary_tsv)

    if (!gate$pass) {
      warning_banner <- annotation_warning_banner(gate$summary, min_reps = comparison_row$min_biological_replicates[[1]])
      write_tsv(
        data.frame(
          comparison_id = comparison_row$comparison_id[[1]],
          layer_id = layer_id,
          inference_status = "exploratory_only",
          reason = gate$reason,
          warning_banner = warning_banner,
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", comparison_row$comparison_id[[1]]),
        warning_banner,
        "- status: exploratory_only",
        sprintf("- reason: %s", gate$reason),
        sprintf("- gate summary: `%s`", cmp_paths$gate_summary_tsv)
      )
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = comparison_row$comparison_id[[1]],
        inference_status = "exploratory_only",
        proportion_tsv = layer_paths$proportion_tsv,
        formal_results_tsv = "",
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv,
        stringsAsFactors = FALSE
      )
      next
    }

    formal_res <- run_formal_propeller(obj, comparison_row, group_var = group_var, batch_var = batch_var)
    write_tsv(formal_res, cmp_paths$formal_results_tsv)
    write_tsv(
      data.frame(
        comparison_id = comparison_row$comparison_id[[1]],
        layer_id = layer_id,
        inference_status = "formal",
        reason = gate$reason,
        warning_banner = "",
        stringsAsFactors = FALSE
      ),
      cmp_paths$status_tsv
    )
    report_lines <- c(
      report_lines,
      sprintf("### `%s`", comparison_row$comparison_id[[1]]),
      "- status: formal",
      sprintf("- gate summary: `%s`", cmp_paths$gate_summary_tsv),
      sprintf("- formal propeller results: `%s`", cmp_paths$formal_results_tsv)
    )
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      comparison_id = comparison_row$comparison_id[[1]],
      inference_status = "formal",
      proportion_tsv = layer_paths$proportion_tsv,
      formal_results_tsv = cmp_paths$formal_results_tsv,
      gate_summary_tsv = cmp_paths$gate_summary_tsv,
      status_tsv = cmp_paths$status_tsv,
      stringsAsFactors = FALSE
    )
  }
}

manifest_path <- file.path(cfg$table_dir, "composition", "composition_manifest.tsv")
write_tsv(bind_rows(manifest_rows), manifest_path)
report_lines <- c(report_lines, "", sprintf("- manifest: `%s`", manifest_path))
write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("composition analysis 完成，manifest: ", manifest_path)
