source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "annotation_stats_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

stage_dir <- ensure_eda_stage_dir(cfg, "pseudobulk_ds")
layer_status_df <- built_layer_status(cfg)
comparison_df <- read_comparison_sheet(cfg)

if (nrow(comparison_df) == 0) {
  stop("comparisons.tsv 中没有启用的比较设计，无法执行 pseudobulk DS。", call. = FALSE)
}

require_formal_pkg <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "formal pseudobulk DS 缺少依赖包: %s。请先更新 environment_r_main.yml 和 install_r_main_packages.R 对应环境。",
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
  if (!"annotation_label" %in% colnames(obj@meta.data)) {
    obj$annotation_label <- if ("cell_type" %in% colnames(obj@meta.data)) as.character(obj$cell_type) else obj$cluster_id
  }
  obj
}

run_formal_pbds <- function(seu, comparison_row, group_var, batch_var) {
  require_formal_pkg(c("SingleCellExperiment", "S4Vectors", "limma", "edgeR", "muscat"))

  keep <- seu@meta.data[[group_var]] %in% c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]])
  seu <- subset(seu, cells = colnames(seu)[keep])
  seu <- maybe_join_layers(seu)
  meta_df <- standardize_design_metadata_df(seu@meta.data)
  meta_df$group_id <- trim_character(meta_df[[group_var]], "")
  meta_df$sample_id <- trim_character(meta_df$sample_id, "")
  meta_df$cluster_id <- trim_character(seu$cluster_id, "")
  if (batch_var %in% colnames(meta_df)) {
    meta_df$batch <- trim_character(meta_df[[batch_var]], "default")
  } else {
    meta_df$batch <- "default"
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = get_assay_matrix(seu, assay = "RNA", type = "counts")),
    colData = S4Vectors::DataFrame(meta_df)
  )
  sce <- muscat::prepSCE(sce, kid = "cluster_id", sid = "sample_id", gid = "group_id", drop = FALSE)
  pb <- muscat::aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))

  sample_info <- unique(meta_df[, c("sample_id", "group_id", "batch"), drop = FALSE])
  sample_info <- sample_info[match(colnames(pb), sample_info$sample_id), , drop = FALSE]
  sample_info$group_id <- factor(sample_info$group_id, levels = c(comparison_row$ident_1[[1]], comparison_row$ident_2[[1]]))
  rownames(sample_info) <- sample_info$sample_id

  design_formula <- if (length(unique(sample_info$batch)) >= 2) "~ 0 + group_id + batch" else "~ 0 + group_id"
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_info)
  rownames(design) <- rownames(sample_info)
  contrast <- limma::makeContrasts(
    contrasts = sprintf("group_id%s-group_id%s", make.names(comparison_row$ident_1[[1]]), make.names(comparison_row$ident_2[[1]])),
    levels = design
  )

  res <- muscat::pbDS(
    pb = pb,
    method = "limma-trend",
    design = design,
    contrast = contrast,
    min_cells = 10,
    filter = "both",
    verbose = FALSE
  )
  tidy_res <- muscat::resDS(sce, res, bind = "row")
  tidy_res$comparison_id <- comparison_row$comparison_id[[1]]
  tidy_res
}

report_lines <- c(
  "# Pseudobulk DS Report",
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

  layer_paths <- pseudobulk_paths(cfg, layer_id)
  ensure_parent_dirs(layer_paths)

  aggregation <- aggregate_cluster_sample_counts(obj, cluster_var = "cluster_id")
  saveRDS(aggregation$counts, layer_paths$aggregation_rds)
  write_tsv(aggregation$metadata, layer_paths$aggregation_tsv)
  write_tsv(
    data.frame(
      layer_id = layer_id,
      aggregation_rds = layer_paths$aggregation_rds,
      aggregation_metadata_tsv = layer_paths$aggregation_tsv,
      pseudobulk_columns = ncol(aggregation$counts),
      stringsAsFactors = FALSE
    ),
    layer_paths$aggregation_manifest_tsv
  )

  report_lines <- c(
    report_lines,
    sprintf("## Layer `%s`", layer_id),
    sprintf("- aggregation counts: `%s`", layer_paths$aggregation_rds),
    sprintf("- aggregation metadata: `%s`", layer_paths$aggregation_tsv)
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

    cmp_paths <- pseudobulk_paths(cfg, layer_id, comparison_row$comparison_id[[1]])
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
          exploratory_marker_discovery = marker_discovery_paths(cfg, layer_id, comparison_row$comparison_id[[1]])$exploratory_tsv,
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", comparison_row$comparison_id[[1]]),
        warning_banner,
        sprintf("- status: exploratory_only"),
        sprintf("- reason: %s", gate$reason),
        sprintf("- gate summary: `%s`", cmp_paths$gate_summary_tsv),
        sprintf("- exploratory fallback: `%s`", marker_discovery_paths(cfg, layer_id, comparison_row$comparison_id[[1]])$exploratory_tsv)
      )
      manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
        layer_id = layer_id,
        comparison_id = comparison_row$comparison_id[[1]],
        inference_status = "exploratory_only",
        aggregation_rds = layer_paths$aggregation_rds,
        ds_results_tsv = "",
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv,
        stringsAsFactors = FALSE
      )
      next
    }

    formal_res <- run_formal_pbds(obj, comparison_row, group_var = group_var, batch_var = batch_var)
    write_tsv(formal_res, cmp_paths$ds_results_tsv)
    saveRDS(formal_res, cmp_paths$result_rds)
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
      sprintf("- formal DS results: `%s`", cmp_paths$ds_results_tsv)
    )
    manifest_rows[[length(manifest_rows) + 1]] <- data.frame(
      layer_id = layer_id,
      comparison_id = comparison_row$comparison_id[[1]],
      inference_status = "formal",
      aggregation_rds = layer_paths$aggregation_rds,
      ds_results_tsv = cmp_paths$ds_results_tsv,
      gate_summary_tsv = cmp_paths$gate_summary_tsv,
      status_tsv = cmp_paths$status_tsv,
      stringsAsFactors = FALSE
    )
  }
}

manifest_path <- file.path(cfg$table_dir, "pseudobulk_ds", "pseudobulk_manifest.tsv")
write_tsv(bind_rows(manifest_rows), manifest_path)
report_lines <- c(report_lines, "", sprintf("- manifest: `%s`", manifest_path))
write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("pseudobulk DS 完成，manifest: ", manifest_path)
