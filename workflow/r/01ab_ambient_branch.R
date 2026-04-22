source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))
source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "plotting_helpers.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "00_raw_objects.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 00_raw_objects.rds", call. = FALSE)
}

stage_dir <- ensure_eda_stage_dir(cfg, "ambient")
ambient_ckpt_dir <- file.path(cfg$checkpoint_dir, "ambient")
dir.create(ambient_ckpt_dir, recursive = TRUE, showWarnings = FALSE)

raw_obj <- readRDS(input_rds)
raw_obj <- maybe_join_layers(raw_obj)
inventory_df <- read_input_inventory(cfg)
readiness_df <- read_branch_readiness(cfg)
split_objs <- split_object_by_sample(raw_obj)

have_soupx <- requireNamespace("SoupX", quietly = TRUE)
have_decontx <- requireNamespace("celda", quietly = TRUE) &&
  requireNamespace("SingleCellExperiment", quietly = TRUE)

empty_df <- function() data.frame(stringsAsFactors = FALSE)

resolve_sample_barcodes <- function(sample_obj) {
  if ("barcode_raw" %in% colnames(sample_obj@meta.data)) {
    barcodes <- as.character(sample_obj$barcode_raw)
  } else {
    barcodes <- colnames(sample_obj)
  }
  barcodes[!nzchar(barcodes)] <- colnames(sample_obj)[!nzchar(barcodes)]
  barcodes
}

align_count_matrices <- function(filtered_counts, raw_counts) {
  common_features <- intersect(rownames(filtered_counts), rownames(raw_counts))
  if (length(common_features) == 0) {
    stop("filtered/raw 矩阵没有共同 feature。", call. = FALSE)
  }
  filtered_counts <- filtered_counts[common_features, , drop = FALSE]
  raw_counts <- raw_counts[common_features, , drop = FALSE]
  list(filtered = filtered_counts, raw = raw_counts)
}

extract_soupx_contamination <- function(sc) {
  candidate_values <- c()
  candidate_values <- c(candidate_values, tryCatch(sc$metaData$rho, error = function(e) numeric(0)))
  candidate_values <- c(candidate_values, tryCatch(sc$fit$rhoEst, error = function(e) numeric(0)))
  candidate_values <- c(candidate_values, tryCatch(sc$fit$rho, error = function(e) numeric(0)))
  candidate_values <- suppressWarnings(as.numeric(candidate_values))
  candidate_values <- candidate_values[is.finite(candidate_values)]
  if (length(candidate_values) == 0) {
    return(NA_real_)
  }
  stats::median(candidate_values)
}

should_apply_ambient <- function(policy, recommended_flag) {
  policy <- tolower(normalize_scalar_value(policy, "manual"))
  if (policy %in% c("auto", "always")) {
    return(TRUE)
  }
  if (policy %in% c("apply_recommended", "recommended")) {
    return(isTRUE(recommended_flag))
  }
  FALSE
}

run_soupx_for_sample <- function(sample_obj, raw_counts, branch_model) {
  filtered_counts <- get_assay_matrix(sample_obj, assay = "RNA", type = "counts")
  sample_barcodes <- resolve_sample_barcodes(sample_obj)
  colnames(filtered_counts) <- sample_barcodes
  aligned <- align_count_matrices(filtered_counts, raw_counts)
  filtered_counts <- aligned$filtered
  raw_counts <- aligned$raw

  cluster_labels <- as.character(branch_model$object$provisional_cluster)
  names(cluster_labels) <- sample_barcodes
  cluster_labels <- cluster_labels[colnames(filtered_counts)]

  sc <- SoupX::SoupChannel(tod = raw_counts, toc = filtered_counts)
  sc <- SoupX::setClusters(sc, clusters = cluster_labels)
  sc <- SoupX::autoEstCont(sc, doPlot = FALSE)
  corrected_counts <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  colnames(corrected_counts) <- colnames(sample_obj)

  list(
    corrected_counts = corrected_counts,
    contamination_fraction = extract_soupx_contamination(sc),
    per_cell_contamination = rep(extract_soupx_contamination(sc), ncol(sample_obj)),
    method = "soupx",
    cluster_labels = as.character(branch_model$object$provisional_cluster)
  )
}

run_decontx_for_sample <- function(sample_obj, raw_counts, branch_model) {
  filtered_counts <- get_assay_matrix(sample_obj, assay = "RNA", type = "counts")
  sample_barcodes <- resolve_sample_barcodes(sample_obj)
  colnames(filtered_counts) <- sample_barcodes
  cluster_labels <- as.character(branch_model$object$provisional_cluster)
  names(cluster_labels) <- sample_barcodes

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = filtered_counts)
  )
  if (!is.null(raw_counts)) {
    aligned <- align_count_matrices(filtered_counts, raw_counts)
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = aligned$filtered)
    )
    cluster_labels <- cluster_labels[colnames(aligned$filtered)]
    bg_sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = aligned$raw)
    )
    sce <- celda::decontX(sce, z = cluster_labels, background = bg_sce)
  } else {
    sce <- celda::decontX(sce, z = cluster_labels)
  }

  corrected_counts <- tryCatch(
    celda::decontXcounts(sce),
    error = function(e) SummarizedExperiment::assay(sce, "decontXcounts")
  )
  corrected_counts <- as.matrix(corrected_counts)
  colnames(corrected_counts) <- colnames(sample_obj)[match(colnames(corrected_counts), sample_barcodes)]
  corrected_counts <- corrected_counts[, colnames(sample_obj), drop = FALSE]

  per_cell <- tryCatch(
    as.numeric(SummarizedExperiment::colData(sce)$decontX_contamination),
    error = function(e) rep(NA_real_, ncol(sample_obj))
  )
  if (length(per_cell) != ncol(sample_obj)) {
    per_cell <- rep(NA_real_, ncol(sample_obj))
  }

  list(
    corrected_counts = corrected_counts,
    contamination_fraction = if (all(is.na(per_cell))) NA_real_ else stats::median(per_cell, na.rm = TRUE),
    per_cell_contamination = per_cell,
    method = "decontx",
    cluster_labels = as.character(branch_model$object$provisional_cluster)
  )
}

sample_rows <- list()
decision_rows <- list()
leakage_rows <- list()
cellbender_rows <- list()
sample_objects <- list()

for (sample_id in names(split_objs)) {
  message("运行 ambient branch: ", sample_id)
  sample_dir <- file.path(ambient_ckpt_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  sample_obj <- maybe_join_layers(split_objs[[sample_id]])
  inventory_row <- resolve_inventory_row(cfg, sample_id, inventory_df = inventory_df)
  readiness_row <- if (nrow(readiness_df) > 0 && "sample_id" %in% colnames(readiness_df)) {
    readiness_df[readiness_df$sample_id == sample_id, , drop = FALSE]
  } else {
    NULL
  }
  method_plan <- choose_ambient_method(readiness_row, cfg)
  branch_model <- build_sample_branch_model(
    sample_obj = sample_obj,
    dims = cfg$ambient_cluster_dims,
    resolution = cfg$ambient_cluster_resolution,
    hvg_nfeatures = cfg$hvg_nfeatures,
    min_cells = cfg$ambient_min_cells
  )

  before_counts <- get_assay_matrix(sample_obj, assay = "RNA", type = "counts")
  corrected_counts <- before_counts
  contamination_fraction <- NA_real_
  per_cell_contamination <- rep(NA_real_, ncol(sample_obj))
  applied_method <- "none"
  ambient_status <- branch_model$status
  execution_note <- normalize_scalar_value(method_plan$ambient_notes)
  raw_counts <- NULL
  raw_read_error <- ""

  raw_path <- if (!is.null(inventory_row) && "raw_matrix_dir" %in% colnames(inventory_row)) {
    normalize_scalar_value(inventory_row$raw_matrix_dir[1])
  } else {
    ""
  }
  if (nzchar(raw_path) && file.exists(raw_path)) {
    raw_counts <- tryCatch(
      read_matrix_from_path(raw_path),
      error = function(e) {
        raw_read_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  cellbender_row <- build_cellbender_stub(sample_id, sample_obj, inventory_row, cfg, raw_counts = raw_counts)
  cellbender_row$cellbender_ready <- method_plan$cellbender_ready
  cellbender_row$cellbender_reason <- if (method_plan$cellbender_ready) "stub_only" else normalize_scalar_value(method_plan$ambient_notes, "not_ready")
  cellbender_rows[[sample_id]] <- cellbender_row

  run_result <- NULL
  if (identical(branch_model$status, "built")) {
    if (identical(method_plan$preferred_method, "soupx")) {
      if (!have_soupx) {
        ambient_status <- "skipped_soupx_package_missing"
      } else if (is.null(raw_counts)) {
        ambient_status <- "skipped_soupx_raw_unavailable"
        execution_note <- collapse_unique_values(c(execution_note, raw_read_error), sep = "; ")
      } else {
        run_result <- tryCatch(
          run_soupx_for_sample(sample_obj, raw_counts, branch_model),
          error = function(e) {
            execution_note <<- collapse_unique_values(c(execution_note, conditionMessage(e)), sep = "; ")
            NULL
          }
        )
        if (!is.null(run_result)) {
          ambient_status <- "completed_soupx"
        }
      }
    }

    if (is.null(run_result) && identical(method_plan$decontx_ready, TRUE)) {
      if (!have_decontx) {
        ambient_status <- if (grepl("^completed", ambient_status)) ambient_status else "skipped_decontx_package_missing"
      } else if (identical(method_plan$preferred_method, "soupx")) {
        run_result <- tryCatch(
          run_decontx_for_sample(sample_obj, raw_counts, branch_model),
          error = function(e) {
            execution_note <<- collapse_unique_values(c(execution_note, conditionMessage(e)), sep = "; ")
            NULL
          }
        )
        ambient_status <- if (!is.null(run_result)) "fallback_decontx" else "failed_soupx_and_decontx"
      } else if (identical(method_plan$preferred_method, "decontx")) {
        run_result <- tryCatch(
          run_decontx_for_sample(sample_obj, raw_counts, branch_model),
          error = function(e) {
            execution_note <<- collapse_unique_values(c(execution_note, conditionMessage(e)), sep = "; ")
            NULL
          }
        )
        if (!is.null(run_result)) {
          ambient_status <- "completed_decontx"
        } else {
          ambient_status <- "failed_decontx"
        }
      }
    }
  }

  if (!is.null(run_result)) {
    corrected_counts <- run_result$corrected_counts
    contamination_fraction <- run_result$contamination_fraction
    per_cell_contamination <- run_result$per_cell_contamination
    applied_method <- run_result$method

    corrected_path <- file.path(sample_dir, sprintf("%s_%s_corrected_counts.rds", sample_id, applied_method))
    saveRDS(corrected_counts, corrected_path)

    marker_df <- safe_find_all_markers(branch_model$object, top_n = cfg$ambient_marker_top_n)
    leakage_df <- calculate_marker_leakage(
      before_counts = before_counts,
      after_counts = corrected_counts,
      cluster_labels = run_result$cluster_labels,
      marker_df = marker_df
    )
    if (nrow(leakage_df) > 0) {
      leakage_df$sample_id <- sample_id
      leakage_df$method <- applied_method
      leakage_rows[[sample_id]] <- leakage_df
      write_tsv(leakage_df, file.path(sample_dir, sprintf("%s_marker_leakage.tsv", sample_id)))
    }
  }

  recommended_to_replace <- ambient_recommend_apply(contamination_fraction, cfg)
  applied_to_main <- FALSE
  if (!is.null(run_result) && should_apply_ambient(cfg$ambient_apply_policy, recommended_to_replace)) {
    sample_obj[["RNA"]] <- CreateAssayObject(counts = corrected_counts)
    sample_obj <- add_basic_qc_metrics(sample_obj, cfg, declared_gene_id_type = normalize_scalar_value(sample_obj$gene_id_type[1], "auto"))$object
    applied_to_main <- TRUE
  }

  sample_obj$ambient_requested_method <- method_plan$preferred_method
  sample_obj$ambient_applied_method <- applied_method
  sample_obj$ambient_status <- ambient_status
  sample_obj$ambient_recommended_to_replace <- recommended_to_replace
  sample_obj$ambient_applied_to_main <- applied_to_main
  sample_obj$ambient_contamination_fraction <- contamination_fraction
  sample_obj$ambient_cell_contamination <- per_cell_contamination
  sample_objects[[sample_id]] <- sample_obj

  sample_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    cells_input = ncol(sample_obj),
    raw_matrix_available = method_plan$raw_matrix_available,
    raw_matrix_kind = method_plan$raw_matrix_kind,
    preferred_method = method_plan$preferred_method,
    fallback_method = method_plan$fallback_method,
    executed_method = applied_method,
    ambient_status = ambient_status,
    contamination_fraction = contamination_fraction,
    recommended_to_replace = recommended_to_replace,
    applied_to_main = applied_to_main,
    branch_model_status = branch_model$status,
    soupx_package_available = have_soupx,
    decontx_package_available = have_decontx,
    execution_note = execution_note,
    stringsAsFactors = FALSE
  )

  decision_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    ambient_any_ready = method_plan$decontx_ready,
    ambient_soupx_ready = method_plan$soupx_ready,
    ambient_decontx_ready = method_plan$decontx_ready,
    ambient_cellbender_ready = method_plan$cellbender_ready,
    preferred_method = method_plan$preferred_method,
    fallback_method = method_plan$fallback_method,
    apply_policy = cfg$ambient_apply_policy,
    recommended_to_replace = recommended_to_replace,
    applied_to_main = applied_to_main,
    ambient_status = ambient_status,
    stringsAsFactors = FALSE
  )
}

final_obj <- merge_named_objects(sample_objects)
readiness_export <- readiness_df
if (nrow(readiness_export) > 0 && "sample_id" %in% colnames(readiness_export)) {
  readiness_export <- readiness_export[readiness_export$sample_id != "__PROJECT__", , drop = FALSE]
}
sample_summary_df <- bind_rows(sample_rows)
decision_df <- bind_rows(decision_rows)
leakage_df <- if (length(leakage_rows) > 0) bind_rows(leakage_rows) else empty_df()
cellbender_df <- bind_rows(cellbender_rows)

saveRDS(final_obj, file.path(cfg$checkpoint_dir, "01_after_ambient_branch.rds"))
write_tsv(readiness_export, file.path(stage_dir, "ambient_readiness.tsv"))
write_tsv(sample_summary_df, file.path(stage_dir, "ambient_sample_summary.tsv"))
write_tsv(decision_df, file.path(stage_dir, "ambient_method_decisions.tsv"))
write_tsv(leakage_df, file.path(stage_dir, "ambient_marker_leakage.tsv"))
write_tsv(cellbender_df, file.path(stage_dir, "cellbender_stub.tsv"))

contam_plot_df <- sample_summary_df %>%
  mutate(sample_id = factor(sample_id, levels = sample_id))
if (sum(is.finite(contam_plot_df$contamination_fraction)) > 0) {
  contam_plot <- ggplot(contam_plot_df, aes(x = sample_id, y = contamination_fraction, fill = executed_method)) +
    geom_col() +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Ambient contamination summary", x = NULL, y = "Contamination fraction", fill = "Method")
  save_plot_dual(contam_plot, file.path(stage_dir, "ambient_contamination_summary.png"), width = 9, height = 5)
}

if (nrow(leakage_df) > 0) {
  leakage_plot_df <- leakage_df %>%
    group_by(sample_id) %>%
    summarise(
      leakage_before = mean(leakage_before, na.rm = TRUE),
      leakage_after = mean(leakage_after, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c("leakage_before", "leakage_after"), names_to = "stage", values_to = "value")
  leakage_plot <- ggplot(leakage_plot_df, aes(x = sample_id, y = value, fill = stage)) +
    geom_col(position = "dodge") +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Marker leakage before/after ambient correction", x = NULL, y = "Mean off/on cluster leakage")
  save_plot_dual(leakage_plot, file.path(stage_dir, "ambient_marker_leakage.png"), width = 10, height = 5)
}

report_lines <- c(
  "# Ambient Branch Report",
  "",
  sprintf("- 样本数: `%s`", nrow(sample_summary_df)),
  sprintf("- SoupX package available: `%s`", have_soupx),
  sprintf("- DecontX package available: `%s`", have_decontx),
  sprintf("- apply_policy: `%s`", cfg$ambient_apply_policy),
  "",
  "## Key Files",
  sprintf("- `ambient_readiness.tsv`: `%s`", file.path(stage_dir, "ambient_readiness.tsv")),
  sprintf("- `ambient_sample_summary.tsv`: `%s`", file.path(stage_dir, "ambient_sample_summary.tsv")),
  sprintf("- `ambient_method_decisions.tsv`: `%s`", file.path(stage_dir, "ambient_method_decisions.tsv")),
  sprintf("- `ambient_marker_leakage.tsv`: `%s`", file.path(stage_dir, "ambient_marker_leakage.tsv")),
  sprintf("- `cellbender_stub.tsv`: `%s`", file.path(stage_dir, "cellbender_stub.tsv")),
  "",
  "## Review Focus",
  "- 明确每个样本是否有 raw droplets、实际运行了哪种 ambient 方法、以及是否建议替换主流程 counts。",
  "- `CellBender` 本轮仅输出 readiness/config/command stub，不阻塞主线。"
)

if (nrow(sample_summary_df) > 0) {
  report_lines <- c(report_lines, "", "## Sample Decisions")
  for (i in seq_len(nrow(sample_summary_df))) {
    row <- sample_summary_df[i, , drop = FALSE]
    report_lines <- c(
      report_lines,
      sprintf(
        "- `%s`: preferred=`%s`; executed=`%s`; status=`%s`; raw_available=`%s`; contamination_fraction=`%s`; recommended_to_replace=`%s`; applied_to_main=`%s`",
        row$sample_id,
        row$preferred_method,
        row$executed_method,
        row$ambient_status,
        row$raw_matrix_available,
        ifelse(is.finite(row$contamination_fraction), sprintf("%.4f", row$contamination_fraction), "NA"),
        row$recommended_to_replace,
        row$applied_to_main
      )
    )
    if (nzchar(normalize_scalar_value(row$execution_note))) {
      report_lines <- c(report_lines, sprintf("  note=`%s`", row$execution_note))
    }
  }
}

write_markdown(report_lines, file.path(stage_dir, "report.md"))
message("ambient branch 已输出到: ", stage_dir)
