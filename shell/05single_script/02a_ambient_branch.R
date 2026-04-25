#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  }
)

source_utf8 <- function(path) {
  source(path, encoding = "UTF-8")
}

source_utf8(file.path(.script_dir, "helpers", "runtime_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "config.R"))
source_utf8(file.path(.script_dir, "helpers", "project_paths_02.R"))
source_utf8(file.path(.script_dir, "helpers", "gtf_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "plotting_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tidyr", "tibble", "scales", "jsonlite"))

cfg <- get_single_script_config_02()
module_name <- "02a_ambient_branch"
prepare_dirs_02(cfg)
set.seed(cfg$random_seed)

manifest_01a <- read_manifest_local(cfg$module_01a_manifest_path)
input_rds <- resolve_output_local(manifest_01a, "raw_object")
if (!file.exists(input_rds)) {
  stop(sprintf("缺少输入检查点: %s", input_rds), call. = FALSE)
}

stage_dir <- cfg$ambient_report_dir
ambient_ckpt_dir <- file.path(cfg$checkpoint_dir, "ambient")
ensure_dir(ambient_ckpt_dir)

raw_obj <- readRDS(input_rds)
raw_obj <- maybe_join_layers(raw_obj)
inventory_df <- read_input_inventory_local(cfg)
branch_readiness_df <- read_branch_readiness_local(cfg)
require_input_inventory_contract_local(cfg, inventory_df)
require_branch_readiness_contract_local(cfg, branch_readiness_df)

split_objs <- split_object_by_sample(raw_obj)
have_soupx <- requireNamespace("SoupX", quietly = TRUE)
have_decontx <- requireNamespace("celda", quietly = TRUE) &&
  requireNamespace("SingleCellExperiment", quietly = TRUE) &&
  requireNamespace("SummarizedExperiment", quietly = TRUE)

sample_rows <- list()
decision_rows <- list()
leakage_rows <- list()
cellbender_rows <- list()
sample_objects <- list()

for (sample_id in names(split_objs)) {
  message("运行 ambient branch: ", sample_id)
  sample_dir <- file.path(ambient_ckpt_dir, sample_id)
  ensure_dir(sample_dir)

  sample_obj <- maybe_join_layers(split_objs[[sample_id]])
  inventory_row <- resolve_inventory_row_local(cfg, sample_id, inventory_df = inventory_df)
  readiness_row <- if (nrow(branch_readiness_df) > 0 && "sample_id" %in% colnames(branch_readiness_df)) {
    branch_readiness_df[branch_readiness_df$sample_id == sample_id, , drop = FALSE]
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
  execution_notes <- c(normalize_scalar_value(method_plan$ambient_notes))
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
    if (nzchar(raw_read_error)) {
      execution_notes <- c(execution_notes, raw_read_error)
    }
  }

  cellbender_row <- build_cellbender_stub(sample_id, sample_obj, inventory_row, cfg, raw_counts = raw_counts)
  cellbender_row$cellbender_ready <- method_plan$cellbender_ready
  cellbender_row$cellbender_reason <- if (method_plan$cellbender_ready) "stub_only" else normalize_scalar_value(method_plan$ambient_notes, "not_ready")
  cellbender_rows[[sample_id]] <- cellbender_row

  run_result <- NULL
  if (identical(branch_model$status, "built")) {
    method_sequence <- unique(tolower(c(method_plan$preferred_method, method_plan$fallback_method)))
    method_sequence <- method_sequence[nzchar(method_sequence) & !method_sequence %in% c("none", "na")]
    if (length(method_sequence) == 0) {
      ambient_status <- "skipped_no_method"
    }
    for (method_name in method_sequence) {
      if (!ambient_method_ready(method_name, method_plan)) {
        execution_notes <- c(execution_notes, sprintf("%s_not_ready", method_name))
        next
      }
      dispatch <- dispatch_ambient_method(method_name, sample_obj, raw_counts, branch_model)
      if (nzchar(normalize_scalar_value(dispatch$note))) {
        execution_notes <- c(execution_notes, dispatch$note)
      }
      ambient_status <- dispatch$status
      if (!is.null(dispatch$result)) {
        run_result <- dispatch$result
        applied_method <- dispatch$result$method
        if (!identical(method_name, method_plan$preferred_method)) {
          ambient_status <- paste0("fallback_", applied_method)
        }
        break
      }
    }
    if (is.null(run_result) && length(method_sequence) > 0 && !grepl("^skipped|^failed|^fallback|^completed", ambient_status)) {
      ambient_status <- "skipped_no_available_method"
    }
  }

  if (!is.null(run_result)) {
    corrected_counts <- run_result$corrected_counts
    contamination_fraction <- run_result$contamination_fraction
    per_cell_contamination <- run_result$per_cell_contamination

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
      write_tsv_local(leakage_df, file.path(sample_dir, sprintf("%s_marker_leakage.tsv", sample_id)))
    }
  }

  recommended_to_replace <- ambient_recommend_apply(contamination_fraction, cfg)
  applied_to_main <- FALSE
  if (!is.null(run_result) && should_apply_ambient(method_plan$apply_policy, recommended_to_replace)) {
    sample_obj[["RNA"]] <- Seurat::CreateAssayObject(counts = corrected_counts)
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
    apply_policy = method_plan$apply_policy,
    branch_model_status = branch_model$status,
    soupx_package_available = have_soupx,
    decontx_package_available = have_decontx,
    execution_note = collapse_unique_values(execution_notes, sep = "; "),
    stringsAsFactors = FALSE
  )

  decision_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    ambient_any_ready = isTRUE(method_plan$soupx_ready) || isTRUE(method_plan$decontx_ready),
    ambient_soupx_ready = method_plan$soupx_ready,
    ambient_decontx_ready = method_plan$decontx_ready,
    ambient_cellbender_ready = method_plan$cellbender_ready,
    preferred_method = method_plan$preferred_method,
    fallback_method = method_plan$fallback_method,
    apply_policy = method_plan$apply_policy,
    recommended_to_replace = recommended_to_replace,
    applied_to_main = applied_to_main,
    ambient_status = ambient_status,
    stringsAsFactors = FALSE
  )
}

final_obj <- merge_named_objects(sample_objects)
readiness_export <- branch_readiness_df
if (nrow(readiness_export) > 0 && "sample_id" %in% colnames(readiness_export)) {
  readiness_export <- readiness_export[readiness_export$sample_id != "__PROJECT__", , drop = FALSE]
}
sample_summary_df <- dplyr::bind_rows(sample_rows)
decision_df <- dplyr::bind_rows(decision_rows)
leakage_df <- if (length(leakage_rows) > 0) dplyr::bind_rows(leakage_rows) else data.frame(
  gene = character(0),
  source_cluster = character(0),
  mean_in_source_before = numeric(0),
  mean_outside_before = numeric(0),
  leakage_before = numeric(0),
  mean_in_source_after = numeric(0),
  mean_outside_after = numeric(0),
  leakage_after = numeric(0),
  sample_id = character(0),
  method = character(0),
  stringsAsFactors = FALSE
)
cellbender_df <- dplyr::bind_rows(cellbender_rows)

ambient_object_path <- file.path(cfg$checkpoint_dir, "02a_after_ambient_branch.rds")
ambient_readiness_path <- file.path(stage_dir, "ambient_readiness.tsv")
ambient_summary_path <- file.path(stage_dir, "ambient_sample_summary.tsv")
ambient_decisions_path <- file.path(stage_dir, "ambient_method_decisions.tsv")
ambient_leakage_path <- file.path(stage_dir, "ambient_marker_leakage.tsv")
cellbender_path <- file.path(stage_dir, "cellbender_stub.tsv")
contam_png <- file.path(stage_dir, "ambient_contamination_summary.png")
leakage_png <- file.path(stage_dir, "ambient_marker_leakage.png")
report_path <- file.path(stage_dir, "report.md")

saveRDS(final_obj, ambient_object_path)
write_tsv_local(readiness_export, ambient_readiness_path)
write_tsv_local(sample_summary_df, ambient_summary_path)
write_tsv_local(decision_df, ambient_decisions_path)
write_tsv_local(leakage_df, ambient_leakage_path)
write_tsv_local(cellbender_df, cellbender_path)

contam_plot_df <- sample_summary_df %>%
  dplyr::mutate(
    sample_id = factor(sample_id, levels = sample_id),
    contamination_fraction_plot = ifelse(is.finite(contamination_fraction), contamination_fraction, 0)
  )
contam_plot <- ggplot2::ggplot(contam_plot_df, ggplot2::aes(x = sample_id, y = contamination_fraction_plot, fill = executed_method)) +
  ggplot2::geom_col() +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(title = "Ambient contamination summary", x = NULL, y = "Contamination fraction", fill = "Method")
save_plot_dual(contam_plot, contam_png, width = 9, height = 5)

if (nrow(leakage_df) > 0) {
  leakage_plot_df <- leakage_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      leakage_before = mean(leakage_before, na.rm = TRUE),
      leakage_after = mean(leakage_after, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = c("leakage_before", "leakage_after"), names_to = "stage", values_to = "value")
  leakage_plot <- ggplot2::ggplot(leakage_plot_df, ggplot2::aes(x = sample_id, y = value, fill = stage)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = "Marker leakage before/after ambient correction", x = NULL, y = "Mean off/on cluster leakage")
} else {
  leakage_plot <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No marker leakage table available") +
    ggplot2::theme_void()
}
save_plot_dual(leakage_plot, leakage_png, width = 10, height = 5)

report_lines <- c(
  "# Ambient Branch Report",
  "",
  sprintf("- 样本数: `%s`", nrow(sample_summary_df)),
  sprintf("- SoupX package available: `%s`", have_soupx),
  sprintf("- DecontX package available: `%s`", have_decontx),
  sprintf("- apply_policy: `%s`", cfg$ambient_apply_policy),
  "",
  "## Key Files",
  sprintf("- `ambient_object`: `%s`", ambient_object_path),
  sprintf("- `ambient_readiness.tsv`: `%s`", ambient_readiness_path),
  sprintf("- `ambient_sample_summary.tsv`: `%s`", ambient_summary_path),
  sprintf("- `ambient_method_decisions.tsv`: `%s`", ambient_decisions_path),
  sprintf("- `ambient_marker_leakage.tsv`: `%s`", ambient_leakage_path),
  sprintf("- `cellbender_stub.tsv`: `%s`", cellbender_path),
  "",
  "## Sample Decisions"
)

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

report_lines <- c(
  report_lines,
  "",
  "## Review Focus",
  "- 明确每个样本是否有 raw droplets、实际运行了哪种 ambient 方法、以及是否建议替换主流程 counts。",
  "- `CellBender` 本轮仅输出 readiness/config/command stub，不阻塞主线。"
)

ensure_dir(dirname(report_path))
write_markdown_local(report_lines, report_path)

write_manifest_local(
  manifest_path = cfg$module_02a_manifest_path,
  new_outputs = list(
    ambient_object = build_output_entry(ambient_object_path, "rds", module_name, "merged Seurat object after ambient branch metadata/correction policy", base_dir = cfg$project_root),
    ambient_readiness = build_output_entry(ambient_readiness_path, "tsv", module_name, "one row per sample ambient readiness input", base_dir = cfg$project_root, schema = infer_schema_from_df(readiness_export)),
    ambient_sample_summary = build_output_entry(ambient_summary_path, "tsv", module_name, "one row per sample ambient execution summary", base_dir = cfg$project_root, schema = infer_schema_from_df(sample_summary_df)),
    ambient_method_decisions = build_output_entry(ambient_decisions_path, "tsv", module_name, "one row per sample ambient method decision", base_dir = cfg$project_root, schema = infer_schema_from_df(decision_df)),
    ambient_marker_leakage = build_output_entry(ambient_leakage_path, "tsv", module_name, "one row per marker/cluster/sample leakage comparison", base_dir = cfg$project_root, schema = infer_schema_from_df(leakage_df)),
    cellbender_stub = build_output_entry(cellbender_path, "tsv", module_name, "one row per sample CellBender command stub", base_dir = cfg$project_root, schema = infer_schema_from_df(cellbender_df)),
    ambient_contamination_summary_png = build_output_entry(contam_png, "png", module_name, "ambient contamination summary plot", base_dir = cfg$project_root),
    ambient_contamination_summary_pdf = build_output_entry(sub("\\.png$", ".pdf", contam_png), "pdf", module_name, "ambient contamination summary plot", base_dir = cfg$project_root),
    ambient_marker_leakage_png = build_output_entry(leakage_png, "png", module_name, "ambient marker leakage plot", base_dir = cfg$project_root),
    ambient_marker_leakage_pdf = build_output_entry(sub("\\.png$", ".pdf", leakage_png), "pdf", module_name, "ambient marker leakage plot", base_dir = cfg$project_root),
    report = build_output_entry(report_path, "md", module_name, "human-readable ambient branch report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    manifest_01a = cfg$module_01a_manifest_path,
    raw_object = input_rds,
    input_inventory_file = cfg$input_inventory_file,
    branch_readiness_file = cfg$branch_readiness_file
  ),
  version = cfg$module_version,
  depends_on = list(
    module_01a = cfg$module_01a_manifest_path
  )
)

message("02a 完成。ambient branch 已输出到: ", stage_dir)
