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
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "layer_config_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "comparison_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "deg_utils.R"))

load_required_packages(c("Seurat", "dplyr", "ggplot2", "tibble", "jsonlite"))

cfg <- get_single_script_config_05()
module_name <- "05c_composition"
prepare_dirs_05(cfg)
set.seed(cfg$random_seed)

require_formal_pkg_05c <- function(pkgs, label) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(sprintf("%s missing required R packages: %s", label, paste(missing, collapse = ", ")), call. = FALSE)
  }
}

plot_proportion_summary_05 <- function(summary_df, output_png) {
  ensure_dir(dirname(output_png))
  if (nrow(summary_df) == 0) {
    plot_obj <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No composition data")
  } else {
    plot_obj <- ggplot2::ggplot(summary_df, ggplot2::aes(x = sample_id, y = proportion, fill = cluster_id)) +
      ggplot2::geom_col() +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::labs(title = "Sample-level cluster proportions", x = NULL, y = "Proportion", fill = "Cluster")
  }
  save_plot_local(plot_obj, output_png, width = 8, height = 5)
}

run_formal_propeller_05 <- function(seu, vars) {
  require_formal_pkg_05c(c("speckle", "limma"), "formal composition")

  meta_df <- standardize_design_metadata_df_05(seu@meta.data)
  meta_df$deg_composition_group <- as.character(seu$deg_composition_group)
  meta_df$group_id <- trimws(as.character(meta_df[[vars$group_var]]))
  meta_df$batch <- if (vars$batch_var %in% colnames(meta_df)) trimws(as.character(meta_df[[vars$batch_var]])) else "default"
  meta_df$batch[!nzchar(meta_df$batch)] <- "default"

  props <- speckle::getTransformedProps(
    clusters = meta_df$deg_composition_group,
    sample = meta_df$sample_id,
    transform = "logit"
  )
  sample_info <- unique(meta_df[, c("sample_id", "group_id", "batch"), drop = FALSE])
  sample_info <- sample_info[match(colnames(props$Proportions), sample_info$sample_id), , drop = FALSE]
  sample_info$group_id <- factor(sample_info$group_id, levels = c(vars$ident_1, vars$ident_2))

  design_formula <- if (length(unique(sample_info$batch)) >= 2) "~ 0 + group_id + batch" else "~ 0 + group_id"
  design <- stats::model.matrix(stats::as.formula(design_formula), data = sample_info)
  contrast <- limma::makeContrasts(
    contrasts = sprintf("group_id%s-group_id%s", make.names(vars$ident_1), make.names(vars$ident_2)),
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
  res$comparison_id <- vars$comparison_id
  res
}

empty_composition_table_05c <- function() {
  empty_df_05(c("sample_id", "group_id", "cluster_id", "cell_number", "sample_total", "proportion"))
}

normalize_output_path_05c <- function(path) {
  path <- normalize_scalar_value(path)
  if (!nzchar(path)) "" else normalizePath(path, winslash = "/", mustWork = FALSE)
}

qc_composition_notice_05c <- "该结果只表示 scRNA 捕获到的 GC/TC 细胞组成。不能解释为组织中 GC/TC 的真实比例变化。真实组织比例由 ST region 注释或 ST 去卷积评估。"

write_qc_composition_outputs_05c <- function(cfg, vars, layer_id, cmp_paths, composition_df,
                                             inference_status, reason, formal_results_tsv = "") {
  if (!identical(vars$analysis_mode, "qc_composition")) {
    return(list(
      composition_tsv = normalize_output_path_05c(cmp_paths$proportion_tsv),
      qc_mirror_tsv = "",
      qc_report_md = "",
      qc_mirror_written = "no",
      biological_conclusion_allowed = "yes",
      tissue_abundance_interpretation_allowed = "yes"
    ))
  }

  if (is.null(composition_df) || !is.data.frame(composition_df)) {
    composition_df <- empty_composition_table_05c()
  }
  qc_paths <- qc_composition_paths_05(cfg, vars$output_alias)
  write_tsv_local(composition_df, cmp_paths$proportion_tsv)
  write_tsv_local(composition_df, cmp_paths$comparison_tsv)
  write_tsv_local(composition_df, qc_paths$table_tsv)

  sample_n <- if ("sample_id" %in% colnames(composition_df)) length(unique(composition_df$sample_id[nzchar(as.character(composition_df$sample_id))])) else 0L
  group_n <- if ("group_id" %in% colnames(composition_df)) length(unique(composition_df$group_id[nzchar(as.character(composition_df$group_id))])) else 0L
  report_lines <- c(
    sprintf("# %s", vars$report_title),
    "",
    sprintf("- canonical_question_id: `%s`", vars$source_question_id),
    sprintf("- display_question_id: `%s`", vars$display_question_id),
    sprintf("- comparison_id: `%s`", vars$comparison_id),
    sprintf("- layer_id: `%s`", layer_id),
    "- analysis_mode: `qc_composition`",
    "- gene_program_role: `qc_only`",
    "- biological_conclusion_allowed: `no`",
    "- tissue_abundance_interpretation_allowed: `no`",
    sprintf("- inference_status: `%s`", inference_status),
    sprintf("- reason: %s", reason),
    sprintf("- canonical_table: `%s`", cmp_paths$comparison_tsv),
    sprintf("- qc_mirror_table: `%s`", qc_paths$table_tsv),
    sprintf("- formal_results_tsv: `%s`", formal_results_tsv),
    sprintf("- sample_n: `%s`", sample_n),
    sprintf("- group_n: `%s`", group_n),
    "",
    "## Interpretation Guardrail",
    "",
    qc_composition_notice_05c,
    "",
    "E03 is a scRNA captured-cell balance QC result.",
    "It is not a spatial or tissue abundance estimate.",
    "Do not use E03 to claim GC/TC tissue proportion changes.",
    "Use ST region annotation, ST deconvolution, or histology/image quantification for tissue abundance."
  )
  ensure_dir(dirname(qc_paths$report_md))
  write_markdown_local(report_lines, qc_paths$report_md)

  list(
    composition_tsv = normalize_output_path_05c(cmp_paths$comparison_tsv),
    qc_mirror_tsv = normalize_output_path_05c(qc_paths$table_tsv),
    qc_report_md = normalize_output_path_05c(qc_paths$report_md),
    qc_mirror_written = "yes",
    biological_conclusion_allowed = "no",
    tissue_abundance_interpretation_allowed = "no"
  )
}

composition_manifest_row_05c <- function(layer_id, vars, inference_status,
                                         proportion_tsv = "", formal_results_tsv = "",
                                         gate_summary_tsv = "", status_tsv = "",
                                         qc_info = NULL) {
  if (is.null(qc_info)) {
    qc_info <- list(
      composition_tsv = normalize_output_path_05c(proportion_tsv),
      qc_mirror_tsv = "",
      qc_report_md = "",
      qc_mirror_written = "no",
      biological_conclusion_allowed = if (identical(vars$analysis_mode, "qc_composition")) "no" else "yes",
      tissue_abundance_interpretation_allowed = if (identical(vars$analysis_mode, "qc_composition")) "no" else "yes"
    )
  }
  data.frame(
    layer_id = layer_id,
    comparison_id = vars$comparison_id,
    source_question_id = vars$source_question_id,
    display_question_id = vars$display_question_id,
    output_alias = vars$output_alias,
    report_title = vars$report_title,
    analysis_mode = vars$analysis_mode,
    gene_program_role = vars$gene_program_role,
    is_qc_composition = if (identical(vars$analysis_mode, "qc_composition")) "yes" else "no",
    composition_group_var = vars$composition_group_var,
    inference_status = inference_status,
    proportion_tsv = normalize_output_path_05c(proportion_tsv),
    composition_tsv = normalize_scalar_value(qc_info$composition_tsv),
    formal_results_tsv = normalize_output_path_05c(formal_results_tsv),
    gate_summary_tsv = normalize_output_path_05c(gate_summary_tsv),
    status_tsv = normalize_output_path_05c(status_tsv),
    qc_mirror_tsv = normalize_scalar_value(qc_info$qc_mirror_tsv),
    qc_report_md = normalize_scalar_value(qc_info$qc_report_md),
    biological_conclusion_allowed = normalize_scalar_value(qc_info$biological_conclusion_allowed),
    tissue_abundance_interpretation_allowed = normalize_scalar_value(qc_info$tissue_abundance_interpretation_allowed),
    qc_mirror_written = normalize_scalar_value(qc_info$qc_mirror_written, "no"),
    stringsAsFactors = FALSE
  )
}

qc_report_summary_lines_05c <- function(vars, qc_info) {
  if (!identical(vars$analysis_mode, "qc_composition")) {
    return(character(0))
  }
  c(
    "- qc_composition: `yes`",
    "- biological_conclusion_allowed: `no`",
    "- tissue_abundance_interpretation_allowed: `no`",
    sprintf("- qc_mirror_tsv: `%s`", normalize_scalar_value(qc_info$qc_mirror_tsv)),
    sprintf("- qc_report_md: `%s`", normalize_scalar_value(qc_info$qc_report_md)),
    sprintf("- interpretation_guardrail: %s", qc_composition_notice_05c)
  )
}

layer_status_df <- deg_layer_status(cfg)
comparison_df <- read_deg_comparison_sheet(cfg)

report_lines <- c(
  "# 05c Composition",
  "",
  sprintf("- comparison_sheet: `%s`", cfg$comparison_sheet),
  sprintf("- min_biological_replicates_default: `%s`", cfg$min_biological_replicates),
  ""
)
manifest_rows <- list()

for (idx in seq_len(nrow(layer_status_df))) {
  layer_row <- layer_status_df[idx, , drop = FALSE]
  layer_id <- layer_row$layer_id[[1]]
  message("05c composition layer: ", layer_id)
  obj <- load_layer_for_deg(cfg, layer_row)

  layer_paths <- composition_paths_05(cfg, layer_id)
  proportion_df <- sample_level_proportion_summary_05(obj, cluster_var = "cluster_id", sample_var = "sample_id", group_var = "group_id")
  write_tsv_local(proportion_df, layer_paths$proportion_tsv)
  tryCatch(plot_proportion_summary_05(proportion_df, layer_paths$proportion_plot_png), error = function(e) warning(conditionMessage(e), call. = FALSE))

  report_lines <- c(
    report_lines,
    sprintf("## Layer `%s`", layer_id),
    sprintf("- sample_level_proportions: `%s`", layer_paths$proportion_tsv),
    sprintf("- proportion_plot: `%s`", layer_paths$proportion_plot_png)
  )

  for (cmp_idx in seq_len(nrow(comparison_df))) {
    comparison_row <- comparison_df[cmp_idx, , drop = FALSE]
    if (!comparison_applies_to_layer(comparison_row, layer_id)) {
      next
    }

    vars <- resolve_comparison_vars_05(obj, comparison_row)
    cmp_paths <- composition_paths_05(cfg, layer_id, vars$comparison_id)
    subset_result <- subset_cells_for_comparison(obj, vars)

    if (!vars$analysis_mode %in% c("composition", "qc_composition")) {
      write_tsv_local(empty_gate_summary_05(), cmp_paths$gate_summary_tsv)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          analysis_mode = vars$analysis_mode,
          composition_group_var = vars$composition_group_var,
          inference_status = "skipped_non_composition",
          reason = "05c composition is only run for composition or qc_composition rows",
          warning_banner = "",
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      manifest_rows[[length(manifest_rows) + 1]] <- composition_manifest_row_05c(
        layer_id, vars, "skipped_non_composition",
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv
      )
      next
    }

    if (!identical(subset_result$status, "ok")) {
      write_tsv_local(empty_gate_summary_05(), cmp_paths$gate_summary_tsv)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          analysis_mode = vars$analysis_mode,
          composition_group_var = vars$composition_group_var,
          inference_status = "skipped",
          reason = subset_result$reason,
          warning_banner = "",
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      qc_info <- write_qc_composition_outputs_05c(
        cfg, vars, layer_id, cmp_paths, empty_composition_table_05c(),
        "skipped", subset_result$reason
      )
      manifest_rows[[length(manifest_rows) + 1]] <- composition_manifest_row_05c(
        layer_id, vars, "skipped",
        proportion_tsv = if (identical(vars$analysis_mode, "qc_composition")) cmp_paths$proportion_tsv else "",
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv,
        qc_info = qc_info
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", vars$comparison_id),
        sprintf("- status: skipped; reason: %s", subset_result$reason),
        qc_report_summary_lines_05c(vars, qc_info)
      )
      next
    }

    composition_result <- resolve_composition_group_var_05(subset_result$object, vars)
    if (!identical(composition_result$status, "ok")) {
      write_tsv_local(empty_gate_summary_05(), cmp_paths$gate_summary_tsv)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          analysis_mode = vars$analysis_mode,
          composition_group_var = vars$composition_group_var,
          inference_status = "skipped",
          reason = composition_result$reason,
          warning_banner = "",
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      qc_info <- write_qc_composition_outputs_05c(
        cfg, vars, layer_id, cmp_paths, empty_composition_table_05c(),
        "skipped", composition_result$reason
      )
      manifest_rows[[length(manifest_rows) + 1]] <- composition_manifest_row_05c(
        layer_id, vars, "skipped",
        proportion_tsv = if (identical(vars$analysis_mode, "qc_composition")) cmp_paths$proportion_tsv else "",
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv,
        qc_info = qc_info
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", vars$comparison_id),
        sprintf("- status: skipped; reason: %s", composition_result$reason),
        qc_report_summary_lines_05c(vars, qc_info)
      )
      next
    }

    obj_sub <- composition_result$object
    obj_sub$deg_composition_group <- as.character(obj_sub@meta.data[[composition_result$group_var]])
    cmp_proportion_df <- sample_level_proportion_summary_05(
      obj_sub,
      cluster_var = "deg_composition_group",
      sample_var = "sample_id",
      group_var = vars$group_var
    )
    write_tsv_local(cmp_proportion_df, cmp_paths$proportion_tsv)

    gate <- replicate_gate_summary_05(obj_sub@meta.data, vars)
    write_tsv_local(gate$summary, cmp_paths$gate_summary_tsv)

    if (!gate$pass) {
      status <- if (vars$force_exploratory) "exploratory_forced" else "exploratory_only"
      warning_banner <- annotation_warning_banner_05(gate$summary, vars, forced = vars$force_exploratory)
      write_tsv_local(
        data.frame(
          comparison_id = vars$comparison_id,
          layer_id = layer_id,
          analysis_mode = vars$analysis_mode,
          composition_group_var = vars$composition_group_var,
          inference_status = status,
          reason = gate$reason,
          warning_banner = warning_banner,
          stringsAsFactors = FALSE
        ),
        cmp_paths$status_tsv
      )
      qc_info <- write_qc_composition_outputs_05c(
        cfg, vars, layer_id, cmp_paths, cmp_proportion_df,
        status, gate$reason
      )
      manifest_rows[[length(manifest_rows) + 1]] <- composition_manifest_row_05c(
        layer_id, vars, status,
        proportion_tsv = cmp_paths$proportion_tsv,
        gate_summary_tsv = cmp_paths$gate_summary_tsv,
        status_tsv = cmp_paths$status_tsv,
        qc_info = qc_info
      )
      report_lines <- c(
        report_lines,
        sprintf("### `%s`", vars$comparison_id),
        warning_banner,
        sprintf("- status: %s", status),
        sprintf("- reason: %s", gate$reason),
        qc_report_summary_lines_05c(vars, qc_info)
      )
      next
    }

    formal_res <- run_formal_propeller_05(obj_sub, vars)
    write_tsv_local(formal_res, cmp_paths$formal_results_tsv)
    write_tsv_local(
      data.frame(
        comparison_id = vars$comparison_id,
        layer_id = layer_id,
        analysis_mode = vars$analysis_mode,
        composition_group_var = vars$composition_group_var,
        inference_status = "formal",
        reason = gate$reason,
        warning_banner = "",
        stringsAsFactors = FALSE
      ),
      cmp_paths$status_tsv
    )
    qc_info <- write_qc_composition_outputs_05c(
      cfg, vars, layer_id, cmp_paths, cmp_proportion_df,
      "formal", gate$reason, formal_results_tsv = cmp_paths$formal_results_tsv
    )
    manifest_rows[[length(manifest_rows) + 1]] <- composition_manifest_row_05c(
      layer_id, vars, "formal",
      proportion_tsv = cmp_paths$proportion_tsv,
      formal_results_tsv = cmp_paths$formal_results_tsv,
      gate_summary_tsv = cmp_paths$gate_summary_tsv,
      status_tsv = cmp_paths$status_tsv,
      qc_info = qc_info
    )
    report_lines <- c(
      report_lines,
      sprintf("### `%s`", vars$comparison_id),
      "- status: formal",
      sprintf("- formal_propeller: `%s`", cmp_paths$formal_results_tsv),
      qc_report_summary_lines_05c(vars, qc_info)
    )
  }
}

manifest_df <- if (length(manifest_rows) > 0) dplyr::bind_rows(manifest_rows) else empty_df_05(c(
  "layer_id", "comparison_id", "source_question_id", "display_question_id",
  "output_alias", "report_title", "analysis_mode", "gene_program_role",
  "is_qc_composition", "composition_group_var", "inference_status",
  "proportion_tsv", "composition_tsv", "formal_results_tsv",
  "gate_summary_tsv", "status_tsv", "qc_mirror_tsv", "qc_report_md",
  "biological_conclusion_allowed", "tissue_abundance_interpretation_allowed",
  "qc_mirror_written"
))
manifest_tsv <- file.path(cfg$composition_table_dir, "composition_manifest.tsv")
write_tsv_local(manifest_df, manifest_tsv)
report_lines <- c(report_lines, "", sprintf("- manifest: `%s`", manifest_tsv))
report_path <- file.path(cfg$composition_report_dir, "report.md")
write_markdown_local(report_lines, report_path)

if (file.exists(cfg$module_05c_manifest_path)) {
  unlink(cfg$module_05c_manifest_path)
}
write_manifest_local(
  manifest_path = cfg$module_05c_manifest_path,
  new_outputs = list(
    composition_manifest_tsv = build_output_entry(manifest_tsv, "tsv", module_name, "one row per layer/comparison composition status", base_dir = cfg$project_root, schema = infer_schema_from_df(manifest_df)),
    report = build_output_entry(report_path, "md", module_name, "composition report", base_dir = cfg$project_root)
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(layer_status_tsv = cfg$layer_status_file, comparison_sheet = cfg$comparison_sheet),
  version = cfg$module_version,
  depends_on = list(module_05b = cfg$module_05b_manifest_path)
)

message("05c completed. manifest: ", manifest_tsv)
