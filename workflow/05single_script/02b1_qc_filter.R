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
source_utf8(file.path(.script_dir, "helpers", "ambient_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_02()
module_name <- "02b1_qc_filter"
prepare_dirs_02(cfg)
set.seed(cfg$random_seed)

manifest_02a <- read_manifest_local(cfg$module_02a_manifest_path)
ambient_rds <- resolve_output_local(manifest_02a, "ambient_object")
if (!file.exists(ambient_rds)) {
  stop(sprintf("缺少 ambient 输入对象: %s", ambient_rds), call. = FALSE)
}

ambient_obj <- readRDS(ambient_rds)
ambient_obj <- maybe_join_layers(ambient_obj)
threshold_overrides <- read_qc_threshold_overrides_local(cfg)
split_objs <- split_object_by_sample(ambient_obj)

qc_fail_reasons <- function(nfeature_fail, ncount_fail, log10umi_fail, mito_fail) {
  reasons <- character(length(nfeature_fail))
  for (i in seq_along(reasons)) {
    parts <- c(
      if (isTRUE(nfeature_fail[[i]])) "low_nfeature" else "",
      if (isTRUE(ncount_fail[[i]])) "low_ncount" else "",
      if (isTRUE(log10umi_fail[[i]])) "low_log10umi" else "",
      if (isTRUE(mito_fail[[i]])) "high_mito" else ""
    )
    parts <- parts[nzchar(parts)]
    reasons[[i]] <- if (length(parts) == 0) "" else paste(parts, collapse = ";")
  }
  reasons
}

filtered_list <- list()
summary_rows <- list()
diagnostic_rows <- list()

for (sample_id in names(split_objs)) {
  message("运行基础 QC filter: ", sample_id)
  sample_obj <- maybe_join_layers(split_objs[[sample_id]])
  thresholds <- get_sample_qc_thresholds_local(cfg, sample_id, threshold_overrides)

  if (!"percent.mito" %in% colnames(sample_obj@meta.data)) {
    sample_obj$percent.mito <- 0
  }
  if (!"percent.ribo" %in% colnames(sample_obj@meta.data)) {
    sample_obj$percent.ribo <- 0
  }
  if (!"percent.rbc" %in% colnames(sample_obj@meta.data)) {
    sample_obj$percent.rbc <- 0
  }
  if (!"log10GenesPerUMI" %in% colnames(sample_obj@meta.data)) {
    sample_obj$log10GenesPerUMI <- log10(sample_obj$nFeature_RNA + 1) / log10(sample_obj$nCount_RNA + 1)
    sample_obj$log10GenesPerUMI[!is.finite(sample_obj$log10GenesPerUMI)] <- 0
  }

  meta_df <- sample_obj@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      sample_id = sample_id,
      barcode = if ("barcode_raw" %in% colnames(.)) barcode_raw else cell_id,
      fail_low_nfeature = nFeature_RNA < thresholds$qc_min_nfeature,
      fail_low_ncount = nCount_RNA < thresholds$qc_min_ncount,
      fail_low_log10umi = log10GenesPerUMI < thresholds$qc_min_log10umi,
      fail_high_mito = percent.mito > thresholds$qc_max_mito_pct,
      qc_pass = !(fail_low_nfeature | fail_low_ncount | fail_low_log10umi | fail_high_mito),
      qc_fail_reason = qc_fail_reasons(fail_low_nfeature, fail_low_ncount, fail_low_log10umi, fail_high_mito),
      qc_min_nfeature = thresholds$qc_min_nfeature,
      qc_min_ncount = thresholds$qc_min_ncount,
      qc_min_log10umi = thresholds$qc_min_log10umi,
      qc_max_mito_pct = thresholds$qc_max_mito_pct
    )

  rownames(meta_df) <- meta_df$cell_id
  sample_obj@meta.data[colnames(sample_obj), "qc_pass"] <- meta_df[colnames(sample_obj), "qc_pass"]
  sample_obj@meta.data[colnames(sample_obj), "qc_fail_reason"] <- meta_df[colnames(sample_obj), "qc_fail_reason"]
  sample_obj@meta.data[colnames(sample_obj), "qc_min_nfeature"] <- thresholds$qc_min_nfeature
  sample_obj@meta.data[colnames(sample_obj), "qc_min_ncount"] <- thresholds$qc_min_ncount
  sample_obj@meta.data[colnames(sample_obj), "qc_min_log10umi"] <- thresholds$qc_min_log10umi
  sample_obj@meta.data[colnames(sample_obj), "qc_max_mito_pct"] <- thresholds$qc_max_mito_pct
  filtered_list[[sample_id]] <- sample_obj

  summary_rows[[sample_id]] <- data.frame(
    sample_id = sample_id,
    cells_in = nrow(meta_df),
    cells_pass = sum(meta_df$qc_pass, na.rm = TRUE),
    cells_fail = sum(!meta_df$qc_pass, na.rm = TRUE),
    fail_low_nfeature = sum(meta_df$fail_low_nfeature, na.rm = TRUE),
    fail_low_ncount = sum(meta_df$fail_low_ncount, na.rm = TRUE),
    fail_low_log10umi = sum(meta_df$fail_low_log10umi, na.rm = TRUE),
    fail_high_mito = sum(meta_df$fail_high_mito, na.rm = TRUE),
    qc_min_nfeature = thresholds$qc_min_nfeature,
    qc_min_ncount = thresholds$qc_min_ncount,
    qc_min_log10umi = thresholds$qc_min_log10umi,
    qc_max_mito_pct = thresholds$qc_max_mito_pct,
    stringsAsFactors = FALSE
  )

  diagnostic_rows[[sample_id]] <- meta_df %>%
    dplyr::select(
      sample_id,
      barcode,
      cell_id,
      qc_pass,
      qc_fail_reason,
      nCount_RNA,
      nFeature_RNA,
      percent.mito,
      percent.ribo,
      percent.rbc,
      log10GenesPerUMI,
      qc_min_nfeature,
      qc_min_ncount,
      qc_min_log10umi,
      qc_max_mito_pct
    )
}

qc_filter_object_path <- file.path(cfg$checkpoint_dir, "02b1_after_qc_filter.rds")
qc_filter_summary_path <- file.path(cfg$table_dir, "qc_filter_sample_summary.tsv")
qc_filter_diagnostics_path <- file.path(cfg$table_dir, "qc_filter_diagnostics.csv")

summary_df <- dplyr::bind_rows(summary_rows)
diagnostic_df <- dplyr::bind_rows(diagnostic_rows)

saveRDS(filtered_list, qc_filter_object_path)
write_tsv_local(summary_df, qc_filter_summary_path)
write_csv_local(diagnostic_df, qc_filter_diagnostics_path)

write_manifest_local(
  manifest_path = cfg$module_02b1_manifest_path,
  new_outputs = list(
    qc_filter_object = build_output_entry(qc_filter_object_path, "rds", module_name, "per-sample Seurat object list with qc_pass and qc_fail_reason metadata", base_dir = cfg$project_root),
    qc_filter_sample_summary = build_output_entry(qc_filter_summary_path, "tsv", module_name, "one row per sample QC filter summary", base_dir = cfg$project_root, schema = infer_schema_from_df(summary_df)),
    qc_filter_diagnostics = build_output_entry(qc_filter_diagnostics_path, "csv", module_name, "one row per cell QC diagnostics", base_dir = cfg$project_root, schema = infer_schema_from_df(diagnostic_df))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    manifest_02a = cfg$module_02a_manifest_path,
    ambient_object = ambient_rds,
    qc_threshold_file = cfg$qc_threshold_file
  ),
  version = cfg$module_version,
  depends_on = list(
    module_02a = cfg$module_02a_manifest_path
  )
)

message("02b1 完成。QC filter 对象: ", qc_filter_object_path)
