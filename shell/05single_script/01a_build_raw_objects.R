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
source_utf8(file.path(.script_dir, "helpers", "project_paths_01.R"))
source_utf8(file.path(.script_dir, "helpers", "gtf_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "manifest_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "report_utils.R"))
source_utf8(file.path(.script_dir, "helpers", "metadata_io.R"))
source_utf8(file.path(.script_dir, "helpers", "qc_utils.R"))

load_required_packages(c("Seurat", "dplyr", "tibble", "jsonlite"))

cfg <- get_single_script_config_01()
module_name <- "01a_build_raw_objects"
prepare_dirs_01(cfg)
set.seed(cfg$random_seed)

sample_sheet <- main_sample_sheet_local(cfg)
inventory_df <- read_input_inventory_local(cfg)
require_input_inventory_contract_local(cfg, inventory_df)
cc_genes_path <- require_ortholog_cc_genes_local(cfg)
sample_ids <- unique(sample_sheet$sample_id)
sample_ids <- sample_ids[nzchar(sample_ids)]

if (length(sample_ids) == 0) {
  stop("samples.tsv 中没有启用 run_main=yes 的样本。", call. = FALSE)
}

resolve_matrix_dir <- function(cfg, sample_id, inventory_row) {
  if (is.null(inventory_row)) {
    stop(sprintf("input_inventory.tsv 缺少样本 `%s` 的记录。请先运行 shell/03stages/audit_inputs.sh。", sample_id), call. = FALSE)
  }
  matrix_dir <- normalize_scalar_value(inventory_row$filtered_matrix_dir[1])
  if (!nzchar(matrix_dir)) {
    stop(sprintf("input_inventory.tsv 中样本 `%s` 的 filtered_matrix_dir 为空。请先运行 shell/03stages/audit_inputs.sh 和 standardize_inputs.sh。", sample_id), call. = FALSE)
  }
  if (!dir.exists(matrix_dir)) {
    stop(sprintf("input_inventory.tsv 中样本 `%s` 的 filtered_matrix_dir 不存在: %s", sample_id, matrix_dir), call. = FALSE)
  }
  if (!is_mex_matrix_dir_local(matrix_dir)) {
    stop(sprintf("input_inventory.tsv 中样本 `%s` 的 filtered_matrix_dir 不是有效 MEX 目录: %s", sample_id, matrix_dir), call. = FALSE)
  }
  matrix_dir
}

read_sample_counts <- function(matrix_dir, gene_id_type_resolved, feature_name_profile) {
  use_gene_ids <- tolower(gene_id_type_resolved) == "ensembl" ||
    tolower(feature_name_profile) == "ensembl"
  gene_column <- if (use_gene_ids) 1L else 2L
  pick_expression_matrix(Seurat::Read10X(data.dir = matrix_dir, gene.column = gene_column))
}

sample_builds <- lapply(sample_ids, function(sample_id) {
  sample_row <- sample_sheet[sample_sheet$sample_id == sample_id, , drop = FALSE][1, ]
  inventory_row <- resolve_inventory_row_local(cfg, sample_id, inventory_df = inventory_df)
  matrix_dir <- resolve_matrix_dir(cfg, sample_id, inventory_row)

  platform_resolved <- safe_get_field(inventory_row, sample_row, "platform_resolved", "platform", "generic_mex")
  gene_id_type_resolved <- safe_get_field(inventory_row, sample_row, "gene_id_type_resolved", "gene_id_type", "unknown")
  feature_name_profile <- safe_get_field(inventory_row, NULL, "feature_name_profile", default = "unknown")
  reference_version <- safe_get_field(inventory_row, sample_row, "reference_version")

  counts <- read_sample_counts(matrix_dir, gene_id_type_resolved, feature_name_profile)
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.features = 0
  )

  obj$sample_id <- sample_id
  obj$barcode_raw <- colnames(obj)
  resolved_group <- normalize_scalar_value(sample_row$group_id[1])
  if (!nzchar(resolved_group)) {
    resolved_group <- if (nzchar(sample_row$condition[1])) sample_row$condition[1] else sample_id
  }
  obj$group_id <- resolved_group
  obj$group <- resolved_group
  obj$analysis_group <- resolved_group
  obj$condition <- if (nzchar(sample_row$condition[1])) sample_row$condition[1] else resolved_group
  obj$biological_replicate <- if (nzchar(sample_row$biological_replicate[1])) sample_row$biological_replicate[1] else sample_id
  obj$technical_replicate <- if (nzchar(sample_row$technical_replicate[1])) sample_row$technical_replicate[1] else "1"
  obj$batch <- if (nzchar(sample_row$batch[1])) sample_row$batch[1] else "default"
  obj$platform <- platform_resolved
  obj$gene_id_type <- gene_id_type_resolved
  obj$feature_name_profile <- feature_name_profile
  obj$reference_version <- reference_version
  obj$timepoint <- normalize_scalar_value(sample_row$timepoint[1])
  obj$tissue <- normalize_scalar_value(sample_row$tissue[1])
  obj$chemistry <- normalize_scalar_value(sample_row$chemistry[1])

  qc_payload <- add_basic_qc_metrics(obj, cfg, declared_gene_id_type = gene_id_type_resolved)
  obj <- qc_payload$object
  obj$mito_feature_count <- qc_payload$feature_context$mito_feature_count
  obj$mito_detection_method <- qc_payload$feature_context$mito_detection_method
  obj$mito_detection_detail <- qc_payload$feature_context$mito_detection_detail
  obj$mito_warning_codes <- qc_payload$feature_context$mito_warning_codes_text
  obj$mito_warning_messages <- qc_payload$feature_context$mito_warning_messages_text
  obj$species_guess <- qc_payload$feature_context$species_guess
  obj$ribo_feature_count <- qc_payload$feature_context$ribo_feature_count
  obj$cell_cycle_s_feature_count <- qc_payload$feature_context$cell_cycle_s_feature_count
  obj$cell_cycle_g2m_feature_count <- qc_payload$feature_context$cell_cycle_g2m_feature_count
  obj$cell_cycle_gene_source <- qc_payload$feature_context$cell_cycle_gene_source

  list(
    object = obj,
    feature_context = qc_payload$feature_context,
    matrix_dir = matrix_dir
  )
})
names(sample_builds) <- sample_ids

sample_objects <- lapply(sample_builds, function(x) x$object)
sample_feature_contracts <- dplyr::bind_rows(lapply(sample_ids, function(sample_id) {
  feature_context <- sample_builds[[sample_id]]$feature_context
  data.frame(
    sample_id = sample_id,
    matrix_dir = sample_builds[[sample_id]]$matrix_dir,
    mito_detection_method = feature_context$mito_detection_method,
    mito_detection_detail = feature_context$mito_detection_detail,
    mito_warning_codes = feature_context$mito_warning_codes_text,
    mito_warning_messages = feature_context$mito_warning_messages_text,
    mito_detected_gene_names = feature_context$mito_detected_gene_names_text,
    mito_detected_feature_names = feature_context$mito_detected_feature_names_text,
    mito_reference_seqnames = feature_context$mito_reference_seqnames_text,
    species_guess = feature_context$species_guess,
    stringsAsFactors = FALSE
  )
}))

raw_obj <- sample_objects[[1]]
if (length(sample_objects) > 1) {
  raw_obj <- merge(
    x = sample_objects[[1]],
    y = sample_objects[2:length(sample_objects)],
    add.cell.ids = sample_ids
  )
}
raw_obj <- maybe_join_layers(raw_obj)
raw_obj@misc$sample_feature_contracts <- sample_feature_contracts

sample_summary <- raw_obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::group_by(sample_id, analysis_group, batch, platform, gene_id_type, feature_name_profile, reference_version, timepoint, tissue, chemistry) %>%
  dplyr::summarise(
    cells_raw = dplyr::n(),
    median_ncount = stats::median(nCount_RNA),
    median_nfeature = stats::median(nFeature_RNA),
    median_percent_mito = stats::median(percent.mito),
    median_percent_ribo = stats::median(percent.ribo),
    median_log10GenesPerUMI = stats::median(log10GenesPerUMI),
    mito_feature_count = dplyr::first(mito_feature_count),
    ribo_feature_count = dplyr::first(ribo_feature_count),
    cell_cycle_s_feature_count = dplyr::first(cell_cycle_s_feature_count),
    cell_cycle_g2m_feature_count = dplyr::first(cell_cycle_g2m_feature_count),
    cell_cycle_gene_source = dplyr::first(cell_cycle_gene_source),
    .groups = "drop"
  ) %>%
  dplyr::left_join(sample_feature_contracts, by = "sample_id")

raw_obj_path <- file.path(cfg$checkpoint_dir, "01a_raw_objects.rds")
sample_summary_path <- file.path(cfg$table_dir, "raw_sample_summary.csv")
feature_contracts_path <- file.path(cfg$table_dir, "sample_feature_contracts.tsv")

ensure_dir(dirname(raw_obj_path))
saveRDS(raw_obj, raw_obj_path)
write_csv_local(sample_summary, sample_summary_path)
write_tsv_local(sample_feature_contracts, feature_contracts_path)

write_manifest_local(
  manifest_path = cfg$module_01a_manifest_path,
  new_outputs = list(
    raw_object = build_output_entry(raw_obj_path, "rds", module_name, "merged raw Seurat object before QC filtering", base_dir = cfg$project_root),
    raw_sample_summary = build_output_entry(sample_summary_path, "csv", module_name, "one row per sample raw QC summary", base_dir = cfg$project_root, schema = infer_schema_from_df(sample_summary)),
    sample_feature_contracts = build_output_entry(feature_contracts_path, "tsv", module_name, "one row per sample feature/QC contract", base_dir = cfg$project_root, schema = infer_schema_from_df(sample_feature_contracts))
  ),
  module_name = module_name,
  base_dir = cfg$project_root,
  inputs = list(
    sample_sheet = cfg$sample_sheet,
    canonical_sample_sheet = cfg$canonical_sample_sheet,
    input_inventory_file = cfg$input_inventory_file,
    data_dir = cfg$data_dir,
    clean_gtf = cfg$clean_gtf,
    reference_gtf = cfg$reference_gtf,
    ortholog_manifest_path = cfg$ortholog_manifest_path,
    cell_cycle_genes = cc_genes_path
  ),
  version = cfg$module_version,
  depends_on = list(
    ortholog_cache = cfg$ortholog_manifest_path
  )
)

message("01a 完成。输出: ", raw_obj_path)
