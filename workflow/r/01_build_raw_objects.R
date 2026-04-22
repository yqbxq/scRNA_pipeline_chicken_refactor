source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

sample_sheet <- main_sample_sheet(cfg)
inventory_df <- read_input_inventory(cfg)
sample_ids <- unique(sample_sheet$sample_id)
sample_ids <- sample_ids[nzchar(sample_ids)]

if (length(sample_ids) == 0) {
  stop("samples.tsv 中没有启用 run_main=yes 的样本。", call. = FALSE)
}

pick_expression_matrix <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  if ("Gene Expression" %in% names(x)) {
    return(x[["Gene Expression"]])
  }
  x[[1]]
}

sample_builds <- lapply(sample_ids, function(sample_id) {
  sample_row <- sample_sheet[sample_sheet$sample_id == sample_id, , drop = FALSE][1, ]
  inventory_row <- resolve_inventory_row(cfg, sample_id, inventory_df = inventory_df)
  matrix_dir <- file.path(cfg$data_dir, sample_id)
  if (!dir.exists(matrix_dir)) {
    stop(sprintf("缺少主流程输入矩阵目录: %s", matrix_dir), call. = FALSE)
  }

  counts <- pick_expression_matrix(Read10X(data.dir = matrix_dir))
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.features = 0
  )
  platform_resolved <- safe_get_field(inventory_row, sample_row, "platform_resolved", "platform", "generic_mex")
  gene_id_type_resolved <- safe_get_field(inventory_row, sample_row, "gene_id_type_resolved", "gene_id_type", "unknown")
  feature_name_profile <- safe_get_field(inventory_row, NULL, "feature_name_profile", default = "unknown")
  reference_version <- safe_get_field(inventory_row, sample_row, "reference_version")
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
  obj$biological_replicate <- if (nzchar(sample_row$biological_replicate)) sample_row$biological_replicate else sample_id
  obj$technical_replicate <- if (nzchar(sample_row$technical_replicate)) sample_row$technical_replicate else "1"
  obj$batch <- if (nzchar(sample_row$batch)) sample_row$batch else "default"
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
  list(object = obj, feature_context = qc_payload$feature_context)
})
names(sample_builds) <- sample_ids

sample_objects <- lapply(sample_builds, function(x) x$object)
sample_feature_contracts <- bind_rows(lapply(sample_ids, function(sample_id) {
  feature_context <- sample_builds[[sample_id]]$feature_context
  data.frame(
    sample_id = sample_id,
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
  group_by(sample_id, analysis_group, batch, platform, gene_id_type, feature_name_profile, reference_version, timepoint, tissue, chemistry) %>%
  summarise(
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
    .groups = "drop"
  ) %>%
  left_join(sample_feature_contracts, by = "sample_id")

saveRDS(raw_obj, file.path(cfg$checkpoint_dir, "00_raw_objects.rds"))
write.csv(sample_summary, file.path(cfg$table_dir, "raw_sample_summary.csv"), row.names = FALSE)

message("已保存 00_raw_objects.rds 和 raw_sample_summary.csv")
