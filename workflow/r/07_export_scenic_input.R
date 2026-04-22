source(file.path(Sys.getenv("PIPELINE_ROOT", Sys.getenv("PROJECT_ROOT")), "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(biomaRt)
})

cfg <- read_cfg()
prepare_dirs(cfg)
set.seed(cfg$random_seed)

input_rds <- file.path(cfg$checkpoint_dir, "03_after_annotation.rds")
if (!file.exists(input_rds)) {
  stop("缺少输入检查点: 03_after_annotation.rds", call. = FALSE)
}

obj <- readRDS(input_rds)
obj <- maybe_join_layers(obj)
counts <- get_assay_matrix(obj, type = "counts")

ortholog_map_file <- Sys.getenv("SCENIC_ORTHOLOG_MAP_FILE", "")

load_orthologs_from_file <- function(path) {
  orthologs_df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("external_gene_name", "hsapiens_homolog_associated_gene_name")
  if (!all(required_cols %in% colnames(orthologs_df))) {
    stop(
      sprintf(
        "SCENIC_ORTHOLOG_MAP_FILE 缺少必需列: %s",
        paste(setdiff(required_cols, colnames(orthologs_df)), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  orthologs_df %>%
    filter(
      external_gene_name != "",
      !is.na(external_gene_name),
      hsapiens_homolog_associated_gene_name != "",
      !is.na(hsapiens_homolog_associated_gene_name)
    ) %>%
    distinct(external_gene_name, .keep_all = TRUE)
}

query_orthologs_from_biomart <- function(genes, mirror) {
  chicken_mart <- useEnsembl(
    biomart = "ensembl",
    dataset = "ggallus_gene_ensembl",
    mirror = mirror
  )

  getBM(
    attributes = c("external_gene_name", "hsapiens_homolog_associated_gene_name"),
    filters = "external_gene_name",
    values = genes,
    mart = chicken_mart
  ) %>%
    filter(hsapiens_homolog_associated_gene_name != "") %>%
    distinct(external_gene_name, .keep_all = TRUE)
}

orthologs <- NULL
if (nzchar(ortholog_map_file) && file.exists(ortholog_map_file)) {
  message("使用本地 ortholog 映射文件: ", ortholog_map_file)
  orthologs <- load_orthologs_from_file(ortholog_map_file)
} else {
  orthologs <- query_orthologs_from_biomart(rownames(counts), cfg$ensembl_mirror)
}

map_vec <- setNames(orthologs$hsapiens_homolog_associated_gene_name, orthologs$external_gene_name)
mapped_idx <- which(rownames(counts) %in% names(map_vec))

counts_human <- counts[mapped_idx, , drop = FALSE]
rownames(counts_human) <- map_vec[rownames(counts_human)]
counts_human_agg <- rowsum(as.matrix(counts_human), group = rownames(counts_human))
counts_export <- t(counts_human_agg)

write.csv(orthologs, file.path(cfg$scenic_input_dir, "chicken_to_human_orthologs.csv"), row.names = FALSE)
write.csv(counts_export, file.path(cfg$scenic_input_dir, "expr_mat_human.csv"), row.names = TRUE)

message("SCENIC 输入矩阵已导出")
