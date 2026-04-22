args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript smoke_feature_context.R /path/to/pipeline_root", call. = FALSE)
}

pipeline_root <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(pipeline_root, "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

tmp_dir <- tempfile("scrna_feature_contract_")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

make_cfg <- function(reference_gtf = "", clean_gtf = "", mito_gene_list_file = "", reference_version = "") {
  list(
    reference_gtf = reference_gtf,
    clean_gtf = clean_gtf,
    mito_gene_list_file = mito_gene_list_file,
    reference_version = reference_version
  )
}

gtf_path <- file.path(tmp_dir, "mini.gtf")
writeLines(
  c(
    'MT\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000001"; gene_name "MT-CO1"; gene_biotype "protein_coding";',
    '1\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000002"; gene_name "RPLP0"; gene_biotype "protein_coding";',
    '1\tmock\tgene\t1\t1000\t.\t+\t.\tgene_id "ENSGALG00000000003"; gene_name "FOO1"; gene_biotype "protein_coding";'
  ),
  con = gtf_path
)

counts <- Matrix::Matrix(
  c(
    8, 0,
    3, 5,
    2, 9
  ),
  nrow = 3,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts) <- c("ENSGALG00000000001", "ENSGALG00000000002", "ENSGALG00000000003")
colnames(counts) <- c("cellA", "cellB")

obj <- CreateSeuratObject(counts = counts, project = "smoke", min.features = 0)
payload <- add_basic_qc_metrics(
  obj,
  cfg = make_cfg(reference_gtf = gtf_path),
  declared_gene_id_type = "ensembl"
)

stopifnot(identical(payload$feature_context$feature_name_profile, "ensembl"))
stopifnot(identical(payload$feature_context$mito_detection_method, "gtf"))
stopifnot(payload$feature_context$mito_feature_count >= 1)
stopifnot(payload$feature_context$ribo_feature_count >= 1)
stopifnot(payload$object$percent.mito[1] > 0)
stopifnot(payload$object$percent.ribo[1] > 0)

mito_list_path <- file.path(tmp_dir, "mito_gene_list.txt")
writeLines(c("# one gene per line", "CUSTOM_MT"), con = mito_list_path)

counts_user <- Matrix::Matrix(
  c(
    10, 0,
    2, 3,
    1, 7
  ),
  nrow = 3,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts_user) <- c("CUSTOM-MT", "RPLP0", "FOO1")
colnames(counts_user) <- c("cellC", "cellD")

obj_user <- CreateSeuratObject(counts = counts_user, project = "smoke_user", min.features = 0)
payload_user <- add_basic_qc_metrics(
  obj_user,
  cfg = make_cfg(mito_gene_list_file = mito_list_path),
  declared_gene_id_type = "symbol"
)

stopifnot(identical(payload_user$feature_context$mito_detection_method, "user_list"))
stopifnot(payload_user$feature_context$mito_feature_count == 1)
stopifnot(payload_user$object$percent.mito[1] > 0)
stopifnot(payload_user$object$percent.ribo[1] > 0)

counts_prefix <- Matrix::Matrix(
  c(
    6, 0,
    4, 8,
    1, 5
  ),
  nrow = 3,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts_prefix) <- c("ND1", "COX3", "ATP6AP1")
colnames(counts_prefix) <- c("cellE", "cellF")

obj_prefix <- CreateSeuratObject(counts = counts_prefix, project = "smoke_prefix", min.features = 0)
payload_prefix <- add_basic_qc_metrics(
  obj_prefix,
  cfg = make_cfg(reference_version = "GRCg7b"),
  declared_gene_id_type = "symbol"
)

stopifnot(identical(payload_prefix$feature_context$mito_detection_method, "prefix_fallback"))
stopifnot(payload_prefix$feature_context$mito_feature_count == 2)
stopifnot(!"ATP6AP1" %in% payload_prefix$feature_context$mito_features)
stopifnot(grepl("mito_prefix_fallback", payload_prefix$feature_context$mito_warning_codes_text, fixed = TRUE))
stopifnot(payload_prefix$object$percent.mito[1] > 0)

cat("smoke_feature_context: PASS\n")
