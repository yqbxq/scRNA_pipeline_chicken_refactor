#!/usr/bin/env Rscript

.script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
  }
)

PIPELINE_ROOT <- Sys.getenv("PIPELINE_ROOT", unset = normalizePath(file.path(.script_dir, "..", "..", "..", ".."), winslash = "/", mustWork = FALSE))
source(file.path(PIPELINE_ROOT, "workflow/02lib/r/r_runtime_bootstrap.R"), encoding = "UTF-8")
source(file.path(PIPELINE_ROOT, "workflow/05single_script/spatial/helpers/spatial_common.R"), encoding = "UTF-8")

if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
  message("skip: Seurat and Matrix are required for smoke_spatial_loader.R")
  quit(status = 0)
}

tmp <- tempfile("spatial_loader_")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

counts <- Matrix::Matrix(c(5, 0, 2, 3, 4, 0), nrow = 3, sparse = TRUE)
Matrix::writeMM(counts, file.path(tmp, "matrix.mtx"))
write.table(
  data.frame(gene_id = c("g1", "g2", "g3"), gene_name = c("MT-CO1", "ACTB", "ND1"), type = "Gene Expression"),
  file.path(tmp, "features.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
writeLines(c("spot1", "spot2"), file.path(tmp, "barcodes.tsv"))
write.csv(
  data.frame(barcode = c("spot1", "spot2"), row = c(10, 11), col = c(20, 21), in_tissue = c(TRUE, TRUE)),
  file.path(tmp, "coords.csv"),
  row.names = FALSE
)

sample_row <- data.frame(
  sample_id = "st_smoke",
  condition = "syf",
  timepoint = "D0",
  platform = "generic",
  section_id = "smoke_section",
  chip_id = "smoke_chip",
  bundle_layout = "generic_spatial_matrix",
  standardized_outs = tmp,
  image_path = "",
  stringsAsFactors = FALSE
)
section_row <- data.frame(
  section_id = "smoke_section",
  chip_id = "smoke_chip",
  condition = "syf",
  platform = "generic",
  bundle_root = tmp,
  enabled = "yes",
  stringsAsFactors = FALSE
)
contract_row <- data.frame(platform = "generic", bundle_layout = "generic_spatial_matrix", stringsAsFactors = FALSE)

obj <- load_spatial_object(sample_row, section_row, contract_row)
stopifnot(Seurat::DefaultAssay(obj) == "Spatial")
stopifnot(length(obj$section_id) == ncol(obj))
stopifnot(identical(unique(obj$section_id), "smoke_section"))

mito <- detect_spatial_mito_features(obj)
stopifnot(isTRUE(mito$valid))
stopifnot(!is.na(mito$source))
obj <- inject_spatial_mito_qc(obj, mito)
metrics <- compute_spot_qc_metrics(obj)
stopifnot(nrow(metrics) == ncol(obj))
stopifnot(any(metrics$mito_qc_valid))

message("smoke_spatial_loader_ok")
