args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript smoke_ambient_doublet_logic.R /path/to/pipeline_root", call. = FALSE)
}

pipeline_root <- normalizePath(args[[1]], mustWork = TRUE)
source(file.path(pipeline_root, "workflow", "r", "pipeline_common.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cfg <- list(
  ambient_apply_policy = "manual",
  ambient_recommend_min_contamination = 0.05,
  cellbender_mode = "stub",
  cellbender_fpr = 0.01,
  cellbender_cuda = "yes",
  cellbender_extra_args = "",
  checkpoint_dir = tempdir()
)

soupx_ready_row <- data.frame(
  ambient_preferred_method = "soupx",
  ambient_fallback_method = "decontx",
  ambient_apply_default = "manual",
  ambient_soupx_ready = "true",
  ambient_decontx_ready = "true",
  ambient_cellbender_ready = "true",
  raw_matrix_available = "true",
  raw_matrix_kind = "mex_dir",
  ambient_notes = "cellbender_stub_only",
  stringsAsFactors = FALSE
)

fallback_row <- data.frame(
  ambient_preferred_method = "decontx",
  ambient_fallback_method = "none",
  ambient_apply_default = "manual",
  ambient_soupx_ready = "false",
  ambient_decontx_ready = "true",
  ambient_cellbender_ready = "false",
  raw_matrix_available = "false",
  raw_matrix_kind = "",
  ambient_notes = "soupx_requires_raw_matrix;decontx_fallback_only",
  stringsAsFactors = FALSE
)

chosen_soupx <- choose_ambient_method(soupx_ready_row, cfg)
stopifnot(identical(chosen_soupx$preferred_method, "soupx"))
stopifnot(isTRUE(chosen_soupx$soupx_ready))
stopifnot(isTRUE(chosen_soupx$cellbender_ready))

chosen_fallback <- choose_ambient_method(fallback_row, cfg)
stopifnot(identical(chosen_fallback$preferred_method, "decontx"))
stopifnot(!chosen_fallback$soupx_ready)
stopifnot(chosen_fallback$decontx_ready)

stopifnot(abs(estimate_expected_doublet_rate(1000, 0.008) - 0.008) < 1e-9)
stopifnot(abs(estimate_expected_doublet_rate(10000, 0.008) - 0.08) < 1e-9)
stopifnot(estimate_expected_doublet_count(1000, 0.008) >= 1)

counts <- Matrix::Matrix(
  c(
    3, 0, 2,
    1, 4, 0
  ),
  nrow = 2,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts) <- c("GENE1", "GENE2")
colnames(counts) <- c("cellA", "cellB", "cellC")
obj <- CreateSeuratObject(counts = counts, project = "smoke", min.features = 0)
inventory_row <- data.frame(raw_matrix_dir = "/tmp/mock_raw.h5", stringsAsFactors = FALSE)
stub_row <- build_cellbender_stub("S1", obj, inventory_row, cfg, raw_counts = counts)

stopifnot(stub_row$expected_cells[[1]] == 3)
stopifnot(stub_row$total_droplets_included[[1]] == 3)
stopifnot(grepl("cellbender", stub_row$command_stub[[1]], fixed = TRUE))
stopifnot(grepl("--expected-cells", stub_row$command_stub[[1]], fixed = TRUE))

stopifnot(ambient_recommend_apply(0.10, cfg))
stopifnot(!ambient_recommend_apply(0.01, cfg))

cat("smoke_ambient_doublet_logic: PASS\n")
