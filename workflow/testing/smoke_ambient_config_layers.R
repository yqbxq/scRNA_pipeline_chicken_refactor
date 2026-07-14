args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("usage: Rscript smoke_ambient_config_layers.R /path/to/pipeline_root", call. = FALSE)
}

pipeline_root <- normalizePath(args[[1]], mustWork = TRUE)
helper_dir <- file.path(pipeline_root, "workflow", "05single_script", "helpers")
source(file.path(helper_dir, "runtime_utils.R"), encoding = "UTF-8")
source(file.path(helper_dir, "project_paths_02.R"), encoding = "UTF-8")
source(file.path(helper_dir, "gtf_utils.R"), encoding = "UTF-8")
source(file.path(helper_dir, "metadata_io.R"), encoding = "UTF-8")
source(file.path(helper_dir, "qc_utils.R"), encoding = "UTF-8")
source(file.path(helper_dir, "ambient_utils.R"), encoding = "UTF-8")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

readiness_row <- data.frame(
  ambient_preferred_method = "soupx",
  ambient_fallback_method = "decontx",
  ambient_apply_default = "manual",
  ambient_soupx_ready = "true",
  ambient_decontx_ready = "true",
  ambient_cellbender_ready = "false",
  raw_matrix_available = "true",
  raw_matrix_kind = "mex_dir",
  ambient_notes = "",
  stringsAsFactors = FALSE
)

explicit_cfg <- list(
  ambient_primary_method = "soupx",
  ambient_fallback_method = "decontx",
  ambient_apply_policy = "recommended",
  ambient_apply_policy_explicit = TRUE
)
chosen_explicit <- choose_ambient_method(readiness_row, explicit_cfg)
stopifnot(identical(chosen_explicit$apply_policy, "recommended"))

implicit_cfg <- explicit_cfg
implicit_cfg$ambient_apply_policy <- "recommended"
implicit_cfg$ambient_apply_policy_explicit <- FALSE
chosen_implicit <- choose_ambient_method(readiness_row, implicit_cfg)
stopifnot(identical(chosen_implicit$apply_policy, "manual"))

counts <- Matrix::Matrix(
  c(
    1, 0, 2, 0,
    0, 3, 0, 4
  ),
  nrow = 2,
  byrow = TRUE,
  sparse = TRUE
)
rownames(counts) <- c("GENE1", "GENE2")
colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
obj <- CreateSeuratObject(counts = counts, project = "smoke", min.features = 0)
obj[["RNA"]] <- CreateAssayObject(counts = counts)
joined <- maybe_join_layers(obj)
stopifnot(inherits(joined[["RNA"]], "Assay"))
stopifnot(identical(ncol(joined), ncol(obj)))

cat("smoke_ambient_config_layers: PASS\n")
