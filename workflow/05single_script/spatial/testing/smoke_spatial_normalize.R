#!/usr/bin/env Rscript

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT is required", call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

assert <- function(ok, msg) {
  if (!isTRUE(ok)) {
    stop(msg, call. = FALSE)
  }
}

same_bytes <- function(a, b) {
  if (!file.exists(a) || !file.exists(b) || file.info(a)$size != file.info(b)$size) {
    return(FALSE)
  }
  identical(readBin(a, "raw", n = file.info(a)$size), readBin(b, "raw", n = file.info(b)$size))
}

qc_summary <- read_tsv(file.path(project_root, "results/spatial/tables/spatial_post_qc/qc_filter_summary.tsv"))
assert(nrow(qc_summary) == 3L, "expected three sections in qc_filter_summary.tsv")
assert(all(qc_summary$retention > 0.8), "fixture should retain more than 80% of spots")
assert(all(file.exists(file.path(project_root, qc_summary$post_qc_rds))), "post-QC RDS files missing")
assert(length(list.files(file.path(project_root, "reports/eda/spatial_post_qc/figures"), pattern = "^qc_filter_mask_.*\\.png$")) == 3L, "mask figure count mismatch")

triage <- read_tsv(file.path(project_root, "results/spatial/tables/spatial_post_qc/post_filter_triage.tsv"))
assert(identical(colnames(triage), c("section_id", "excessive_drop", "boundary_drop", "drift_from_pre_qc", "unbalanced_sections")), "post_filter_triage.tsv column order changed")
assert(!any(tolower(triage$excessive_drop) == "true"), "fixture should not trigger excessive_drop")
assert(!any(tolower(triage$unbalanced_sections) == "true"), "fixture should not trigger unbalanced_sections")
assert(file.exists(file.path(project_root, "reports/eda/spatial_post_qc/report.md")), "post-QC report missing")

selected <- read_tsv(file.path(project_root, "reports/spatial/normalization_compare/selected_method.tsv"))
assert(nrow(selected) == 3L, "expected selected method row per section")
assert(all(selected$final_choice == "m3_sct_v2"), "default selected method should be m3_sct_v2")
for (idx in seq_len(nrow(selected))) {
  canonical <- file.path(project_root, selected$canonical_rds[[idx]])
  variant <- file.path(
    project_root,
    "results/spatial/checkpoints/normalization_variants",
    selected$sample_id[[idx]],
    paste0(selected$final_choice[[idx]], "_post_norm.rds")
  )
  assert(same_bytes(canonical, variant), sprintf("canonical normalized object is not bytewise equal for %s", selected$sample_id[[idx]]))
}

hvg <- read_tsv(file.path(project_root, "reports/spatial/normalization_compare/hvg_iou_matrix.tsv"))
for (section_id in unique(hvg$section_id)) {
  sec <- hvg[hvg$section_id == section_id, , drop = FALSE]
  methods <- sec$method
  mat <- as.matrix(sec[, methods, drop = FALSE])
  storage.mode(mat) <- "numeric"
  assert(isTRUE(all.equal(mat, t(mat), check.attributes = FALSE)), sprintf("HVG IoU matrix is not symmetric for %s", section_id))
  assert(all(diag(mat) == 1), sprintf("HVG IoU matrix diagonal is not 1 for %s", section_id))
}

timing <- read_tsv(file.path(project_root, "reports/spatial/normalization_compare/method_timing.tsv"))
assert(any(timing$method == "m4_pearson_residuals" & timing$status == "skipped"), "m4 should be marked skipped when not requested")
assert(file.exists(file.path(project_root, "reports/spatial/normalization_compare/report.md")), "normalization report missing")

message("smoke_spatial_normalize_r_ok")
