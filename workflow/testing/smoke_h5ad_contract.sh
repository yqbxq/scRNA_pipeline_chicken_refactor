#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d)"

cleanup() {
  rm -rf "${TMP_ROOT}"
}
trap cleanup EXIT

CONTRACT_FILE="${PIPELINE_ROOT}/metadata/h5ad_export_contract.tsv.template"

Rscript - "${PIPELINE_ROOT}" "${CONTRACT_FILE}" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
pipeline_root <- args[[1]]
contract_file <- args[[2]]
source(file.path(pipeline_root, "workflow/05single_script/helpers/h5ad_contract_utils.R"))

assert_true <- function(value, label) {
  if (!isTRUE(value)) stop(sprintf("assertion failed: %s", label), call. = FALSE)
}

expect_error <- function(expr, label) {
  ok <- FALSE
  tryCatch({
    force(expr)
  }, error = function(e) {
    ok <<- TRUE
  })
  assert_true(ok, label)
}

expect_warning <- function(expr, label) {
  warned <- FALSE
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  assert_true(warned, label)
  value
}

base_scrna_mock <- function() {
  list(
    obs = data.frame(
      sample_id = c("S1", "S2"),
      condition = c("syf", "f5"),
      celltype_l1 = c("GC", "TC"),
      stringsAsFactors = FALSE
    ),
    var = data.frame(
      gene_id = c("ENSG1", "ENSG2"),
      gene_symbol = c("G1", "G2"),
      stringsAsFactors = FALSE
    ),
    obsm = list(
      X_pca = matrix(c(1, 2, 3, 4), nrow = 2),
      X_umap = matrix(c(1, 2, 3, 4), nrow = 2)
    ),
    layers = list(
      counts = matrix(as.integer(c(1, 0, 2, 3)), nrow = 2)
    ),
    uns = list(
      module = "smoke",
      export_timestamp = "2026-05-25T00:00:00Z"
    )
  )
}

base_spatial_mock <- function() {
  obj <- base_scrna_mock()
  obj$obs$section_id <- c("sec1", "sec1")
  obj$obsm$spatial <- matrix(c(10, 20, 11, 21), nrow = 2)
  obj$spatial <- list(
    tissue_positions = data.frame(spot_id = c("s1", "s2"), x = c(10, 11), y = c(20, 21)),
    scale_factors = list(spot_diameter_fullres = 65)
  )
  obj
}

scrna_contract <- load_h5ad_contract(contract_file, modality = "scrna")
spatial_contract <- load_h5ad_contract(contract_file, modality = "spatial")

c1 <- check_h5ad_contract_seurat(base_scrna_mock(), scrna_contract, fail_on = "error")
assert_true(c1$passed, "C1 complete scRNA passed")

c2 <- base_scrna_mock()
c2$obs$sample_id <- NULL
expect_error(check_h5ad_contract_seurat(c2, scrna_contract, fail_on = "error"), "C2 missing sample_id stops")
c2_warn <- expect_warning(check_h5ad_contract_seurat(c2, scrna_contract, fail_on = "warn"), "C2 missing sample_id warns")
assert_true(!c2_warn$passed, "C2 warn mode returns failed primary")

c3 <- base_scrna_mock()
c3$obsm$X_umap <- NULL
c3_warn <- expect_warning(check_h5ad_contract_seurat(c3, scrna_contract, fail_on = "warn"), "C3 missing exploratory UMAP warns")
assert_true(c3_warn$passed, "C3 exploratory missing still passed")

c4 <- base_scrna_mock()
c4$layers$counts <- matrix(c("a", "b", "c", "d"), nrow = 2)
expect_error(check_h5ad_contract_seurat(c4, scrna_contract, fail_on = "error"), "C4 bad counts dtype stops")
c4_warn <- expect_warning(check_h5ad_contract_seurat(c4, scrna_contract, fail_on = "warn"), "C4 bad counts dtype warns")
assert_true(!c4_warn$passed, "C4 warn mode returns failed primary")

c5 <- check_h5ad_contract_seurat(base_spatial_mock(), spatial_contract, fail_on = "error")
assert_true(c5$passed, "C5 complete spatial passed")

c6 <- base_spatial_mock()
c6$spatial$tissue_positions <- NULL
expect_error(check_h5ad_contract_seurat(c6, spatial_contract, fail_on = "error"), "C6 missing tissue_positions stops")
c6_warn <- expect_warning(check_h5ad_contract_seurat(c6, spatial_contract, fail_on = "warn"), "C6 missing tissue_positions warns")
assert_true(!c6_warn$passed, "C6 warn mode returns failed primary")

manifest <- contract_violations_to_manifest(c6_warn)
assert_true(identical(manifest$contract_passed, FALSE), "manifest passed flag")
assert_true(manifest$contract_violations_count > 0, "manifest violation count")

cat("r_h5ad_contract_smoke_ok\n")
RS

PYTHONPATH="${PIPELINE_ROOT}/workflow/04python" python3 - "${CONTRACT_FILE}" <<'PY'
import sys
from helpers.h5ad_contract_utils import load_h5ad_contract, check_h5ad_contract_anndata

contract = load_h5ad_contract(sys.argv[1], modality="scrna")
mock = {
    "obs": {"sample_id": ["S1"], "condition": ["syf"], "celltype_l1": ["GC"]},
    "var": {"gene_id": ["g1"], "gene_symbol": ["G1"]},
    "obsm": {"X_pca": [[0.0, 1.0]], "X_umap": [[0.0, 1.0]]},
    "layers": {"counts": [1]},
    "uns": {"module": "smoke", "export_timestamp": "2026-05-25T00:00:00Z"},
}
result = check_h5ad_contract_anndata(mock, contract, fail_on="error")
assert result["passed"]
mock["obs"].pop("sample_id")
try:
    check_h5ad_contract_anndata(mock, contract, fail_on="error")
except ValueError:
    pass
else:
    raise SystemExit("expected missing sample_id to fail")
print("py_h5ad_contract_smoke_ok")
PY

test -s "${CONTRACT_FILE}"
echo "smoke_h5ad_contract_ok"
