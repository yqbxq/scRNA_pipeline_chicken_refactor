#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript - <<'RS'
source("workflow/05single_script/helpers/runtime_utils.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/config.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/multinichenet_wrapper.R", encoding = "UTF-8")
source("workflow/05single_script/helpers/nichenet_legacy_wrapper.R", encoding = "UTF-8")

cfg <- list(sample_sheet = tempfile(fileext = ".tsv"))
writeLines(c("sample_id\tcondition", "s1\tsyf", "s2\tsyf", "s3\tf5", "s4\tf5"), cfg$sample_sheet)
stopifnot(identical(compute_min_samples_per_group_07c(cfg$sample_sheet), 2L))

Sys.setenv(NICHENET_MODE = "nichenet_legacy")
stopifnot(identical(determine_07c_mode(cfg), "nichenet_legacy"))
stopifnot(identical(legacy_nichenet_mode_label_07c(), "nichenet_legacy"))

Sys.setenv(NICHENET_MODE = "skip")
stopifnot(identical(determine_07c_mode(cfg), "skip"))

links <- data.frame(ligand = c("L1", "L2"), stringsAsFactors = FALSE)
mat <- matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
rownames(mat) <- c("L1", "L2")
colnames(mat) <- c("T1", "T2")
out <- attach_receiver_de_target_count_07c(links, mat, receiver_de_genes = "T1", top_n = 1L)
stopifnot(out$n_targets_in_receiver_de[1] == 1L)
RS

rg -q "multinichenet_mode_used" workflow/05single_script/07c_nichenet.R
rg -q "method_evidence_class" workflow/05single_script/07c_nichenet.R
rg -q "NICHENET_MODE" workflow/03stages/07_communication.sh

echo "smoke_07c_nichenet_fallback_ok"
