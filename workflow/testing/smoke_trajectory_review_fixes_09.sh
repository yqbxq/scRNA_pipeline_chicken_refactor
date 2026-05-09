#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
R_BIN="${R_BIN:-Rscript}"

if ! command -v "${R_BIN}" >/dev/null 2>&1; then
  echo "Missing Rscript for smoke_trajectory_review_fixes_09: ${R_BIN}" >&2
  exit 1
fi

cd "${PIPELINE_ROOT}"

"${R_BIN}" -e 'files <- c(
  "workflow/05single_script/helpers/load_helpers_09.R",
  "workflow/05single_script/helpers/project_paths_09.R",
  "workflow/05single_script/helpers/trajectory_utils.R",
  "workflow/05single_script/09a_trajectory_inputs.R",
  "workflow/05single_script/09b_trajectory_inputs_eda.R",
  "workflow/05single_script/09c_trajectory_root.R",
  "workflow/05single_script/09d_trajectory_slingshot.R",
  "workflow/05single_script/09e_trajectory_monocle3.R",
  "workflow/05single_script/09e2_trajectory_monocle2.R",
  "workflow/05single_script/09f_trajectory_paga_dpt.R",
  "workflow/05single_script/09g_trajectory_palantir.R",
  "workflow/05single_script/09h_trajectory_tradeseq.R",
  "workflow/05single_script/09i_trajectory_consensus.R",
  "workflow/05single_script/09j_trajectory_velocity_link.R",
  "workflow/05single_script/09k_trajectory_figures.R",
  "workflow/05single_script/09l_trajectory_eda.R",
  "workflow/05single_script/09m_trajectory_split_compare.R"
); for (f in files) parse(file = f)'

python3 -m py_compile workflow/04python/paga_dpt.py workflow/04python/palantir.py workflow/04python/extract_velocity_obs.py

bridge_out="$("${R_BIN}" - <<'EOF' 2>&1
source("workflow/05single_script/helpers/metadata_io.R")
source("workflow/05single_script/helpers/trajectory_utils.R")
tmp <- tempfile("bridge09_")
dir.create(tmp)
script <- file.path(tmp, "bridge.py")
writeLines(c(
  "#!/usr/bin/env python3",
  "print('BRIDGE_STDOUT_VISIBLE', flush=True)"
), script, useBytes = TRUE)
Sys.chmod(script, mode = "0755")
Sys.setenv(VELOCITY_ENV_PREFIX = "")
ok <- trajectory_run_python_bridge_09(script, file.path(tmp, "jobs.tsv"), file.path(tmp, "index.tsv"))
unlink(tmp, recursive = TRUE)
stopifnot(isTRUE(ok))
EOF
)"
if [[ "${bridge_out}" != *"BRIDGE_STDOUT_VISIBLE"* ]]; then
  echo "Python bridge stdout was not visible" >&2
  echo "${bridge_out}" >&2
  exit 1
fi

"${R_BIN}" - <<'EOF'
.script_dir <- normalizePath("workflow/05single_script", winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")
source(file.path(.script_dir, "helpers", "metadata_fanout_trajectory.R"), encoding = "UTF-8")
stopifnot(
  length(trajectory_pairs_cols_09) == 17L,
  identical(trajectory_pairs_cols_09, m3_trajectory_cols)
)
current_pairs <- read_tsv_optional("metadata/trajectory_pairs.tsv")
stopifnot(
  nrow(current_pairs) > 0,
  identical(colnames(current_pairs)[seq_along(trajectory_pairs_cols_09)], trajectory_pairs_cols_09)
)
row <- data.frame(methods_extra = "", tools_to_run = "slingshot,monocle3", stringsAsFactors = FALSE)
stopifnot(
  isTRUE(trajectory_method_enabled_09(row, "slingshot")),
  isTRUE(trajectory_method_enabled_09(row, "monocle3")),
  !isTRUE(trajectory_method_enabled_09(row, "paga_dpt")),
  identical(trajectory_method_disabled_reason_09(row, "paga_dpt"), "tools_to_run does not include paga_dpt")
)
row$methods_extra <- "+paga_dpt"
stopifnot(isTRUE(trajectory_method_enabled_09(row, "paga_dpt")))
row$methods_extra <- "-slingshot"
stopifnot(!isTRUE(trajectory_method_enabled_09(row, "slingshot")))
row$methods_extra <- ""
row$tools_to_run <- ""
stopifnot(!isTRUE(trajectory_method_enabled_09(row, "monocle2")))
row$methods_extra <- "+monocle2"
stopifnot(isTRUE(trajectory_method_enabled_09(row, "monocle2")))
EOF

"${R_BIN}" - <<'EOF'
.script_dir <- normalizePath("workflow/05single_script", winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")
exprs <- parse(file = "workflow/05single_script/09j_trajectory_velocity_link.R")
for (expr in exprs) {
  if (is.call(expr) && identical(as.character(expr[[1]]), "<-")) {
    lhs <- as.character(expr[[2]])
    if (lhs %in% c("regex_escape_09j", "velocity_h5ad_token_match_09j", "velocity_h5ad_candidates_09j")) {
      eval(expr, envir = .GlobalEnv)
    }
  }
}
tmp <- tempfile("velh5ad09_")
dir.create(tmp)
invisible(file.create(file.path(tmp, "H01_GC_baseline_extended.h5ad")))
invisible(file.create(file.path(tmp, "H01_GC_baseline.h5ad")))
invisible(file.create(file.path(tmp, "H02_GC_split__syf.h5ad")))
invisible(file.create(file.path(tmp, "H02_GC_split__f5.h5ad")))
cfg <- list(velocity_output_dir = tmp)
hits <- velocity_h5ad_candidates_09j(cfg, "H01_GC_baseline", "")
stopifnot(length(hits) == 1L, basename(hits[[1]]) == "H01_GC_baseline.h5ad")
hits <- velocity_h5ad_candidates_09j(cfg, "H02_GC_split", "syf")
stopifnot(length(hits) == 1L, basename(hits[[1]]) == "H02_GC_split__syf.h5ad")
unlink(tmp, recursive = TRUE)
EOF

"${R_BIN}" - <<'EOF'
.script_dir <- normalizePath("workflow/05single_script", winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")
tmp <- tempfile("units09_")
dir.create(tmp)
pairs <- as.data.frame(setNames(replicate(length(trajectory_pairs_cols_09), "", simplify = FALSE), trajectory_pairs_cols_09), stringsAsFactors = FALSE)
pairs[1, ] <- ""
pairs$trajectory_id <- "H01"
pairs$method <- "trajectory"
pairs$enabled <- "yes"
pairs$root_group <- "pair_root"
pairs$coarse_label_var <- "pair_label"
pairs$tools_to_run <- "slingshot,monocle3"
idx <- as.data.frame(setNames(replicate(length(trajectory_input_index_cols_09), "", simplify = FALSE), trajectory_input_index_cols_09), stringsAsFactors = FALSE)
idx[1, ] <- ""
idx$pair_id <- "H01"
idx$split_value <- "pooled"
idx$status <- "ok"
idx$root_group <- ""
idx$coarse_label_var <- ""
write_tsv_local(pairs, file.path(tmp, "pairs.tsv"))
write_tsv_local(idx, file.path(tmp, "input.tsv"))
cfg <- list(
  trajectory_pairs_sheet = file.path(tmp, "pairs.tsv"),
  module_09a_manifest_path = file.path(tmp, "missing_manifest.json"),
  trajectory_input_index_tsv = file.path(tmp, "input.tsv")
)
units <- trajectory_execution_units_09(cfg)
stopifnot(
  nrow(units) == 1L,
  identical(units$root_group[[1]], "pair_root"),
  identical(units$coarse_label_var[[1]], "pair_label"),
  identical(units$tools_to_run[[1]], "slingshot,monocle3")
)
stopifnot(!any(grepl("\\.pair$", colnames(units))))
unlink(tmp, recursive = TRUE)
EOF

"${R_BIN}" - <<'EOF'
.script_dir <- normalizePath("workflow/05single_script", winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "helpers", "load_helpers_09.R"), encoding = "UTF-8")
exprs <- parse(file = "workflow/05single_script/09a_trajectory_inputs.R")
for (expr in exprs) {
  if (is.call(expr) && identical(as.character(expr[[1]]), "<-")) {
    lhs <- as.character(expr[[2]])
    if (lhs == "outlier_thresholds_09") {
      eval(expr, envir = .GlobalEnv)
    }
  }
}
cfg <- list(
  trajectory_outlier_low_prob = 0.005,
  trajectory_outlier_high_prob = 0.995,
  trajectory_outlier_neighbor_cutoff = 5,
  trajectory_outlier_strict_low_prob = 0.03,
  trajectory_outlier_strict_high_prob = 0.97,
  trajectory_outlier_strict_neighbor_cutoff = 2
)
standard <- outlier_thresholds_09(cfg, "standard")
strict <- outlier_thresholds_09(cfg, "strict")
stopifnot(identical(standard$low_prob, 0.005), identical(standard$high_prob, 0.995), identical(standard$neighbor_cutoff, 5))
stopifnot(identical(strict$low_prob, 0.03), identical(strict$high_prob, 0.97), identical(strict$neighbor_cutoff, 2))
EOF

"${R_BIN}" - <<'EOF'
source("workflow/05single_script/helpers/metadata_io.R")
exprs <- parse(file = "workflow/05single_script/09c_trajectory_root.R")
for (expr in exprs) {
  if (is.call(expr) && identical(as.character(expr[[1]]), "<-")) {
    lhs <- as.character(expr[[2]])
    if (lhs %in% c("root_vote_summary_09c", "select_root_09c")) {
      eval(expr, envir = .GlobalEnv)
    }
  }
}
inference <- data.frame(
  method = c("prior", "cytotrace", "velocity"),
  recommended_root = c("pGC", "lGC", "lGC"),
  status = c("ok", "ok", "ok"),
  stringsAsFactors = FALSE
)
selected <- select_root_09c(inference, "pGC")
stopifnot(
  identical(selected$root, "pGC"),
  identical(selected$status, "warning_prior_disagreement"),
  identical(selected$prior_disagreement, "yes"),
  identical(selected$non_prior_vote_root, "lGC")
)
inference$recommended_root <- c("pGC", "pGC", "pGC")
selected <- select_root_09c(inference, "pGC")
stopifnot(
  identical(selected$status, "ok"),
  identical(selected$prior_disagreement, "no")
)
EOF

"${R_BIN}" - <<'EOF'
source("workflow/05single_script/helpers/metadata_io.R")
source("workflow/05single_script/helpers/layer_config_utils.R")
source("workflow/05single_script/helpers/project_paths_05.R")
source("workflow/05single_script/helpers/project_paths_06.R")
source("workflow/05single_script/helpers/project_paths_07.R")
source("workflow/05single_script/helpers/project_paths_09.R")
exprs <- parse(file = "workflow/05single_script/09c_trajectory_root.R")
for (expr in exprs) {
  if (is.call(expr) && identical(as.character(expr[[1]]), "<-")) {
    lhs <- as.character(expr[[2]])
    if (lhs == "velocity_root_row_09c") {
      eval(expr, envir = .GlobalEnv)
    }
  }
}
tmp <- tempfile("velroot09_")
dir.create(tmp)
writeLines(
  "cluster\troot_score\nsyf_root\t0.9\nf5_root\t0.1",
  file.path(tmp, "velocity_root_terminal_H02_GC_split__syf.tsv"),
  useBytes = TRUE
)
writeLines(
  "cluster\troot_score\nf5_root\t0.8\nsyf_root\t0.2",
  file.path(tmp, "velocity_root_terminal_H02_GC_split__f5.tsv"),
  useBytes = TRUE
)
cfg <- list(velocity_output_dir = tmp)
row <- velocity_root_row_09c(cfg, "H02_GC_split", "syf", "auto")
stopifnot(identical(row$status[[1]], "ok"), identical(row$recommended_root[[1]], "syf_root"))
unlink(file.path(tmp, "velocity_root_terminal_H02_GC_split__syf.tsv"))
writeLines(
  "cluster\troot_score\npooled_root\t1",
  file.path(tmp, "velocity_root_terminal_H02_GC_split.tsv"),
  useBytes = TRUE
)
row <- velocity_root_row_09c(cfg, "H02_GC_split", "syf", "auto")
stopifnot(identical(row$status[[1]], "skipped_no_velocity"))
row <- velocity_root_row_09c(cfg, "H02_GC_split", "pooled", "auto")
stopifnot(identical(row$status[[1]], "ok"), identical(row$recommended_root[[1]], "pooled_root"))
writeLines(
  "cluster\troot_score\npooled_alt\t1",
  file.path(tmp, "velocity_root_terminal_H02_GC_split__pooled.tsv"),
  useBytes = TRUE
)
row <- velocity_root_row_09c(cfg, "H02_GC_split", "pooled", "auto")
stopifnot(identical(row$status[[1]], "skipped_ambiguous_velocity"))
unlink(tmp, recursive = TRUE)
EOF

if rg -n '<<-' \
  workflow/05single_script/09*.R >/dev/null; then
  echo "Unexpected superassignment remains in module 09 scripts" >&2
  exit 1
fi

if rg -n 'Figure_5[FGH]' workflow/05single_script/09*.R >/dev/null; then
  echo "Legacy Figure_5* trajectory file names remain" >&2
  exit 1
fi

rg -n 'skipped_no_root' workflow/04python/paga_dpt.py workflow/04python/palantir.py >/dev/null

echo "smoke_trajectory_review_fixes_09: PASS"
