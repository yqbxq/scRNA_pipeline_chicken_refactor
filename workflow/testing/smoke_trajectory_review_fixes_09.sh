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
  "workflow/05single_script/helpers/trajectory_utils.R",
  "workflow/05single_script/09b_trajectory_inputs_eda.R",
  "workflow/05single_script/09c_trajectory_root.R",
  "workflow/05single_script/09f_trajectory_paga_dpt.R",
  "workflow/05single_script/09g_trajectory_palantir.R",
  "workflow/05single_script/09l_trajectory_eda.R"
); for (f in files) parse(file = f)'

python3 -m py_compile workflow/04python/paga_dpt.py workflow/04python/palantir.py

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
  workflow/05single_script/09b_trajectory_inputs_eda.R \
  workflow/05single_script/09f_trajectory_paga_dpt.R \
  workflow/05single_script/09g_trajectory_palantir.R >/dev/null; then
  echo "Unexpected superassignment remains in 09b/09f/09g" >&2
  exit 1
fi

rg -n 'skipped_no_root' workflow/04python/paga_dpt.py workflow/04python/palantir.py >/dev/null

echo "smoke_trajectory_review_fixes_09: PASS"
