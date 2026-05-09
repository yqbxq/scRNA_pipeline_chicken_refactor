#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

bash -n workflow/05single_script/10a_run_velocyto.sh
bash -n workflow/03stages/10_velocity.sh

Rscript -e 'invisible(parse(file = "workflow/05single_script/10b_prepare_velocity_reference.R")); invisible(parse(file = "workflow/05single_script/helpers/project_paths_10.R")); invisible(parse(file = "workflow/05single_script/helpers/velocity_utils_10.R"))'

R --slave <<'RS'
.script_dir <- normalizePath("workflow/05single_script", winslash = "/", mustWork = TRUE)
source(file.path(.script_dir, "helpers", "load_helpers_10.R"), encoding = "UTF-8")

cfg <- get_single_script_config_10()
stopifnot(grepl("10a_run_velocyto", cfg$module_10a_manifest_path, fixed = TRUE))
stopifnot(grepl("10b_prepare_velocity_reference", cfg$module_10b_manifest_path, fixed = TRUE))
stopifnot(identical(velocity_unit_file_id_10("H06_GC_velocity_split__velocity", "syf"), "H06_GC_velocity_split__velocity__syf"))
stopifnot(identical(velocity_reference_output_key_10("velocity_umap", "H06", "f5"), "velocity_umap__H06__f5"))

pairs <- velocity_read_pairs_10(cfg)
if (nrow(pairs) > 0) {
  stopifnot(all(pairs$method == "velocity"))
}
RS

tmp_root="$(mktemp -d)"
trap 'rm -rf "${tmp_root}"' EXIT
mkdir -p "${tmp_root}/metadata" "${tmp_root}/results"
printf 'sample_id\trun_velocity\n' > "${tmp_root}/metadata/samples.tsv"

PROJECT_ROOT="${tmp_root}" \
RESULTS_DIR="${tmp_root}/results" \
TABLE_DIR="${tmp_root}/results/tables" \
MANIFEST_DIR="${tmp_root}/results/manifests" \
DATA_DIR="${tmp_root}/data" \
METADATA_DIR="${tmp_root}/metadata" \
SAMPLE_SHEET="${tmp_root}/metadata/samples.tsv" \
CANONICAL_SAMPLE_SHEET="${tmp_root}/metadata/samples.canonical.tsv" \
VELOCITY_SAMPLE_IDS="" \
RAW_SAMPLES="" \
SAMPLE_NAMES="" \
bash workflow/05single_script/10a_run_velocyto.sh

python3 - <<PY
import json
from pathlib import Path

root = Path("${tmp_root}")
manifest = json.loads((root / "results/manifests/10a_run_velocyto/_manifest.json").read_text(encoding="utf-8"))
assert "velocity_loom_index" in manifest["outputs"]
assert (root / "results/tables/velocity/inputs/velocity_loom_index.tsv").exists()
PY

echo "smoke_velocity_inputs_10 passed"
