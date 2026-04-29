#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

MODULE_08A_MANIFEST="${MANIFEST_DIR}/08a_scenic_export/_manifest.json"
MODULE_08B_MANIFEST="${MANIFEST_DIR}/08b_scenic_grn/_manifest.json"
SCENIC_EXPORT_INDEX="${TABLE_DIR}/regulation/scenic/scenic_export_index.tsv"
SCENIC_GRN_INDEX="${TABLE_DIR}/regulation/scenic/scenic_grn_index.tsv"

safe_layer_id_08b() {
  printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

ensure_scenic_resources_08b() {
  local missing=0
  for path in "${SCENIC_TF_LIST}" "${SCENIC_MOTIF_ANN}" "${SCENIC_DB_500BP}" "${SCENIC_DB_10KB}" "${SCENIC_RESOURCE_MANIFEST}"; do
    [[ -s "${path}" ]] || missing=1
  done
  if [[ "${missing}" == "0" ]]; then
    return 0
  fi
  env SCENIC_RESOURCES_ONLY=yes bash "${PIPELINE_ROOT}/workflow/03stages/08_regulation.sh"
}

ensure_scenic_resources_08b
require_manifest_output "${MODULE_08A_MANIFEST}" "scenic_export_index_tsv" >/dev/null
[[ -s "${SCENIC_EXPORT_INDEX}" ]] || die "缺少 08a export index: ${SCENIC_EXPORT_INDEX}"

ensure_dir "$(dirname "${SCENIC_GRN_INDEX}")" "$(dirname "${MODULE_08B_MANIFEST}")"

{
  printf 'layer_id\tlayer_role\texpr_mat_human_csv\tadjacencies_tsv\tstatus\treason\n'
  {
    read -r _header
    while IFS=$'\t' read -r layer_id layer_role object_rds expr_mat_human_csv ortholog_csv cell_n chicken_gene_n mapped_chicken_gene_n human_gene_n rest; do
      [[ -n "${layer_id}" ]] || continue
      [[ -s "${expr_mat_human_csv}" ]] || die "缺少 08a 表达矩阵: ${expr_mat_human_csv}"
      safe_layer="$(safe_layer_id_08b "${layer_id}")"
      out_dir="${SCENIC_OUTPUT_DIR}/${safe_layer}"
      adjacency_path="${out_dir}/adjacencies.tsv"
      ensure_dir "${out_dir}"

      echo "08b pySCENIC GRN layer: ${layer_id}"
      run_pyscenic pyscenic grn \
        --num_workers "${SCENIC_THREADS}" \
        --method grnboost2 \
        -o "${adjacency_path}" \
        "${expr_mat_human_csv}" \
        "${SCENIC_TF_LIST}"

      if [[ "${layer_id}" == "${PANORAMA_LAYER_ID}" ]]; then
        cp -f "${adjacency_path}" "${SCENIC_OUTPUT_DIR}/adjacencies.tsv"
      fi

      printf '%s\t%s\t%s\t%s\tcompleted\t\n' \
        "${layer_id}" \
        "${layer_role}" \
        "${expr_mat_human_csv}" \
        "${adjacency_path}"
    done
  } < "${SCENIC_EXPORT_INDEX}"
} > "${SCENIC_GRN_INDEX}"

python_bin="$(detect_python)"
env \
  MODULE_08B_MANIFEST="${MODULE_08B_MANIFEST}" \
  PROJECT_ROOT="${PROJECT_ROOT}" \
  MODULE_VERSION="${MODULE_08_VERSION}" \
  SCENIC_GRN_INDEX="${SCENIC_GRN_INDEX}" \
  SCENIC_RESOURCE_MANIFEST="${SCENIC_RESOURCE_MANIFEST}" \
  MODULE_08A_MANIFEST="${MODULE_08A_MANIFEST}" \
  "${python_bin}" - <<'PY'
import csv
import json
import os
from datetime import datetime
from pathlib import Path

manifest_path = Path(os.environ["MODULE_08B_MANIFEST"])
project_root = Path(os.environ["PROJECT_ROOT"]).resolve()
index_path = Path(os.environ["SCENIC_GRN_INDEX"]).resolve()
resource_manifest = Path(os.environ["SCENIC_RESOURCE_MANIFEST"]).resolve()

def rel(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(project_root))
    except ValueError:
        return str(path)

outputs = {
    "scenic_grn_index_tsv": {
        "path": rel(index_path),
        "type": "tsv",
        "produced_by": "08b_scenic_grn",
        "row_semantics": "one row per pySCENIC GRN task",
    },
    "scenic_resources_manifest": {
        "path": rel(resource_manifest),
        "type": "json",
        "produced_by": "08b_scenic_grn",
        "row_semantics": "SCENIC resource file fingerprint manifest",
    },
}

with index_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        layer_id = row.get("layer_id", "")
        adj = row.get("adjacencies_tsv", "")
        if not layer_id or not adj:
            continue
        safe_layer = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in layer_id)
        outputs[f"adjacencies_tsv_{safe_layer}"] = {
            "path": rel(adj),
            "type": "tsv",
            "produced_by": "08b_scenic_grn",
            "row_semantics": "pySCENIC TF-target adjacency table",
        }

manifest = {
    "module": "08b_scenic_grn",
    "version": os.environ["MODULE_VERSION"],
    "timestamp": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
    "base_dir": str(project_root),
    "inputs": {
        "module_08a": os.environ["MODULE_08A_MANIFEST"],
        "scenic_tf_list": os.environ.get("SCENIC_TF_LIST", ""),
        "scenic_threads": os.environ.get("SCENIC_THREADS", ""),
    },
    "outputs": outputs,
    "depends_on": {
        "module_08a": os.environ["MODULE_08A_MANIFEST"],
        "scenic_resources_manifest": str(resource_manifest),
    },
}
manifest_path.parent.mkdir(parents=True, exist_ok=True)
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY

echo "08b completed. index: ${SCENIC_GRN_INDEX}"
