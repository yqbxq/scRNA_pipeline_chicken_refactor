#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKFLOW_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

usage() {
  cat <<'EOF'
Usage:
  bash workflow/03stages/90_export_h5ad.sh \
    --upstream-manifest <manifest.json> \
    --output-key <manifest output key> \
    --modality scrna|spatial \
    --target-dir <output directory>
EOF
}

upstream_manifest=""
output_key=""
modality="scrna"
target_dir=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --upstream-manifest)
      upstream_manifest="$2"; shift 2 ;;
    --output-key)
      output_key="$2"; shift 2 ;;
    --modality)
      modality="$2"; shift 2 ;;
    --target-dir)
      target_dir="$2"; shift 2 ;;
    -h|--help)
      usage; exit 0 ;;
    *)
      echo "unknown argument: $1" >&2
      usage >&2
      exit 2 ;;
  esac
done

[[ -s "${upstream_manifest}" ]] || die "missing upstream manifest: ${upstream_manifest}"
[[ -n "${output_key}" ]] || die "--output-key is required"
[[ "${modality}" == "scrna" || "${modality}" == "spatial" ]] || die "--modality must be scrna or spatial"
[[ -n "${target_dir}" ]] || die "--target-dir is required"

source_rds="$(manifest_output_path "${upstream_manifest}" "${output_key}")"
[[ -s "${source_rds}" ]] || die "upstream output is missing: ${source_rds}"
ensure_dir "${target_dir}"

target_manifest="${target_dir}/_manifest.json"
source_module="$(
  python_bin="$(detect_python)"
  "${python_bin}" - "${upstream_manifest}" <<'PY'
import json
import sys
from pathlib import Path
path = Path(sys.argv[1])
try:
    data = json.loads(path.read_text(encoding="utf-8"))
    print(data.get("module") or path.parent.name)
except Exception:
    print(path.parent.name)
PY
)"
if [[ "${modality}" == "spatial" ]]; then
  export_script="${WORKFLOW_ROOT}/05single_script/spatial/90b_export_h5ad_spatial.R"
else
  export_script="${WORKFLOW_ROOT}/05single_script/90a_export_h5ad_scrna.R"
fi

run_h5ad_export() {
  local script_path="$1"
  export H5AD_UPSTREAM_MANIFEST="${upstream_manifest}"
  export H5AD_UPSTREAM_RDS="${source_rds}"
  export H5AD_OUTPUT_KEY="${output_key}"
  export H5AD_MODALITY="${modality}"
  export H5AD_TARGET_DIR="${target_dir}"
  export H5AD_TARGET_MANIFEST="${target_manifest}"
  export H5AD_SOURCE_MODULE="${source_module}"
  export H5AD_EXPORT_MODULE="90_export_h5ad"
  if [[ "${H5AD_EXPORT_DIRECT_RSCRIPT:-no}" == "yes" ]]; then
    Rscript "${script_path}"
    return
  fi
  if [[ "${modality}" == "spatial" ]]; then
    run_r_spatial "${script_path}"
  else
    run_r_main "${script_path}"
  fi
}

run_stage_if_stale_with_runner \
  run_h5ad_export \
  --script "${export_script}" \
  --manifest "${target_manifest}" \
  --inputs "${upstream_manifest}" "${source_rds}" "${H5AD_CONTRACT_FILE}" \
  --params H5AD_EXPORT_BACKEND,H5AD_EXPORT_ASSAY,H5AD_EXPORT_LAYERS,H5AD_EXPORT_REDUCTIONS,H5AD_EXPORT_COMPRESSION,H5AD_EXPORT_SPATIAL_PER_SECTION,H5AD_EXPORT_CONTRACT_FAIL_ON,H5AD_CONTRACT_FAIL_ON

require_manifest_output "${target_manifest}" "export_summary_tsv" >/dev/null
echo "90_export_h5ad_ok: ${target_manifest}"
