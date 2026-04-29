#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${WORKFLOW_ROOT}/.." && pwd)}"

source "${WORKFLOW_ROOT}/02lib/common.sh"

NICHENET_RESOURCE_DIR="${NICHENET_RESOURCE_DIR:-${RESOURCE_DIR}/nichenet}"
ZENODO_BASE_URL="${NICHENET_ZENODO_BASE_URL:-https://zenodo.org/record/7074291/files}"

ensure_dir "${NICHENET_RESOURCE_DIR}"

download_one() {
  local file_name="$1"
  local dest="${NICHENET_RESOURCE_DIR}/${file_name}"
  local url="${ZENODO_BASE_URL}/${file_name}?download=1"

  if [[ -s "${dest}" ]]; then
    echo "已存在，跳过: ${dest}"
    return 0
  fi

  echo "下载 ${file_name}"
  if command -v curl >/dev/null 2>&1; then
    curl -L --retry 3 --retry-delay 5 -o "${dest}.tmp" "${url}"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "${dest}.tmp" "${url}"
  else
    die "缺少 curl/wget，无法下载 NicheNet 资源。"
  fi
  mv -f "${dest}.tmp" "${dest}"
}

download_one "lr_network_human_21122021.rds"
download_one "ligand_target_matrix_nsga2r_final.rds"
download_one "weighted_networks_nsga2r_final.rds"

cat > "${NICHENET_RESOURCE_DIR}/README.md" <<EOF
# NicheNet Resources

Downloaded from Zenodo record 7074291:
${ZENODO_BASE_URL}

Files:
- lr_network_human_21122021.rds
- ligand_target_matrix_nsga2r_final.rds
- weighted_networks_nsga2r_final.rds
EOF

echo "NicheNet 资源准备完成: ${NICHENET_RESOURCE_DIR}"
