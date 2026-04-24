#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")/../02lib" && pwd)/common.sh"

usage() {
  cat <<'EOF'
用法：
  bash shell/03stages/register_delivery.sh \
    --source-path /path/to/delivery \
    [--delivery-id delivery_20260411] \
    [--mode symlink|copy]

说明：
  1. 将原始交付物登记到 ${RAW_RECEIVED_DIR}。
  2. 默认使用软链接，避免重复拷贝大文件。
  3. 会追加 delivery manifest，并刷新 workflow 状态。
EOF
}

SOURCE_PATH=""
DELIVERY_ID=""
REGISTER_MODE="symlink"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --source-path)
      SOURCE_PATH="$2"
      shift 2
      ;;
    --delivery-id)
      DELIVERY_ID="$2"
      shift 2
      ;;
    --mode)
      REGISTER_MODE="$2"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "未知参数: $1"
      ;;
  esac
done

[[ -n "${SOURCE_PATH}" ]] || die "必须提供 --source-path"
[[ -e "${SOURCE_PATH}" ]] || die "交付物路径不存在: ${SOURCE_PATH}"
[[ "${REGISTER_MODE}" == "symlink" || "${REGISTER_MODE}" == "copy" ]] || die "--mode 仅支持 symlink 或 copy"

prepare_project_state_dirs

SOURCE_PATH="$(cd "$(dirname "${SOURCE_PATH}")" && pwd)/$(basename "${SOURCE_PATH}")"
if [[ -z "${DELIVERY_ID}" ]]; then
  DELIVERY_ID="delivery_$(date +%Y%m%d_%H%M%S)"
fi

TARGET_PATH="${RAW_RECEIVED_DIR}/${DELIVERY_ID}"
[[ ! -e "${TARGET_PATH}" ]] || die "目标交付物 ID 已存在: ${TARGET_PATH}"

case "${REGISTER_MODE}" in
  symlink)
    ln -sfn "${SOURCE_PATH}" "${TARGET_PATH}"
    ;;
  copy)
    cp -a "${SOURCE_PATH}" "${TARGET_PATH}"
    ;;
esac

if [[ ! -f "${DELIVERY_MANIFEST}" ]]; then
  printf 'delivery_id\tsource_path\tlink_target\tregistered_at\tmode\n' > "${DELIVERY_MANIFEST}"
fi
if [[ ! -f "${RECEIVED_FILES_MANIFEST}" ]]; then
  printf 'delivery_id\trelative_path\tfile_type\tregistered_at\n' > "${RECEIVED_FILES_MANIFEST}"
fi

registered_at="$(now_iso)"
printf '%s\t%s\t%s\t%s\t%s\n' \
  "${DELIVERY_ID}" \
  "${SOURCE_PATH}" \
  "${TARGET_PATH}" \
  "${registered_at}" \
  "${REGISTER_MODE}" >> "${DELIVERY_MANIFEST}"

if [[ -d "${SOURCE_PATH}" ]]; then
  while IFS= read -r -d '' entry_path; do
    relative_path="${entry_path#${SOURCE_PATH}/}"
    file_type="file"
    if [[ -d "${entry_path}" ]]; then
      file_type="dir"
    fi
    printf '%s\t%s\t%s\t%s\n' \
      "${DELIVERY_ID}" \
      "${relative_path}" \
      "${file_type}" \
      "${registered_at}" >> "${RECEIVED_FILES_MANIFEST}"
  done < <(find "${SOURCE_PATH}" -mindepth 1 -print0)
else
  printf '%s\t%s\t%s\t%s\n' \
    "${DELIVERY_ID}" \
    "$(basename "${SOURCE_PATH}")" \
    "file" \
    "${registered_at}" >> "${RECEIVED_FILES_MANIFEST}"
fi

update_workflow_status \
  "delivery_registered" \
  "bash ${PIPELINE_ROOT}/shell/03stages/validate_metadata.sh"

echo "交付物已登记: ${DELIVERY_ID}"
echo "原始挂载位置: ${TARGET_PATH}"
