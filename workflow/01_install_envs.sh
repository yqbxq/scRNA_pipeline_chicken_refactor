#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

INSTALL_TARGETS="${INSTALL_TARGETS:-r_main,r_legacy,r_scenic,velocity,pyscenic}"
INSTALL_SCENIC_RESOURCES="${INSTALL_SCENIC_RESOURCES:-yes}"

ensure_dir \
  "${PROJECT_ROOT}" \
  "${ENV_DIR}" \
  "${LOG_DIR}" \
  "${RESULTS_DIR}" \
  "${CHECKPOINT_DIR}" \
  "${FIGURE_DIR}" \
  "${TABLE_DIR}" \
  "${VELOCITY_DIR}" \
  "${VELOCITY_INPUT_DIR}" \
  "${VELOCITY_LOOM_DIR}" \
  "${VELOCITY_OUTPUT_DIR}" \
  "${SCENIC_INPUT_DIR}" \
  "${SCENIC_OUTPUT_DIR}" \
  "${R_LIBS_MAIN}" \
  "${R_LIBS_SCENIC}"

run_with_log() {
  local log_file="$1"
  shift
  "$@" 2>&1 | tee "${log_file}"
}

normalize_target_name() {
  case "$1" in
    py_main)
      echo "velocity"
      ;;
    py_scenic)
      echo "pyscenic"
      ;;
    r_heart_legacy)
      echo "r_legacy"
      ;;
    *)
      echo "$1"
      ;;
  esac
}

target_enabled() {
  local requested="${INSTALL_TARGETS// /}"
  local needle
  needle="$(normalize_target_name "$1")"
  [[ ",${requested}," == *",all,"* ]] || [[ ",${requested}," == *",${needle},"* ]]
}

echo "安装目标: ${INSTALL_TARGETS}"
echo "SCENIC 资源安装: ${INSTALL_SCENIC_RESOURCES}"
echo "说明: 现有 prefix 会原地 update，不会新建重复环境；只有不存在的 prefix 才会 create。"

if target_enabled "r_main"; then
  run_with_log "${LOG_DIR}/install_r_main_env.log" \
    create_or_update_conda_env "${R_MAIN_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_main.yml"
  run_with_log "${LOG_DIR}/install_r_main_pkgs.log" \
    run_r_main "${PIPELINE_ROOT}/envs/install_r_main_packages.R"
fi

if target_enabled "r_scenic"; then
  run_with_log "${LOG_DIR}/install_r_scenic_env.log" \
    create_or_update_conda_env "${R_SCENIC_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_scenic.yml"
  run_with_log "${LOG_DIR}/install_r_scenic_pkgs.log" \
    run_r_scenic "${PIPELINE_ROOT}/envs/install_r_scenic_packages.R"
fi

if target_enabled "r_legacy"; then
  run_with_log "${LOG_DIR}/install_r_legacy_env.log" \
    create_or_update_conda_env "${R_LEGACY_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_legacy.yml"
  run_with_log "${LOG_DIR}/install_r_legacy_pkgs.log" \
    run_r_legacy "${PIPELINE_ROOT}/envs/install_r_legacy_packages.R"
fi

if target_enabled "velocity"; then
  run_with_log "${LOG_DIR}/install_velocity_env.log" \
    create_or_update_conda_env "${VELOCITY_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_velocity.yml"
fi

if target_enabled "pyscenic"; then
  run_with_log "${LOG_DIR}/install_pyscenic_env.log" \
    create_or_update_conda_env "${PYSCENIC_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_pyscenic.yml"
fi

if target_enabled "r_interaction"; then
  run_with_log "${LOG_DIR}/install_r_interaction_env.log" \
    create_or_update_conda_env "${R_INTERACTION_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_interaction.yml"
  run_with_log "${LOG_DIR}/install_r_interaction_pkgs.log" \
    run_in_conda_prefix "${R_INTERACTION_ENV_PREFIX}" Rscript "${PIPELINE_ROOT}/envs/install_r_interaction_packages.R"
fi

if target_enabled "r_spatial"; then
  run_with_log "${LOG_DIR}/install_r_spatial_env.log" \
    create_or_update_conda_env "${R_SPATIAL_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_spatial.yml"
  run_with_log "${LOG_DIR}/install_r_spatial_pkgs.log" \
    run_in_conda_prefix "${R_SPATIAL_ENV_PREFIX}" Rscript "${PIPELINE_ROOT}/envs/install_r_spatial_packages.R"
fi

if target_enabled "py_spatial"; then
  run_with_log "${LOG_DIR}/install_py_spatial_env.log" \
    create_or_update_conda_env "${PY_SPATIAL_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_py_spatial.yml"
  run_with_log "${LOG_DIR}/install_py_spatial_pkgs.log" \
    run_in_conda_prefix "${PY_SPATIAL_ENV_PREFIX}" bash "${PIPELINE_ROOT}/envs/install_py_spatial_packages.sh"
fi

if target_enabled "py_spatial_legacy"; then
  run_with_log "${LOG_DIR}/install_py_spatial_legacy_env.log" \
    create_or_update_conda_env "${PY_SPATIAL_LEGACY_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_py_spatial_legacy.yml"
  run_with_log "${LOG_DIR}/install_py_spatial_legacy_pkgs.log" \
    run_in_conda_prefix "${PY_SPATIAL_LEGACY_ENV_PREFIX}" bash "${PIPELINE_ROOT}/envs/install_py_spatial_legacy_packages.sh"
fi

if target_enabled "py_cell2location"; then
  run_with_log "${LOG_DIR}/install_py_cell2location_env.log" \
    create_or_update_conda_env "${PY_CELL2LOCATION_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_py_cell2location.yml"
  run_with_log "${LOG_DIR}/install_py_cell2location_pkgs.log" \
    run_in_conda_prefix "${PY_CELL2LOCATION_ENV_PREFIX}" bash "${PIPELINE_ROOT}/envs/install_py_cell2location_packages.sh"
fi

if target_enabled "r_validation"; then
  run_with_log "${LOG_DIR}/install_r_validation_env.log" \
    create_or_update_conda_env "${R_VALIDATION_ENV_PREFIX}" "${PIPELINE_ROOT}/envs/environment_r_validation.yml"
  run_with_log "${LOG_DIR}/install_r_validation_pkgs.log" \
    run_in_conda_prefix "${R_VALIDATION_ENV_PREFIX}" Rscript "${PIPELINE_ROOT}/envs/install_r_validation_packages.R"
fi

if [[ "${INSTALL_SCENIC_RESOURCES}" == "yes" ]] && { target_enabled "r_scenic" || target_enabled "pyscenic"; }; then
  run_with_log "${LOG_DIR}/install_scenic_resources.log" \
    "${PIPELINE_ROOT}/workflow/41_download_scenic_resources.sh"
fi

echo "服务器环境安装完成。"
echo "安装日志目录: ${LOG_DIR}"
