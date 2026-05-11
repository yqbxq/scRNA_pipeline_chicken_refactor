detect_container_runtime() {
  if command -v apptainer >/dev/null 2>&1; then
    echo "apptainer"
    return
  fi
  if command -v singularity >/dev/null 2>&1; then
    echo "singularity"
    return
  fi
  echo ""
}

detect_conda_frontend() {
  if command -v micromamba >/dev/null 2>&1; then
    echo "micromamba"
    return
  fi
  if command -v mamba >/dev/null 2>&1; then
    echo "mamba"
    return
  fi
  if command -v conda >/dev/null 2>&1; then
    echo "conda"
    return
  fi
  echo ""
}

if ! declare -F append_env_vars >/dev/null 2>&1; then
  # shellcheck disable=SC1090
  source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/env_registry.sh"
fi

CONTAINER_RUNTIME="${CONTAINER_RUNTIME:-$(detect_container_runtime)}"
CONDA_FRONTEND="${CONDA_FRONTEND:-$(detect_conda_frontend)}"

CONTAINER_BIND_ARGS=()
BIND_PATHS=(
  "${PROJECT_ROOT}"
  "${PAPER_MATRIX_SOURCE:-}"
  "${CELLRANGER_OUT_DIR:-}"
  "${DNBC4TOOLS_OUT_DIR:-}"
  "${ORTHOLOG_CACHE_DIR:-}"
  "${DECOUPLER_RESOURCE_DIR:-}"
  "${REFERENCE_DIR:-}"
  "${ENV_DIR:-}"
)

if [[ -n "${EXISTING_R_SIF:-}" ]]; then
  BIND_PATHS+=("$(dirname "${EXISTING_R_SIF}")")
fi

for bind_path in "${BIND_PATHS[@]}"; do
  [[ -n "${bind_path}" ]] || continue
  if [[ -e "${bind_path}" ]]; then
    CONTAINER_BIND_ARGS+=(--bind "${bind_path}:${bind_path}")
  fi
done

run_r_main() {
  local script_path="$1"
  shift || true
  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_R_MAIN[@]}"

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env "${ENV_ARGS[@]}" Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行主流程 R 环境。"
  [[ -d "${R_MAIN_ENV_PREFIX}" ]] || die "主流程 R 环境不存在: ${R_MAIN_ENV_PREFIX}"

  env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${R_MAIN_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_r_scenic() {
  local script_path="$1"
  shift || true
  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_R_SCENIC[@]}"

  if [[ "${USE_EXISTING_SIF}" == "yes" ]]; then
    [[ -n "${CONTAINER_RUNTIME}" ]] || die "USE_EXISTING_SIF=yes，但系统中没有 apptainer/singularity。"
    [[ -f "${EXISTING_R_SIF}" ]] || die "找不到旧 R 容器: ${EXISTING_R_SIF}"

    "${CONTAINER_RUNTIME}" exec "${CONTAINER_BIND_ARGS[@]}" "${EXISTING_R_SIF}" \
      env "${ENV_ARGS[@]}" Rscript "${script_path}" "$@"
    return
  fi

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 SCENIC R 环境。"
  [[ -d "${R_SCENIC_ENV_PREFIX}" ]] || die "SCENIC R 环境不存在: ${R_SCENIC_ENV_PREFIX}"

  env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${R_SCENIC_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_in_conda_prefix() {
  local prefix="$1"
  shift
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"
  "${CONDA_FRONTEND}" run -p "${prefix}" "$@"
}

create_or_update_conda_env() {
  local prefix="$1"
  local yaml_path="$2"
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"

  if [[ -d "${prefix}" ]]; then
    "${CONDA_FRONTEND}" env update -y -p "${prefix}" -f "${yaml_path}" --prune
  else
    "${CONDA_FRONTEND}" env create -y -p "${prefix}" -f "${yaml_path}"
  fi
}

run_pyscenic() {
  [[ -d "${PYSCENIC_ENV_PREFIX}" ]] || die "pySCENIC 环境不存在: ${PYSCENIC_ENV_PREFIX}"
  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_PYSCENIC[@]}"
  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda。"
  env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${PYSCENIC_ENV_PREFIX}" "$@"
}

run_velocity() {
  local command_path="${1:-}"
  [[ -n "${command_path}" ]] || die "run_velocity 需要命令或脚本路径。"
  shift || true

  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_VELOCITY[@]}"

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 velocity 环境。"
  [[ -d "${VELOCITY_ENV_PREFIX}" ]] || die "velocity 环境不存在: ${VELOCITY_ENV_PREFIX}"

  case "${command_path}" in
    *.py)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${VELOCITY_ENV_PREFIX}" python "${command_path}" "$@"
      ;;
    *.R)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${VELOCITY_ENV_PREFIX}" Rscript "${command_path}" "$@"
      ;;
    *.sh)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${VELOCITY_ENV_PREFIX}" bash "${command_path}" "$@"
      ;;
    *)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${VELOCITY_ENV_PREFIX}" "${command_path}" "$@"
      ;;
  esac
}

run_scvelo() {
  local command_path="${1:-}"
  [[ -n "${command_path}" ]] || die "run_scvelo 需要命令或脚本路径。"
  shift || true

  local scvelo_prefix="${SCVELO_ENV_PREFIX:-}"
  if [[ -z "${scvelo_prefix}" || ! -d "${scvelo_prefix}" ]]; then
    scvelo_prefix="${VELOCITY_ENV_PREFIX}"
  fi

  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_SCVELO[@]}"

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 scVelo 环境。"
  [[ -d "${scvelo_prefix}" ]] || die "scVelo/velocity 环境不存在: ${SCVELO_ENV_PREFIX:-<unset>} 或 ${VELOCITY_ENV_PREFIX}"

  case "${command_path}" in
    *.py)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${scvelo_prefix}" python "${command_path}" "$@"
      ;;
    *.R)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${scvelo_prefix}" Rscript "${command_path}" "$@"
      ;;
    *.sh)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${scvelo_prefix}" bash "${command_path}" "$@"
      ;;
    *)
      env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${scvelo_prefix}" "${command_path}" "$@"
      ;;
  esac
}

run_r_legacy() {
  [[ -d "${R_LEGACY_ENV_PREFIX}" ]] || die "legacy R 环境不存在: ${R_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${R_LEGACY_ENV_PREFIX}" Rscript "$@"
}

run_r_interaction() {
  local script_path="$1"
  shift || true
  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_R_INTERACTION[@]}"

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行通讯分析 R 环境。"
  [[ -d "${R_INTERACTION_ENV_PREFIX}" ]] || die "通讯分析 R 环境不存在: ${R_INTERACTION_ENV_PREFIX}"
  env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${R_INTERACTION_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_r_decoupler() {
  local script_path="$1"
  shift || true
  local -a ENV_ARGS=()
  append_env_vars "${ENV_SPECS_R_DECOUPLER[@]}"

  [[ -n "${CONDA_FRONTEND}" ]] || die "系统中找不到 micromamba/mamba/conda，无法运行 decoupleR 环境。"
  [[ -d "${R_DECOUPLER_ENV_PREFIX}" ]] || die "decoupleR R 环境不存在: ${R_DECOUPLER_ENV_PREFIX}"
  env "${ENV_ARGS[@]}" "${CONDA_FRONTEND}" run -p "${R_DECOUPLER_ENV_PREFIX}" Rscript "${script_path}" "$@"
}

run_r_spatial() {
  [[ -d "${R_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 R 环境不存在: ${R_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${R_SPATIAL_ENV_PREFIX}" Rscript "$@"
}

run_py_spatial() {
  [[ -d "${PY_SPATIAL_ENV_PREFIX}" ]] || die "空间分析 Python 环境不存在: ${PY_SPATIAL_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_ENV_PREFIX}" "$@"
}

run_py_spatial_legacy() {
  [[ -d "${PY_SPATIAL_LEGACY_ENV_PREFIX}" ]] || die "空间 legacy Python 环境不存在: ${PY_SPATIAL_LEGACY_ENV_PREFIX}"
  run_in_conda_prefix "${PY_SPATIAL_LEGACY_ENV_PREFIX}" "$@"
}

run_py_cell2location() {
  [[ -d "${PY_CELL2LOCATION_ENV_PREFIX}" ]] || die "cell2location Python 环境不存在: ${PY_CELL2LOCATION_ENV_PREFIX}"
  run_in_conda_prefix "${PY_CELL2LOCATION_ENV_PREFIX}" "$@"
}

run_r_validation() {
  [[ -d "${R_VALIDATION_ENV_PREFIX}" ]] || die "验证分析 R 环境不存在: ${R_VALIDATION_ENV_PREFIX}"
  run_in_conda_prefix "${R_VALIDATION_ENV_PREFIX}" Rscript "$@"
}
