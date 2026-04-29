detect_project_input_mode() {
  local file_path="${CANONICAL_SAMPLE_SHEET}"
  if [[ ! -f "${file_path}" ]]; then
    file_path="${SAMPLE_SHEET}"
  fi
  if [[ ! -f "${file_path}" ]]; then
    echo ""
    return
  fi

  local col_idx
  if ! col_idx="$(tsv_get_col_index "${file_path}" "input_mode" 2>/dev/null)"; then
    echo ""
    return
  fi

  awk -F '\t' -v col="${col_idx}" '
    NR > 1 && $col != "" {
      seen[$col] = 1
    }
    END {
      count = 0
      for (value in seen) {
        values[++count] = value
      }
      if (count == 1) {
        print values[1]
      } else if (count > 1) {
        print "mixed"
      }
    }
  ' "${file_path}"
}

has_fastq_files() {
  local dir_path="$1"
  [[ -d "${dir_path}" ]] || return 1
  find "${dir_path}" -maxdepth 1 -type f \
    \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) \
    -print -quit | grep -q .
}

is_10x_matrix_dir() {
  local dir_path="$1"
  [[ -d "${dir_path}" ]] || return 1
  [[ -f "${dir_path}/matrix.mtx.gz" || -f "${dir_path}/matrix.mtx" ]] || return 1
  [[ -f "${dir_path}/features.tsv.gz" || -f "${dir_path}/features.tsv" || -f "${dir_path}/genes.tsv.gz" || -f "${dir_path}/genes.tsv" ]] || return 1
  [[ -f "${dir_path}/barcodes.tsv.gz" || -f "${dir_path}/barcodes.tsv" ]] || return 1
}

resolve_cellranger_sample_dir() {
  local sample_id="$1"
  local input_mode="$2"
  local source_path="$3"

  case "${input_mode}" in
    fastq)
      if [[ -d "${CELLRANGER_OUT_DIR}/${sample_id}/outs" ]]; then
        echo "${CELLRANGER_OUT_DIR}/${sample_id}"
      fi
      ;;
    cellranger_out)
      if [[ -d "${source_path}/outs" ]]; then
        echo "${source_path}"
      elif [[ -d "${source_path}/${sample_id}/outs" ]]; then
        echo "${source_path}/${sample_id}"
      elif [[ -d "${CELLRANGER_OUT_DIR}/${sample_id}/outs" ]]; then
        echo "${CELLRANGER_OUT_DIR}/${sample_id}"
      fi
      ;;
    *)
      ;;
  esac
}

resolve_dnbelab_sample_dir() {
  local sample_id="$1"
  local source_path="$2"
  local candidate

  for candidate in \
    "${DNBC4TOOLS_OUT_DIR}/${sample_id}" \
    "${source_path}" \
    "${source_path}/${sample_id}" \
    "${DATA_DIR}/${sample_id}"; do
    [[ -n "${candidate}" ]] || continue
    if [[ -d "${candidate}/filter_matrix" || -d "${candidate}/raw_matrix" ]]; then
      echo "${candidate}"
      return 0
    fi
  done
}

resolve_matrix_source_dir() {
  local sample_id="$1"
  local input_mode="$2"
  local source_path="$3"
  local sample_dir=""
  local dnbelab_dir=""

  dnbelab_dir="$(resolve_dnbelab_sample_dir "${sample_id}" "${source_path}")"
  if [[ -n "${dnbelab_dir}" ]] && is_10x_matrix_dir "${dnbelab_dir}/filter_matrix"; then
    echo "${dnbelab_dir}/filter_matrix"
    return 0
  fi

  case "${input_mode}" in
    fastq|cellranger_out)
      sample_dir="$(resolve_cellranger_sample_dir "${sample_id}" "${input_mode}" "${source_path}")"
      if [[ -n "${sample_dir}" ]]; then
        echo "${sample_dir}/outs/filtered_feature_bc_matrix"
      fi
      ;;
    matrix)
      if is_10x_matrix_dir "${source_path}"; then
        echo "${source_path}"
      elif is_10x_matrix_dir "${source_path}/filtered_feature_bc_matrix"; then
        echo "${source_path}/filtered_feature_bc_matrix"
      elif is_10x_matrix_dir "${source_path}/${sample_id}"; then
        echo "${source_path}/${sample_id}"
      elif is_10x_matrix_dir "${DATA_DIR}/${sample_id}"; then
        echo "${DATA_DIR}/${sample_id}"
      fi
      ;;
    *)
      ;;
  esac
}
