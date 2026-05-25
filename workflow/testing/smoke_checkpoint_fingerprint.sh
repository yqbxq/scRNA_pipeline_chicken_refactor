#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d)"

cleanup() {
  rm -rf "${TMP_ROOT}"
}
trap cleanup EXIT

# shellcheck disable=SC1090
source "${PIPELINE_ROOT}/workflow/02lib/common.sh"

run_test_stage() {
  bash "$1"
}

write_test_script() {
  local script_path="$1"
  cat > "${script_path}" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
mkdir -p "$(dirname "${MANIFEST_PATH}")"
if [[ -s "${RUN_COUNT_FILE}" ]]; then
  count="$(cat "${RUN_COUNT_FILE}")"
else
  count=0
fi
count="$((count + 1))"
printf '%s\n' "${count}" > "${RUN_COUNT_FILE}"
python3 - "${MANIFEST_PATH}" "${OUTPUT_DIR}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

manifest_path = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)
result_path = output_dir / "result.tsv"
result_path.write_text("status\nok\n", encoding="utf-8")
manifest = {
    "module": "checkpoint_smoke",
    "version": "1.0",
    "timestamp": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    "base_dir": str(output_dir),
    "inputs": {"input": "input.tsv"},
    "outputs": {
        "result_tsv": {
            "path": "result.tsv",
            "type": "tsv",
            "produced_by": "checkpoint_smoke.sh",
            "row_semantics": "status",
        }
    },
    "depends_on": [],
}
manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
PY
EOF
  chmod +x "${script_path}"
}

setup_case() {
  local case_dir="$1"
  rm -rf "${case_dir}"
  mkdir -p "${case_dir}/out"
  printf 'value\nold\n' > "${case_dir}/input.tsv"
  write_test_script "${case_dir}/stage.sh"
  : > "${case_dir}/run_count"
}

run_count() {
  local case_dir="$1"
  if [[ -s "${case_dir}/run_count" ]]; then
    cat "${case_dir}/run_count"
  else
    echo 0
  fi
}

run_checkpoint() {
  local mode="$1"
  local case_dir="$2"
  local param_value="${3:-alpha}"
  local manifest_path="${case_dir}/out/_manifest.json"
  local before after
  before="$(run_count "${case_dir}")"
  CHECKPOINT_MODE="${mode}" \
    TEST_PARAM="${param_value}" \
    MANIFEST_PATH="${manifest_path}" \
    OUTPUT_DIR="${case_dir}/out" \
    RUN_COUNT_FILE="${case_dir}/run_count" \
    run_stage_if_stale_with_runner \
      run_test_stage \
      --script "${case_dir}/stage.sh" \
      --manifest "${manifest_path}" \
      --inputs "${case_dir}/input.tsv" \
      --params TEST_PARAM >/dev/null
  after="$(run_count "${case_dir}")"
  echo "$((after - before))"
}

assert_delta() {
  local label="$1"
  local expected="$2"
  local actual="$3"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "[FAIL] ${label}: expected delta ${expected}, got ${actual}" >&2
    exit 1
  fi
}

exercise_mode() {
  local mode="$1"
  local case_dir="${TMP_ROOT}/${mode}"
  setup_case "${case_dir}"

  assert_delta "${mode} C1 baseline" 1 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  assert_delta "${mode} C2 reproduce" 0 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"

  sleep 1
  printf 'new\n' >> "${case_dir}/input.tsv"
  assert_delta "${mode} C3 input content changed" 1 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"

  sleep 1
  touch -m "${case_dir}/input.tsv"
  if [[ "${mode}" == "mtime" ]]; then
    assert_delta "${mode} C4 input touched only" 1 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  else
    assert_delta "${mode} C4 input touched only" 0 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  fi

  printf 'value\npreserved_mtime\n' > "${case_dir}/input.tsv"
  touch -d '@1000000000' "${case_dir}/input.tsv"
  if [[ "${mode}" == "mtime" ]]; then
    assert_delta "${mode} C5 input replaced preserve mtime" 0 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  else
    assert_delta "${mode} C5 input replaced preserve mtime" 1 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  fi

  sleep 1
  printf '\n# script hash change\n' >> "${case_dir}/stage.sh"
  if [[ "${mode}" == "mtime" ]]; then
    assert_delta "${mode} C6 script changed" 0 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  else
    assert_delta "${mode} C6 script changed" 1 "$(run_checkpoint "${mode}" "${case_dir}" alpha)"
  fi

  if [[ "${mode}" == "mtime" ]]; then
    assert_delta "${mode} C7 params changed" 0 "$(run_checkpoint "${mode}" "${case_dir}" beta)"
  else
    assert_delta "${mode} C7 params changed" 1 "$(run_checkpoint "${mode}" "${case_dir}" beta)"
  fi
}

exercise_mode mtime
exercise_mode fingerprint

echo "smoke_checkpoint_fingerprint_ok"
