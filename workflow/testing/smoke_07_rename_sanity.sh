#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

fail() {
  echo "[FAIL] $*" >&2
  exit 1
}

require_file() {
  [[ -f "$1" ]] || fail "missing file: $1"
}

reject_file() {
  [[ ! -e "$1" ]] || fail "old file should be absent: $1"
}

require_grep() {
  local pattern="$1"
  local path="$2"
  rg -q "${pattern}" "${path}" || fail "missing pattern ${pattern} in ${path}"
}

reject_grep() {
  local pattern="$1"
  shift
  if rg -q "${pattern}" "$@"; then
    fail "unexpected legacy pattern: ${pattern}"
  fi
}

require_file workflow/04python/07b_liana_consensus.py
require_file workflow/05single_script/07c_nichenet.R
require_file workflow/05single_script/07d_communication_consensus.R
require_file workflow/05single_script/07e_communication_eda.R
require_file workflow/05single_script/helpers/communication_consensus_utils.R
old_nichenet="workflow/05single_script/07b_""nichenet.R"
old_eda="workflow/05single_script/07c_""communication_eda.R"
reject_file "${old_nichenet}"
reject_file "${old_eda}"

require_grep "07b_liana_consensus/_manifest.json" workflow/03stages/07_communication.sh
require_grep "07c_nichenet/_manifest.json" workflow/03stages/07_communication.sh
require_grep "07d_communication_consensus/_manifest.json" workflow/03stages/07_communication.sh
require_grep "07e_communication_eda/_manifest.json" workflow/03stages/07_communication.sh
require_grep "COMMUNICATION_REQUIRE_FULL_PIPELINE" workflow/03stages/07_communication.sh
require_grep "07e_communication" workflow/02lib/common.sh

python3 -m py_compile workflow/04python/07b_liana_consensus.py
if python3 workflow/04python/07b_liana_consensus.py >/tmp/smoke_07b.out 2>&1; then
  fail "07b placeholder should exit non-zero"
fi
rg -q "not yet implemented" /tmp/smoke_07b.out || fail "07b placeholder message missing"

Rscript -e 'parse("workflow/05single_script/07d_communication_consensus.R"); parse("workflow/05single_script/helpers/communication_consensus_utils.R")' >/dev/null
if Rscript workflow/05single_script/07d_communication_consensus.R >/tmp/smoke_07d.out 2>&1; then
  fail "07d placeholder should exit non-zero"
fi
rg -q "not yet implemented" /tmp/smoke_07d.out || fail "07d placeholder message missing"

legacy_pattern="07b_""nichenet\\.R|07c_""communication_eda\\.R|07b_""nichenet/_manifest|07c_""communication_eda/_manifest|MODULE_07B_""VERSION|MODULE_07C_""VERSION"
reject_grep "${legacy_pattern}" \
  workflow docs README.md config SOP_标准流程与论文模式.md

echo "smoke_07_rename_sanity_ok"
