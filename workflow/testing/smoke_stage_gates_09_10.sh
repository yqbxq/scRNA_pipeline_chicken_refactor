#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"

fail() {
  echo "smoke_stage_gates_09_10 failed: $*" >&2
  exit 1
}

assert_contains() {
  local file="$1"
  local pattern="$2"
  grep -F -- "${pattern}" "${file}" >/dev/null || fail "${file} missing pattern: ${pattern}"
}

assert_not_contains() {
  local file="$1"
  local pattern="$2"
  if grep -F -- "${pattern}" "${file}" >/dev/null; then
    fail "${file} should not contain pattern: ${pattern}"
  fi
}

line_number() {
  local file="$1"
  local pattern="$2"
  local line
  line="$(grep -nF -- "${pattern}" "${file}" | head -n 1 | cut -d: -f1 || true)"
  [[ -n "${line}" ]] || fail "${file} missing pattern for order check: ${pattern}"
  printf '%s\n' "${line}"
}

assert_before() {
  local file="$1"
  local first="$2"
  local second="$3"
  local first_line
  local second_line
  first_line="$(line_number "${file}" "${first}")"
  second_line="$(line_number "${file}" "${second}")"
  (( first_line < second_line )) || fail "${file} order mismatch: ${first} should appear before ${second}"
}

for shell_file in \
  workflow/03stages/09_trajectory.sh \
  workflow/03stages/10_velocity.sh \
  workflow/03stages/10_velocity_finalize.sh; do
  bash -n "${shell_file}"
done

for gate_id in trajectory_inputs trajectory_methods trajectory_finalize velocity_inputs velocity_finalize; do
  assert_contains workflow/02lib/common.sh "${gate_id}"
  assert_contains workflow/02lib/stage.sh "status.${gate_id}_gate_passed"
done

assert_contains workflow/02lib/stage.sh "09_trajectory)"
assert_contains workflow/02lib/stage.sh "10_velocity)"
assert_contains workflow/02lib/stage.sh "10_velocity_finalize)"
assert_contains workflow/02lib/stage.sh "status.10_velocity_stage3_completed"
assert_contains workflow/02lib/stage.sh "\"10_velocity_completed\""

assert_before workflow/03stages/09_trajectory.sh "hold_for_gate trajectory_inputs" "09c_trajectory_root.R"
assert_before workflow/03stages/09_trajectory.sh "hold_for_gate trajectory_methods" "09h_trajectory_tradeseq.R"
assert_contains workflow/03stages/09_trajectory.sh "MODULE_09A_MANIFEST"
assert_contains workflow/03stages/09_trajectory.sh "MODULE_09M_MANIFEST"
assert_contains workflow/03stages/09_trajectory.sh "MODULE_09L_MANIFEST"
assert_contains workflow/03stages/09_trajectory.sh "MODULE_09G_MANIFEST"
assert_contains workflow/03stages/09_trajectory.sh "MODULE_09E2_MANIFEST"
assert_contains workflow/03stages/09_trajectory.sh "trajectory_has_monocle2_optin"
assert_contains workflow/03stages/09_trajectory.sh "\"trajectory_finalize\""
assert_contains workflow/03stages/09_trajectory.sh "status.trajectory_finalize_gate_passed"
assert_not_contains workflow/03stages/09_trajectory.sh "RUN_PALANTIR"
assert_not_contains workflow/05single_script/09g_trajectory_palantir.R "RUN_PALANTIR"

assert_before workflow/03stages/10_velocity.sh "hold_for_gate velocity_inputs" "10c_scvelo_dynamical.py"
assert_contains workflow/03stages/10_velocity.sh "10a_run_velocyto.sh"
assert_contains workflow/03stages/10_velocity.sh "10b_prepare_velocity_reference.R"
assert_contains workflow/03stages/10_velocity.sh "10f_cellrank_fate.py"
assert_contains workflow/03stages/10_velocity.sh "10g_velocity_consistency.R"
assert_contains workflow/03stages/10_velocity.sh "10h_velocity_root_terminal.R"
assert_contains workflow/03stages/10_velocity.sh "status.10_velocity_stage3_completed=true"
assert_contains workflow/03stages/10_velocity.sh "status.velocity_inputs_gate_passed"
assert_not_contains workflow/03stages/10_velocity.sh "RUN_CELLRANK"
assert_not_contains workflow/04python/10f_cellrank_fate.py "RUN_CELLRANK"

assert_contains workflow/03stages/10_velocity_finalize.sh "check_stage_deps \"10_velocity_finalize\""
assert_contains workflow/03stages/10_velocity_finalize.sh "10g_velocity_consistency.R"
assert_contains workflow/03stages/10_velocity_finalize.sh "10i_velocity_eda.R"
assert_contains workflow/03stages/10_velocity_finalize.sh "\"velocity_finalize\""
assert_contains workflow/03stages/10_velocity_finalize.sh "MODULE_10I_MANIFEST"
assert_contains workflow/03stages/10_velocity_finalize.sh "status.velocity_finalize_gate_passed"

assert_contains workflow/05single_script/helpers/metadata_fanout_trajectory.R "slingshot,monocle3,paga_dpt,tradeseq,palantir"
assert_contains workflow/05single_script/helpers/metadata_fanout_trajectory.R "scvelo_dynamical,scvelo_stochastic,velocyto,cellrank"
assert_contains metadata/README.md '+monocle2`, `-palantir`, or `-cellrank`'

echo "smoke_stage_gates_09_10 passed"
