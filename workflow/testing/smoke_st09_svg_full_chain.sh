#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

bash -n workflow/03stages/60_run_spatial_extensions.sh
grep -q '09a_spatialde2_svg.R' workflow/03stages/60_run_spatial_extensions.sh
grep -q '09b_sparkx_svg.R' workflow/03stages/60_run_spatial_extensions.sh
grep -q '09c_svg_consensus_report.R' workflow/03stages/60_run_spatial_extensions.sh
grep -q 'SVG_REQUIRE_REVIEW' workflow/03stages/60_run_spatial_extensions.sh

echo "smoke_st09_svg_full_chain_ok"
