#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT_DIR}"

Rscript -e 'invisible(parse("workflow/05single_script/05d_deg_eda.R")); cat("05d_parse_ok\n")'
rg -q "Evidence Tier Distribution" workflow/05single_script/05d_deg_eda.R
rg -q "Downstream Coverage" workflow/05single_script/05d_deg_eda.R
rg -q "evidence_tier_summary_tsv" workflow/05single_script/05d_deg_eda.R
rg -q "downstream_coverage_tsv" workflow/05single_script/05d_deg_eda.R
rg -q "require_manifest_output.*evidence_tier_summary_tsv" workflow/03stages/05_deg.sh

echo "smoke_05d_eda_evidence_tier_ok"
