#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

Rscript -e 'parse("workflow/05single_script/spatial/06e_niche_derivation.R")'
grep -q 'spatial_06e_niche_h5ad' workflow/05single_script/spatial/06e_niche_derivation.R
grep -q '90a_export_h5ad.*spatial_06e_niche' workflow/05single_script/spatial/06e_niche_derivation.R
grep -q '^obs.niche_label' metadata/h5ad_export_contract.tsv

echo "smoke_spatial_niche_h5ad_export_ok"
