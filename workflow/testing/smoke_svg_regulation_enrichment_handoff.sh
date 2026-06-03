#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${ROOT}"

python3 -m py_compile workflow/04python/spatial_decoupler.py
Rscript -e 'files <- c("workflow/05single_script/spatial/08a_spatial_decoupler.R","workflow/05single_script/spatial/06c_region_enrichment_eda.R","workflow/05single_script/spatial/helpers/spatial_svg_utils.R"); for (f in files) parse(f)'
grep -q 'svg_gene_sets' workflow/04python/spatial_decoupler.py
grep -q 'SPATIAL_REGULATION_SVG_GENE_SETS' workflow/05single_script/spatial/08a_spatial_decoupler.R
grep -q 'SVG_ENRICHMENT_HANDOFF_TSV' workflow/05single_script/spatial/06c_region_enrichment_eda.R

echo "smoke_svg_regulation_enrichment_handoff_ok"
