#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p "${TMP_DIR}/metadata" "${TMP_DIR}/results/tables/communication/consensus"
cat > "${TMP_DIR}/metadata/sections.tsv" <<'EOF'
section_id	chip_id	tissue_block	condition	platform	bundle_root	image_lowres	image_hires	tissue_positions	scalefactors	enabled	commot_enabled	notes
sec1	chip1	block1	syf	visium	bundle	lowres	hires	positions	scalefactors	yes	auto	
EOF
cat > "${TMP_DIR}/results/tables/communication/consensus/method_consensus.tsv" <<'EOF'
lr_axis_id	condition_value	pair_id	layer_id	source	target	ligand	receptor	ligand_human	receptor_human	cellchat_hit	liana_consensus_hit	nichenet_hit	commot_spatial_hit	evidence_method_n	evidence_tier	can_be_primary
L1|R1|A->B	syf	P1	panorama	A	B	L1	R1	L1	R1	yes	yes	yes	no	3	primary	yes
EOF

PROJECT_ROOT="${TMP_DIR}" \
PIPELINE_ROOT="${REPO_ROOT}" \
METADATA_DIR="${TMP_DIR}/metadata" \
RESULTS_DIR="${TMP_DIR}/results" \
COMMOT_ENABLED="auto" \
COMMOT_REQUIRE_RUNTIME="no" \
Rscript "${REPO_ROOT}/workflow/05single_script/spatial/08a_spatial_communication_io.R" >/tmp/smoke_08a_commot.log

test -s "${TMP_DIR}/results/spatial/tables/08_commot/commot_input_manifest.tsv"
test -s "${TMP_DIR}/results/spatial/tables/08_commot/commot_lr_candidates.tsv"

PROJECT_ROOT="${TMP_DIR}" \
RESULTS_DIR="${TMP_DIR}/results" \
COMMOT_REQUIRE_RUNTIME="no" \
python3 "${REPO_ROOT}/workflow/04python/07f_commot_spatial.py" >/tmp/smoke_08f_commot.log

test -s "${TMP_DIR}/results/spatial/tables/08_commot/commot_lr.tsv"

PROJECT_ROOT="${TMP_DIR}" \
PIPELINE_ROOT="${REPO_ROOT}" \
METADATA_DIR="${TMP_DIR}/metadata" \
RESULTS_DIR="${TMP_DIR}/results" \
Rscript "${REPO_ROOT}/workflow/05single_script/spatial/08b_spatial_communication_eda.R" >/tmp/smoke_08b_commot.log

test -s "${TMP_DIR}/results/spatial/tables/08_commot/commot_spatial_summary.tsv"
test -s "${TMP_DIR}/reports/eda/spatial_08_commot/report.md"
grep -q 'skipped_' "${TMP_DIR}/results/spatial/tables/08_commot/commot_spatial_summary.tsv"
echo "smoke_08_commot_pipeline: ok"
