#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p \
  "${TMP_DIR}/metadata" \
  "${TMP_DIR}/config" \
  "${TMP_DIR}/results/spatial/tables/08_spatial_communication" \
  "${TMP_DIR}/results/spatial/tables/08_commot" \
  "${TMP_DIR}/results/spatial/tables/09_svg/09c_consensus"

cat >"${TMP_DIR}/results/spatial/tables/08_spatial_communication/spatial_comm_candidates.tsv" <<'TSV'
comm_candidate_id	source_question_id	lr_axis_id	ligand	receptor	sender_cell_type	receiver_cell_type	stage	section_id	region_id	niche_id	pair_id	condition_value	scrna_evidence_tier	liana_support	multinichenet_support	receiver_deg_support	cellchat_hypothesis	status	reason
LR1	I01	LIG1|REC1|A->B	LIG1	REC1	A	B	syf	all			P1	syf	candidate	no	no	no	yes	candidate
TSV
cat >"${TMP_DIR}/results/spatial/tables/09_svg/09c_consensus/svg_lr_support.tsv" <<'TSV'
comm_candidate_id	lr_axis_id	ligand	receptor	sender_cell_type	receiver_cell_type	condition	ligand_svg_tier	receptor_svg_tier	ligand_svg_rank	receptor_svg_rank	receiver_target_svg_support_n	receiver_target_svg_support_genes	svg_spatial_support_level	svg_support_reason
LR1	LIG1|REC1|A->B	LIG1	REC1	A	B	syf	high_confidence_svg	unsupported	1		0		strong	SVG gene-level support: strong
TSV

export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"

Rscript "${ROOT}/workflow/05single_script/spatial/08e_spatial_communication_consensus.R" >/dev/null

CONS="${TMP_DIR}/results/spatial/tables/08_spatial_communication/spatial_communication_consensus.tsv"
test -s "${CONS}"
grep -q 'svg_spatial_support_level' "${CONS}"
grep -q 'strong' "${CONS}"
grep -q 'spatial_hypothesis' "${CONS}"
! grep -q 'spatial_primary' "${CONS}"

echo "smoke_svg_communication_support_ok"
