#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

mkdir -p \
  "${TMP_DIR}/metadata" \
  "${TMP_DIR}/config" \
  "${TMP_DIR}/results/spatial/tables/09_svg/spatialde2" \
  "${TMP_DIR}/results/spatial/tables/09_svg/sparkx" \
  "${TMP_DIR}/results/spatial/tables/08_spatial_communication"

cat >"${TMP_DIR}/results/spatial/tables/09_svg/spatialde2/all_svg.tsv" <<'TSV'
section_id	condition	method	gene_id	gene_symbol	score	pvalue	qvalue	rank	mean_expression	n_spots	n_genes_tested	status	reason	is_ligand_candidate	is_receptor_candidate	is_receiver_target_candidate	is_region_marker_candidate	is_regulator_target_candidate	linked_question_ids
all	syf	SpatialDE2	G001	LIG1	10	0.001	0.01	1	5	20	2	ok		no	no	no	no	no	I06_SVG
all	syf	SpatialDE2	G002	REC1	4	0.02	0.05	2	3	20	2	ok		no	no	no	no	no	I06_SVG
TSV
cat >"${TMP_DIR}/results/spatial/tables/09_svg/sparkx/all_svg.tsv" <<'TSV'
section_id	condition	method	gene_id	gene_symbol	score	pvalue	qvalue	rank	mean_expression	n_spots	n_genes_tested	status	reason	is_ligand_candidate	is_receptor_candidate	is_receiver_target_candidate	is_region_marker_candidate	is_regulator_target_candidate	linked_question_ids
all	syf	SPARK-X	G001	LIG1	9	0.001	0.01	1	5	20	2	ok		no	no	no	no	no	I06_SVG
TSV
cat >"${TMP_DIR}/results/spatial/tables/08_spatial_communication/spatial_comm_candidates.tsv" <<'TSV'
comm_candidate_id	source_question_id	lr_axis_id	ligand	receptor	sender_cell_type	receiver_cell_type	stage	section_id	region_id	niche_id	pair_id	condition_value	scrna_evidence_tier	liana_support	multinichenet_support	receiver_deg_support	cellchat_hypothesis	status	reason
LR1	I01	LIG1|REC1|A->B	LIG1	REC1	A	B	syf	all			P1	syf	candidate	no	no	no	yes	candidate
TSV
cat >"${TMP_DIR}/metadata/scenic_targets.tsv" <<'TSV'
target_id	source_question_id	regulator	target_gene	condition
T1	I02	TF1	LIG1	syf
TSV
cat >"${TMP_DIR}/metadata/gene_program_targets.tsv" <<'TSV'
comparison_id	source_question_id	gene_symbol	condition
TR1	H01	LIG1	syf
TSV
touch "${TMP_DIR}/metadata/trajectory_pairs.tsv"

export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"

Rscript "${ROOT}/workflow/05single_script/spatial/09c_svg_consensus_report.R" >/dev/null

OUT="${TMP_DIR}/results/spatial/tables/09_svg/09c_consensus"
test -s "${OUT}/svg_consensus.tsv"
test -s "${OUT}/svg_lr_support.tsv"
test -s "${OUT}/svg_gene_sets.tsv"
test -s "${OUT}/svg_question_gate_status.tsv"
grep -q 'high_confidence_svg' "${OUT}/svg_consensus.tsv"
grep -q 'strong' "${OUT}/svg_lr_support.tsv"
grep -q 'I06_SVG' "${OUT}/svg_question_gate_status.tsv"

echo "smoke_spatial_svg_consensus_global_ok"
