#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT

export PROJECT_ROOT="${TMP_DIR}"
export PIPELINE_ROOT="${ROOT}"
export RESULTS_DIR="${TMP_DIR}/results"
export REPORT_DIR="${TMP_DIR}/reports"
export METADATA_DIR="${TMP_DIR}/metadata"
export PROJECT_CONFIG_DIR="${TMP_DIR}/config"
export MANIFEST_DIR="${TMP_DIR}/results/manifests"
export SPATIAL_TABLE_DIR="${TMP_DIR}/results/spatial/tables"
export SPATIAL_FIGURE_DIR="${TMP_DIR}/results/spatial/figures"
export SPATIAL_COMMUNICATION_TABLE_DIR="${SPATIAL_TABLE_DIR}/08_spatial_communication"
export SPATIAL_COMMUNICATION_REPORT_DIR="${TMP_DIR}/reports/eda/spatial_08_communication"
export COMMOT_LR_CANDIDATES_TSV="${TMP_DIR}/results/tables/communication/consensus/method_consensus.tsv"
export SPATIAL_VALIDATION_TRUTH_TSV="${TMP_DIR}/truth.tsv"
export SPATIAL_VALIDATION_MODE="external_truth"

mkdir -p \
  "${TMP_DIR}/metadata" \
  "${TMP_DIR}/results/tables/communication/consensus" \
  "${SPATIAL_TABLE_DIR}/spatial_07a_deconvolution_rctd/default_deconv/S1" \
  "${SPATIAL_TABLE_DIR}/spatial_07b_deconvolution_transfer/default_deconv/S1" \
  "${SPATIAL_TABLE_DIR}/spatial_06d_neighborhood" \
  "${SPATIAL_TABLE_DIR}/08_commot"

cat >"${TMP_DIR}/metadata/sections.tsv" <<'TSV'
section_id	chip_id	tissue_block	condition	platform	bundle_root	image_lowres	image_hires	tissue_positions	scalefactors	enabled	commot_enabled	notes
S1	C1	B1	syf	visium						yes	yes	
TSV

cat >"${COMMOT_LR_CANDIDATES_TSV}" <<'TSV'
lr_axis_id	ligand	receptor	source	target	pair_id	condition_value	evidence_tier	cellchat_hit	liana_consensus_hit	nichenet_hit	nichenet_n_targets_in_receiver_de	can_be_primary
L1|R1|Sender->Receiver	L1	R1	Sender	Receiver	P1	syf	primary	yes	yes	yes	5	yes
TSV

for method in rctd transfer; do
  case "${method}" in
    rctd) dir="${SPATIAL_TABLE_DIR}/spatial_07a_deconvolution_rctd"; manifest="rctd_manifest.tsv" ;;
    transfer) dir="${SPATIAL_TABLE_DIR}/spatial_07b_deconvolution_transfer"; manifest="transfer_manifest.tsv" ;;
  esac
  mkdir -p "${dir}/default_deconv/S1"
  cat >"${dir}/default_deconv/S1/spot_celltype_proportions.tsv" <<'TSV'
spot_id	cell_type	proportion
s1	Sender	0.8
s1	Receiver	0.7
s2	Sender	0.2
s2	Receiver	0.1
TSV
  cat >"${dir}/${manifest}" <<TSV
deconv_id	section	tool	status	reason	panorama_input_rds	n_spots	n_celltypes	runtime_sec	proportion_tsv	proportion_wide_tsv	spot_metadata_tsv	summary_tsv	method_object_rds	method_version
default_deconv	S1	${method}	ok			2	2	1	${dir}/default_deconv/S1/spot_celltype_proportions.tsv					1.0
TSV
done

cat >"${TMP_DIR}/truth.tsv" <<'TSV'
spot_id	cell_type	true_proportion
s1	Sender	0.75
s1	Receiver	0.70
s2	Sender	0.25
s2	Receiver	0.10
TSV

cat >"${SPATIAL_TABLE_DIR}/spatial_06d_neighborhood/interaction_matrix.tsv" <<'TSV'
label	Sender	Receiver
Sender	1	1
Receiver	1	1
TSV
cat >"${SPATIAL_TABLE_DIR}/spatial_06d_neighborhood/nhood_enrichment_zscore.tsv" <<'TSV'
label	Sender	Receiver
Sender	1	2
Receiver	2	1
TSV

Rscript "${ROOT}/workflow/05single_script/spatial/07e_deconvolution_compare.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/07f_deconvolution_validation.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/08a_spatial_communication_io.R" >/dev/null

cat >"${SPATIAL_TABLE_DIR}/08_commot/commot_lr.tsv" <<'TSV'
section_id	lr_axis_id	pair_id	condition_value	ligand	receptor	sender	receiver	signal_score	spatial_support	status	reason
S1	L1|R1|Sender->Receiver	P1	syf	L1	R1	Sender	Receiver	1.2	no	ok_commot_run	
TSV

Rscript "${ROOT}/workflow/05single_script/spatial/08b_spatial_communication_eda.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/08c_lr_colocalization.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/08d_communication_neighborhood_consistency.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/08e_spatial_communication_consensus.R" >/dev/null
Rscript "${ROOT}/workflow/05single_script/spatial/08f_spatial_communication_report.R" >/dev/null

test -s "${SPATIAL_TABLE_DIR}/spatial_07e_deconvolution_compare/method_summary.tsv"
test -s "${SPATIAL_TABLE_DIR}/spatial_07f_deconvolution_validation/spatial_question_gate_status.tsv"
test -s "${SPATIAL_COMMUNICATION_TABLE_DIR}/spatial_comm_candidates.tsv"
test -s "${SPATIAL_TABLE_DIR}/08_commot/commot_celltype_scores.tsv"
test -s "${SPATIAL_COMMUNICATION_TABLE_DIR}/lr_colocalization.tsv"
test -s "${SPATIAL_COMMUNICATION_TABLE_DIR}/communication_neighborhood_consistency.tsv"
test -s "${SPATIAL_COMMUNICATION_TABLE_DIR}/spatial_communication_consensus.tsv"
test -s "${SPATIAL_COMMUNICATION_TABLE_DIR}/spatial_communication_question_gate_status.tsv"
grep -q 'I19_comm_in_space' "${SPATIAL_COMMUNICATION_TABLE_DIR}/spatial_communication_question_gate_status.tsv"
grep -q 'spatial_primary' "${SPATIAL_COMMUNICATION_TABLE_DIR}/spatial_communication_consensus.tsv"

echo "smoke_st0708_full_chain_ok"
