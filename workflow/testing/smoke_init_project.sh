#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
TMP_ROOT="$(mktemp -d)"

cleanup() {
  rm -rf "${TMP_ROOT}"
}
trap cleanup EXIT

mkdir -p "${TMP_ROOT}/refs" "${TMP_ROOT}/matrix/S1" "${TMP_ROOT}/matrix/S2"

cat > "${TMP_ROOT}/refs/genome.fa" <<'EOF'
>chr1
ACGT
EOF

cat > "${TMP_ROOT}/refs/genes.gtf" <<'EOF'
chr1	source	gene	1	4	.	+	.	gene_id "g1"; gene_name "G1";
EOF

for sample_id in S1 S2; do
  cat > "${TMP_ROOT}/matrix/${sample_id}/matrix.mtx" <<'EOF'
%%MatrixMarket matrix coordinate integer general
1 1 1
1 1 1
EOF
  cat > "${TMP_ROOT}/matrix/${sample_id}/features.tsv" <<'EOF'
g1	G1	Gene Expression
EOF
  printf '%s_cell1\n' "${sample_id}" > "${TMP_ROOT}/matrix/${sample_id}/barcodes.tsv"
done

"${PIPELINE_ROOT}/workflow/03stages/init_project.sh" \
  --project-root "${TMP_ROOT}/project" \
  --genome-fasta "${TMP_ROOT}/refs/genome.fa" \
  --reference-gtf "${TMP_ROOT}/refs/genes.gtf" \
  --sample-names S1,S2 \
  --group1-name syf \
  --group1-samples S1 \
  --group2-name f5 \
  --group2-samples S2 \
  --external-matrix-source "${TMP_ROOT}/matrix" >/dev/null

export SCRNA_PIPELINE_CONFIG="${TMP_ROOT}/project/config/project_config.sh"

bash "${PIPELINE_ROOT}/workflow/03stages/validate_metadata.sh" >/dev/null
bash "${PIPELINE_ROOT}/workflow/03stages/audit_inputs.sh" >/dev/null
bash "${PIPELINE_ROOT}/workflow/03stages/standardize_inputs.sh" >/dev/null
bash "${PIPELINE_ROOT}/workflow/03stages/input_summary.sh" >/dev/null

[[ -s "${TMP_ROOT}/project/config/project_config.sh" ]]
grep -q '^export RBC_GENE_LIST_FILE=' "${TMP_ROOT}/project/config/project_config.sh"
grep -q 'RBC / hemoglobin marker gene list' "${TMP_ROOT}/project/config/rbc_gene_list.txt"
[[ -s "${TMP_ROOT}/project/metadata/samples.canonical.tsv" ]]
[[ -s "${TMP_ROOT}/project/reports/intake/input_inventory.tsv" ]]
[[ -s "${TMP_ROOT}/project/reports/intake/branch_readiness.tsv" ]]
[[ -s "${TMP_ROOT}/project/reports/intake/intake_summary.md" ]]
grep -q $'^trajectory_methods\tpending' "${TMP_ROOT}/project/config/eda_gates.tsv"
grep -q $'^scdesign3_targets\tpending' "${TMP_ROOT}/project/config/eda_gates.tsv"
grep -q $'^panorama\tpanorama\tyes' "${TMP_ROOT}/project/config/object_layers.tsv"

echo "smoke_init_project_ok"
