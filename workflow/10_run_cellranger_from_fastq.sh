#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

[[ -x "${CELLRANGER_BIN}" ]] || die "Cell Ranger 不可执行: ${CELLRANGER_BIN}"
[[ -d "${FASTQ_DIR}" ]] || die "缺少 FASTQ 目录: ${FASTQ_DIR}"
[[ -f "${REFERENCE_GTF}" ]] || die "缺少参考 GTF: ${REFERENCE_GTF}"

ensure_dir "${REFERENCE_DIR}" "${CELLRANGER_OUT_DIR}"

if [[ ! -f "${CLEAN_GTF}" ]]; then
  awk -F '\t' -v OFS='\t' '/^#/ {print; next} {gsub(/ /, "_", $2); print}' "${REFERENCE_GTF}" > "${CLEAN_GTF}"
fi

if [[ ! -d "${CELLRANGER_REF_DIR}" ]]; then
  [[ -f "${GENOME_FASTA_GZ}" ]] || die "缺少基因组 FASTA 文件: ${GENOME_FASTA_GZ}"
  tmp_fa="${REFERENCE_DIR}/genome.fa"
  if [[ ! -f "${tmp_fa}" ]]; then
    case "${GENOME_FASTA_GZ}" in
      *.gz)
        gunzip -c "${GENOME_FASTA_GZ}" > "${tmp_fa}"
        ;;
      *)
        cp "${GENOME_FASTA_GZ}" "${tmp_fa}"
        ;;
    esac
  fi

  "${CELLRANGER_BIN}" mkref \
    --genome=GRCg7b \
    --fasta="${tmp_fa}" \
    --genes="${CLEAN_GTF}" \
    --nthreads="${CELLRANGER_THREADS}"
fi

IFS=',' read -r -a raw_samples <<< "${RAW_SAMPLES}"
for sample in "${raw_samples[@]}"; do
  if [[ -d "${CELLRANGER_OUT_DIR}/${sample}" ]]; then
    echo "已存在 ${sample}，跳过"
    continue
  fi

  "${CELLRANGER_BIN}" count \
    --id="${sample}" \
    --transcriptome="${CELLRANGER_REF_DIR}" \
    --fastqs="${FASTQ_DIR}" \
    --sample="${sample}" \
    --create-bam=true \
    --localcores="${CELLRANGER_THREADS}" \
    --localmem="${CELLRANGER_MEM_GB}"
done

update_workflow_status \
  "cellranger_completed" \
  "bash ${PIPELINE_ROOT}/workflow/03_audit_inputs.sh"

echo "Cell Ranger 运行完成"
