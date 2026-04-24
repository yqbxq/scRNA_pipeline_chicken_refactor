#!/usr/bin/env bash
set -euo pipefail

STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
SHELL_ROOT="$(cd "${STAGE_DIR}/.." && pwd)"
PIPELINE_ROOT="${PIPELINE_ROOT:-$(cd "${SHELL_ROOT}/.." && pwd)}"

source "${SHELL_ROOT}/02lib/common.sh"
source "${SHELL_ROOT}/02lib/executor.sh"

[[ -d "${FASTQ_DIR}" ]] || die "缺少 FASTQ 目录: ${FASTQ_DIR}"
[[ -f "${REFERENCE_GTF}" ]] || die "缺少参考 GTF: ${REFERENCE_GTF}"
[[ -f "${GENOME_FASTA_GZ}" ]] || die "缺少基因组 FASTA: ${GENOME_FASTA_GZ}"

ensure_dir "${REFERENCE_DIR}" "${DNBC4TOOLS_OUT_DIR:-${DATA_DIR}/dnbc4tools_out}"

if [[ ! -f "${CLEAN_GTF}" ]]; then
  awk -F '\t' -v OFS='\t' '/^#/ {print; next} {gsub(/ /, "_", $2); print}' \
    "${REFERENCE_GTF}" > "${CLEAN_GTF}"
fi

GENOME_FA="${REFERENCE_DIR}/genome.fa"
if [[ ! -f "${GENOME_FA}" ]]; then
  case "${GENOME_FASTA_GZ}" in
    *.gz)
      gunzip -c "${GENOME_FASTA_GZ}" > "${GENOME_FA}"
      ;;
    *)
      cp "${GENOME_FASTA_GZ}" "${GENOME_FA}"
      ;;
  esac
fi

if [[ ! -d "${STAR_INDEX_DIR}" ]]; then
  ensure_dir "${STAR_INDEX_DIR}"
  STAR --runMode genomeGenerate \
    --genomeDir "${STAR_INDEX_DIR}" \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${CLEAN_GTF}" \
    --runThreadN "${STAR_THREADS:-8}" \
    --genomeSAindexNbases "${STAR_SA_INDEX_NBASES:-13}" # Chicken GRCg7b is smaller than human, so 13 is the intended default here.
fi

IFS=',' read -r -a raw_samples <<< "${RAW_SAMPLES}"
for sample in "${raw_samples[@]}"; do
  output_dir="${DNBC4TOOLS_OUT_DIR:-${DATA_DIR}/dnbc4tools_out}/${sample}"
  if [[ -d "${output_dir}" ]]; then
    echo "已存在 ${sample}，跳过"
    continue
  fi

  dnbc4tools count \
    --name "${sample}" \
    --outdir "${output_dir}" \
    --starindex "${STAR_INDEX_DIR}" \
    --gtf "${CLEAN_GTF}" \
    --read1 "${FASTQ_DIR}/${sample}"*R1*.fastq.gz \
    --read2 "${FASTQ_DIR}/${sample}"*R2*.fastq.gz \
    --threads "${DNBC4TOOLS_THREADS:-8}"
done

update_workflow_status "alignment_completed" ""
