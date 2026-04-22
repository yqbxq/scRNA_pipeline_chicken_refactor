#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

require_status_flag_or_warn \
  "status.velocity_upstream_ready" \
  "RNA velocity 上游尚未 ready。请先确保 BAM / Cell Ranger 输出已通过 intake 审计。"

shopt -s nullglob

ensure_dir "${VELOCITY_LOOM_DIR}" "${LOG_DIR}"

run_in_conda_prefix "${VELOCITY_ENV_PREFIX}" velocyto --help >/dev/null

IFS=',' read -r -a raw_samples <<< "${RAW_SAMPLES}"

for sample in "${raw_samples[@]}"; do
  sample_dir="${CELLRANGER_OUT_DIR}/${sample}"
  source_loom="${sample_dir}/velocyto/${sample}.loom"
  target_loom="${VELOCITY_LOOM_DIR}/${sample}.loom"
  cellsorted_bam="${sample_dir}/outs/cellsorted_possorted_genome_bam.bam"
  cellsorted_tmp=("${sample_dir}/outs"/cellsorted_possorted_genome_bam.bam.tmp.*)

  [[ -d "${sample_dir}" ]] || die "Missing Cell Ranger output dir: ${sample_dir}"
  [[ -f "${sample_dir}/outs/possorted_genome_bam.bam" ]] || die "Missing BAM file: ${sample_dir}/outs/possorted_genome_bam.bam"

  if [[ -f "${target_loom}" ]]; then
    echo "Existing target loom found, skip ${sample}: ${target_loom}"
    continue
  fi

  if [[ -f "${source_loom}" ]]; then
    echo "Existing source loom found, link ${sample}: ${source_loom}"
    ln -sf "${source_loom}" "${target_loom}"
    continue
  fi

  [[ -f "${CLEAN_GTF}" ]] || die "Missing cleaned GTF: ${CLEAN_GTF}"

  related_jobs="$(pgrep -af "velocyto run10x .*${sample_dir}|samtools sort .*${sample_dir}/outs/possorted_genome_bam.bam" || true)"
  if [[ -n "${related_jobs}" ]]; then
    die "Detected running velocyto/samtools jobs for ${sample}: ${related_jobs}"
  fi

  if [[ -f "${cellsorted_bam}" || -f "${cellsorted_bam}.bai" || ${#cellsorted_tmp[@]} -gt 0 ]]; then
    echo "Clean stale cellsorted BAM temp files for ${sample}"
    rm -f "${cellsorted_bam}" "${cellsorted_bam}.bai" "${cellsorted_tmp[@]}"
  fi

  echo "Start velocyto for ${sample}"
  run_in_conda_prefix "${VELOCITY_ENV_PREFIX}" \
    velocyto run10x -@ "${VELOCYTO_THREADS}" "${sample_dir}" "${CLEAN_GTF}"

  [[ -f "${source_loom}" ]] || die "velocyto finished but loom not found: ${source_loom}"
  ln -sf "${source_loom}" "${target_loom}"
  echo "Finished ${sample}: ${target_loom}"
done

update_workflow_status \
  "velocyto_completed" \
  "bash ${PIPELINE_ROOT}/workflow/31_run_scvelo.sh"

echo "velocyto run completed"
