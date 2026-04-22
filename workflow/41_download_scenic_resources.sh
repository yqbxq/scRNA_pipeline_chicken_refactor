#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

ensure_dir "${SCENIC_DB_DIR}" "${LOG_DIR}"

download_if_missing() {
  local url="$1"
  local out="$2"
  if [[ -f "${out}" ]]; then
    echo "已存在 ${out}，跳过下载"
    return
  fi

  if command -v curl >/dev/null 2>&1; then
    curl -L --fail --output "${out}" "${url}"
    return
  fi

  if command -v wget >/dev/null 2>&1; then
    wget -O "${out}" "${url}"
    return
  fi

  die "系统中既没有 curl，也没有 wget。"
}

download_if_missing \
  "https://raw.githubusercontent.com/aertslab/SCENIC/master/inst/extdata/hs_hgnc_tfs.txt" \
  "${SCENIC_TF_LIST}"

download_if_missing \
  "https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl" \
  "${SCENIC_MOTIF_ANN}"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
  "${SCENIC_DB_500BP}"

download_if_missing \
  "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather" \
  "${SCENIC_DB_10KB}"

update_workflow_status \
  "scenic_resources_ready" \
  "bash ${PIPELINE_ROOT}/workflow/42_run_scenic.sh"

echo "SCENIC 资源检查/下载完成"
