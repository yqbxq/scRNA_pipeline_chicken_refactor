#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

die "入口 workflow/90_optional_cellranger_from_fastq.sh 已废弃。请改用 workflow/10_run_cellranger_from_fastq.sh。"
