#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

die "入口 workflow/98_download_scenic_resources.sh 已废弃。请改用 workflow/41_download_scenic_resources.sh。"
