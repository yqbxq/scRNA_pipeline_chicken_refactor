#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

die "入口 workflow/00_prepare_author_input.sh 已废弃。请改用 workflow/04_standardize_inputs.sh。"
