#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "$0")" && pwd)/lib/common.sh"

die "入口 workflow/97_build_ortholog_cache.sh 已废弃。请改用 workflow/40_build_ortholog_cache.sh。"
