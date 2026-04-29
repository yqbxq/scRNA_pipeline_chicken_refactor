#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "metadata validator is not implemented yet; replace this placeholder in M2" >&2
echo "expected implementation: Rscript ${SCRIPT_DIR}/r/96_validate_metadata_consistency.R" >&2
exit 64
