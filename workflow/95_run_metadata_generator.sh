#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "metadata generator is not implemented yet; replace this placeholder in M3" >&2
echo "expected implementation: Rscript ${SCRIPT_DIR}/r/95_build_metadata_from_questions.R" >&2
exit 64
