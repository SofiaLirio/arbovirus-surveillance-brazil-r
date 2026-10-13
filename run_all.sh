#!/usr/bin/env bash
# Full rebuild from raw data: packages -> raw->db ETL -> all analyses (R).
# Usage (from any folder):  bash run_all.sh [--skip-etl] [--skip-install]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RSCRIPT="$(command -v Rscript || true)"
if [ -z "$RSCRIPT" ]; then
  RSCRIPT="$(ls -d /c/Program\ Files/R/R-*/bin/Rscript.exe 2>/dev/null | sort | tail -1 || true)"
fi
[ -n "$RSCRIPT" ] || { echo "Rscript not found: install R 4.5+ from https://cloud.r-project.org" >&2; exit 1; }
"$RSCRIPT" "$ROOT/run_all.R" "$@"
