#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
KAT="${1:-source_snapshots/PQCsignKAT_SNOVA_28_5_19_4.txt}"
if [[ ! -f "$KAT" ]]; then
  echo "Missing official KAT response file: $KAT" >&2
  echo "See source_snapshots/README.md for the pinned repository commit and expected filename." >&2
  exit 2
fi
python3 repro/symmetry_attack_validation.py --kat "$KAT" --out results/validation_results.json
python3 repro/cross_column_certificate.py --kat "$KAT" --out results/cross_column_certificate.json
