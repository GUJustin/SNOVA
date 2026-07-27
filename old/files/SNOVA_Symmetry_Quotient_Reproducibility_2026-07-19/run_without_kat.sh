#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python3 repro/symmetry_attack_validation.py --rank-only --out results/validation_rank_only.json
python3 repro/official_estimator.py --out results/cost_profiles.json
