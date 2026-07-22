#!/usr/bin/env bash
# Run the 5 matched random controls at the SAME budget as the official core,
# sequentially (8 GB RAM cannot host two msolve runs at once).
set -uo pipefail
cd /Users/jthaler/Dropbox/NIST-break-papers/SNOVA/experiments
for k in 01 02 03 04; do
  echo "===== control $k ====="
  ./run_msolve_watch.sh "systems/controls/ctrl_50in41_${k}.ms" \
      "logs/ctrl_50in41_${k}.log" 5500 1500
done
echo "===== all controls done ====="
