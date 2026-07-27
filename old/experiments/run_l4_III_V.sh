#!/usr/bin/env bash
# Solving-degree experiment for the III_l4 (70/55) and V_l4 (90/73) square-spec
# cores, each vs 5 matched random controls. Runs AFTER the in-flight I_l4
# controls clear (8 GB RAM hosts only one msolve at a time). Sequential.
set -uo pipefail
cd /Users/jthaler/Dropbox/NIST-break-papers/SNOVA/experiments

# wait for any running msolve (the I_l4 control sweep) to finish
while pgrep -x msolve >/dev/null 2>&1; do sleep 30; done
echo "[driver] prior msolve cleared at $(date)"

MEM=5500
TMO=1800

run() { ./run_msolve_watch.sh "$1" "$2" "$MEM" "$TMO"; }

echo "===== III_l4 official core 70/55 ====="
run systems/III_l4_core_70in55.ms logs/official_III_l4_70in55.log
for k in 00 01 02 03 04; do
  echo "===== III_l4 control $k (70/55) ====="
  run "systems/controls/ctrl_70in55_${k}.ms" "logs/ctrl_70in55_${k}.log"
done

echo "===== V_l4 official core 90/73 ====="
run systems/V_l4_core_90in73.ms logs/official_V_l4_90in73.log
for k in 00 01 02 03 04; do
  echo "===== V_l4 control $k (90/73) ====="
  run "systems/controls/ctrl_90in73_${k}.ms" "logs/ctrl_90in73_${k}.log"
done

echo "===== III_l4 + V_l4 sweep done at $(date) ====="
