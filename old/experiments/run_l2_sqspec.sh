#!/usr/bin/env bash
# SUPPLEMENTARY structural-check sweep: square-specialization cores of the l=2
# official keys (I 48/36, III 72/60, V 96/79) vs 5 matched random controls each.
# These are NOT the paper's l=2 gate-cost cores (those are Hashimoto blocks); they
# are a valid alternative attack path, run only to check for exploitable structure
# in l=2 key material. Chained to start after the l=4 sweep (run_l4_III_V.sh) exits.
set -uo pipefail
cd /Users/jthaler/Dropbox/NIST-break-papers/SNOVA/experiments

# wait until the l=4 driver process is gone (it, in turn, waits out the I_l4 controls)
while pgrep -f run_l4_III_V.sh >/dev/null 2>&1; do sleep 60; done
# and until no msolve is mid-run
while pgrep -x msolve >/dev/null 2>&1; do sleep 30; done
echo "[l2 driver] l=4 sweep cleared at $(date)"

MEM=5500; TMO=1800
run() { ./run_msolve_watch.sh "$1" "$2" "$MEM" "$TMO"; }

emit() {  # <shape_tag> <m> <n>
  local tag=$1 m=$2 n=$3
  echo "===== ${tag} square-spec core ${m}/${n} ====="
  run "systems/${tag}_sqspec_core_${m}in${n}.ms" "logs/sqspec_${tag}_${m}in${n}.log"
  for k in 00 01 02 03 04; do
    echo "===== ${tag} control $k (${m}/${n}) ====="
    run "systems/controls/ctrl_${m}in${n}_${k}.ms" "logs/ctrl_${m}in${n}_${k}.log"
  done
}

emit I_l2   48 36
emit III_l2 72 60
emit V_l2   96 79
echo "===== l=2 square-spec sweep done at $(date) ====="
