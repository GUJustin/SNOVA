#!/usr/bin/env bash
# Run msolve -v 2 on a .ms system under a memory + time budget, capturing the
# F4 round-table degrees. On an 8 GB machine a full degree-~14 GB is infeasible
# (= the attack cost itself), so we record the MAX degree actually reached before
# the budget forces a kill. That max is an honest lower bound, never a fabricated
# completed degree.
#
# usage: run_msolve_watch.sh <system.ms> <logfile> [mem_cap_mb] [timeout_s]
set -uo pipefail
SYS="$1"; LOG="$2"; CAP_MB="${3:-4500}"; TIMEOUT_S="${4:-1800}"
mkdir -p "$(dirname "$LOG")"
echo "[*] msolve -v 2 on $SYS  (cap=${CAP_MB}MB timeout=${TIMEOUT_S}s)" | tee "$LOG"
date "+[*] start %Y-%m-%dT%H:%M:%S" | tee -a "$LOG"

msolve -v 2 -f "$SYS" -o /dev/null >>"$LOG" 2>&1 &
PID=$!
START=$(date +%s)
REASON="completed"
while kill -0 "$PID" 2>/dev/null; do
  RSS_KB=$(ps -o rss= -p "$PID" 2>/dev/null | tr -d ' ')
  [ -z "$RSS_KB" ] && break
  RSS_MB=$(( RSS_KB / 1024 ))
  NOW=$(date +%s); ELAPSED=$(( NOW - START ))
  if [ "$RSS_MB" -gt "$CAP_MB" ]; then
    REASON="killed:memory(${RSS_MB}MB>${CAP_MB}MB)"; kill -9 "$PID" 2>/dev/null; break
  fi
  if [ "$ELAPSED" -gt "$TIMEOUT_S" ]; then
    REASON="killed:timeout(${ELAPSED}s>${TIMEOUT_S}s)"; kill -9 "$PID" 2>/dev/null; break
  fi
  sleep 2
done
wait "$PID" 2>/dev/null
END=$(date +%s)
echo "[*] stop reason: $REASON  wall=$(( END - START ))s" | tee -a "$LOG"
# F4 round-table data rows look like:  "  <deg> <sel> <pairs> <r> x <c> <density>% ..."
# The max degree reached = first column of the last data row that appears.
MAXDEG=$(grep -Eo '^[[:space:]]*[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+ x [0-9]+' "$LOG" \
         | awk '{print $1}' | sort -n | tail -1)
echo "[*] max F4 round degree reached: ${MAXDEG:-none}" | tee -a "$LOG"
