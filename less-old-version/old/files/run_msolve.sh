#!/usr/bin/env bash
# msolve alternative for concern #1 (run where msolve is installed).
# INPUT: system.ms in msolve format:
#   line 1: comma-separated variable names
#   line 2: field characteristic (19)
#   remaining: one polynomial per line, terminated by commas, last has no trailing comma
# msolve -v 2 prints the degrees reached by its F4; grep them out.
set -euo pipefail
SYS="${1:-system.ms}"
echo "[*] running msolve -v 2 on $SYS"
msolve -v 2 -f "$SYS" -o /dev/null 2> msolve.log || true
echo "[*] degrees reached (compare to semi-regular D_reg):"
grep -Ei "degree|D_reg|step" msolve.log || echo "  (inspect msolve.log manually)"
