#!/usr/bin/env python3
"""Regenerate the final all-nine numerical ledger from exact integer formulas."""
from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INTERNAL = HERE / "internal"
subprocess.run(
    [sys.executable, str(INTERNAL / "snova_master_upgrade_v7_primary.py")],
    cwd=INTERNAL,
    check=True,
)
source = INTERNAL / "snova_master_upgrade_v7_primary.json"
target = HERE / "primary_ledger.json"
shutil.copyfile(source, target)
source.unlink(missing_ok=True)
print(f"wrote {target}")
