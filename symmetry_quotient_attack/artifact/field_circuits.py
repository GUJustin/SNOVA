#!/usr/bin/env python3
"""Regenerate and exhaustively verify the charged F19/F19^2/F19^4 circuits."""
from __future__ import annotations
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INTERNAL = HERE / "internal"
script = INTERNAL / "snova_field_tower_circuits_v7.py"
subprocess.run([sys.executable, str(script)], cwd=INTERNAL, check=True)
# The generator retains its historical output names internally.
net_src = INTERNAL / "snova_f19_mul_balanced_v6_netlist.json"
field_src = INTERNAL / "snova_field_tower_circuits_v6.json"
net_dst = HERE / "f19_multiplier_netlist.json"
field_dst = HERE / "field_tower_circuits.json"
shutil.copyfile(net_src, net_dst)
shutil.copyfile(field_src, field_dst)
net_src.unlink(missing_ok=True)
field_src.unlink(missing_ok=True)
net = json.loads(net_dst.read_text())
field = json.loads(field_dst.read_text())
assert net["gate_count"] == 150
assert field["F19_2"]["multiplication"] == 692
assert field["F19_4"]["multiplication"] == 2628
print("field-circuit checks passed: 150 / 692 / 2628 gates")
