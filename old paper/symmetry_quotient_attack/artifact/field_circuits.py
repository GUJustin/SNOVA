#!/usr/bin/env python3
"""Regenerate the scalar netlist and extension-field recurrence ledgers.

The scalar multiplier is exhaustively evaluated by the generator.  The
F19^2/F19^4 values are composition recurrences plus tower-identity checks;
complete extension-field netlists are neither emitted nor exhaustively tested.
Use ``verify_artifact.py`` for nonmutating checks of committed outputs.
"""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INTERNAL = HERE / "internal"
script = INTERNAL / "generate_field_tower_circuits.py"
subprocess.run([sys.executable, str(script)], cwd=INTERNAL, check=True)
net = json.loads((HERE / "f19_multiplier_netlist.json").read_text())
field = json.loads((HERE / "field_tower_circuits.json").read_text())
assert net["gate_count"] == 150
assert field["F19_2"]["multiplication"] == 692
assert field["F19_4"]["multiplication"] == 2628
print("field-circuit checks passed: 150 / 692 / 2628 gates")
