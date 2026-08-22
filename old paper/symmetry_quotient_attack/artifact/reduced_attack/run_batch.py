#!/usr/bin/env python3
"""Run the reduced attack on eight additional deterministic toy keys."""
import json
import time
from pathlib import Path

import numpy as np
import attack as demo

HERE = Path(__file__).resolve().parent
OUT = HERE / "batch_results.json"
records = []
start = time.time()
for trial in range(8):
    seed = demo.shake(f"SNOVA reduced batch key {trial}".encode(), 48)
    pk, _ = demo.keygen_bytes(seed)
    digest = f"reduced SNOVA batch message {trial}".encode()
    rng = np.random.default_rng(0x534E4F560000 + trial)
    t0 = time.time()
    try:
        transcript = demo.attack(pk, digest, rng, max_salts=10000, max_slices=2000)
        records.append({
            "trial": trial,
            "success": True,
            "seconds": time.time() - t0,
            "pk_sha256": transcript["public_key_sha256"],
            "sig_sha256": transcript["signature_sha256"],
            "salt_counter": transcript["salt_counter"],
            "slices_tested": transcript["slices_tested"],
            "quotient_rank": transcript["quotient_rank"],
        })
    except Exception as exc:
        records.append({"trial": trial, "success": False, "seconds": time.time() - t0, "error": repr(exc)})
result = {
    "trials": len(records),
    "successes": sum(row["success"] for row in records),
    "total_seconds": time.time() - start,
    "records": records,
}
OUT.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
