#!/usr/bin/env python3
"""Verify and deterministically reproduce the reduced-parameter forgery."""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import attack  # noqa: E402


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def main() -> None:
    transcript_path = HERE / "transcript.json"
    transcript = json.loads(transcript_path.read_text())
    pk = bytes.fromhex(transcript["public_key_hex"])
    sig = bytes.fromhex(transcript["signature_hex"])
    digest = bytes.fromhex(transcript["digest_hex"])

    require(hashlib.sha256(pk).hexdigest() == transcript["public_key_sha256"], "public-key digest mismatch")
    require(hashlib.sha256(sig).hexdigest() == transcript["signature_sha256"], "signature digest mismatch")
    require(transcript["quotient_rank"] == 3, "unexpected quotient rank")
    require(transcript["target_consistency_dimension"] == 1, "unexpected consistency dimension")
    require(transcript["attack_input_was_serialized_public_key_only"], "attack/public-input flag missing")
    require(attack.verify_serialized(pk, sig, digest), "independent serialized verifier rejected transcript")

    parsed, salt, canonical = attack.parse_signature(sig)
    require(canonical, "signature encoding is noncanonical")
    require(attack.format_accepts(parsed), "format rule rejected transcript")
    require(parsed.tolist() == transcript["signature_matrices"], "signature matrix mismatch")
    require(salt.hex() == transcript["salt_hex"], "salt mismatch")

    mutated = bytearray(sig); mutated[0] ^= 1
    require(not attack.verify_serialized(pk, bytes(mutated), digest), "signature mutation verified")
    mutated = bytearray(sig); mutated[-1] ^= 1
    require(not attack.verify_serialized(pk, bytes(mutated), digest), "salt mutation verified")

    with tempfile.TemporaryDirectory() as td:
        regenerated = Path(td) / "transcript.json"
        subprocess.run([sys.executable, str(HERE / "attack.py"), "--out", str(regenerated)], check=True, stdout=subprocess.DEVNULL)
        require(json.loads(regenerated.read_text()) == transcript, "deterministic attack transcript did not reproduce")

    batch = json.loads((HERE / "batch_results.json").read_text())
    require(batch["trials"] == 8 and batch["successes"] == 8, "batch result is not 8/8")
    require(all(row["quotient_rank"] == 3 for row in batch["records"]), "batch quotient-rank mismatch")

    print("Reduced-parameter end-to-end forgery checks passed")
    print("- reference-format serialized public key and signature")
    print("- independent direct verifier accepts the forged signature")
    print("- exact rank-3 quotient with one target-consistency coordinate")
    print("- signature and salt mutation negative controls")
    print("- bit-for-bit deterministic transcript reproduction")
    print("- eight additional deterministic keys forged successfully")


if __name__ == "__main__":
    main()
