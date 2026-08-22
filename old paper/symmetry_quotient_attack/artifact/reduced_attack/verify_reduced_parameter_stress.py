#!/usr/bin/env python3
"""Fresh-key stress checks for the reduced zero-offset forgery.

This checker deliberately has no committed expected-output ledger.  It creates
fresh deterministic reduced keys, runs the public attack, reconstructs every
verification equation through the KAT-anchored literal evaluator, and fails on
the first discrepancy.  It remains a toy-size Python-transcription test.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np

import reduced_parameter_end_to_end_forgery as toy
import snova_v23_reference as ref


FORGERY_TRIALS = 24
RANK_TRIALS = 100


def main() -> None:
    forgery_records = []
    abq = ref.reconstruct_fixed_abq(toy.PARAMS)
    for trial in range(FORGERY_TRIALS):
        message = f"independent corrected message {trial}".encode()
        public_key, _metadata = toy.reduced_keygen(
            f"independent corrected key {trial}".encode()
        )
        signature, solve = toy.forge(
            public_key,
            message,
            np.random.default_rng(0xA11CE000 + trial),
        )
        digits = ref.decode_base19(
            signature[:-16], toy.PARAMS.variables * toy.PARAMS.r
        )
        candidate = digits.reshape(toy.PARAMS.variables, toy.PARAMS.r)
        output = ref.literal_public_output(
            toy.PARAMS,
            abq,
            ref.expand_public_key(public_key, toy.PARAMS),
            candidate,
        )
        target = ref.message_target(
            public_key[:16], message, signature[-16:], toy.PARAMS.outputs
        )
        if not np.array_equal(output, target):
            raise AssertionError(f"fresh forgery {trial} failed literal evaluation")
        if ref.rejection_count(candidate, toy.PARAMS) > toy.PARAMS.n // 4:
            raise AssertionError(f"fresh forgery {trial} failed the rejection rule")
        changed_key = bytearray(public_key)
        changed_key[0] ^= 1
        if toy.literal_verify(bytes(changed_key), message, signature):
            raise AssertionError(f"fresh forgery {trial} accepted under a changed key")
        forgery_records.append(
            {
                "public_key_sha256": hashlib.sha256(public_key).hexdigest(),
                "salt_counter": solve["salt_counter"],
                "signature_sha256": hashlib.sha256(signature).hexdigest(),
                "slices_tested": solve["slices_tested"],
                "trial": trial,
            }
        )

    ranks = []
    for trial in range(RANK_TRIALS):
        public_key, _metadata = toy.reduced_keygen(
            f"rank census key {trial}".encode()
        )
        coefficients, _pairs = toy.interpolate_restricted_map(
            abq, ref.expand_public_key(public_key, toy.PARAMS)
        )
        ranks.append(ref.rank_mod(coefficients))
    if any(rank != toy.PARAMS.unordered for rank in ranks):
        raise AssertionError(f"unexpected restricted ranks: {sorted(set(ranks))}")

    rng = np.random.default_rng(20260731)
    for length in range(1, 101):
        values = rng.integers(0, ref.Q, size=length, dtype=np.int64)
        packed = ref.encode_base19(values)
        if not np.array_equal(ref.decode_base19(packed, length), values):
            raise AssertionError(f"base-19 round trip failed at length {length}")

    print(
        json.dumps(
            {
                "base19_round_trip_lengths": [1, 100],
                "forgery_successes": len(forgery_records),
                "forgery_trials": FORGERY_TRIALS,
                "maximum_salt_counter": max(
                    record["salt_counter"] for record in forgery_records
                ),
                "maximum_slices_tested": max(
                    record["slices_tested"] for record in forgery_records
                ),
                "rank_3_successes": sum(rank == 3 for rank in ranks),
                "rank_trials": RANK_TRIALS,
                "scope": (
                    "fresh-key stress test at unofficial (2,1,19,2,2); "
                    "not production-size evidence"
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
