#!/usr/bin/env python3
"""Public-only Version-2.3 KAT and common-column reduction correspondence test.

The supplied KAT signature is verified as a known answer.  Separately, the
attack reduction is instantiated on its public key and checked, but the
residual quadratic system is not solved and no forgery is claimed.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

import snova_v23_reference as ref


ARTIFACT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_KAT = ARTIFACT_ROOT / "official/PQCsignKAT_SNOVA_28_5_19_4.txt"
DEFAULT_EXPECTED = Path(__file__).with_name(
    "official_key_reduction_harness_output.json"
)


def _archived_reset_improve(
    matrix: np.ndarray, S: np.ndarray, ell: int
) -> np.ndarray:
    """Reproduce the archived Python bug: retry from the original matrix."""
    original = np.asarray(matrix, dtype=np.int64) % ref.Q
    if original.shape != (ell, ell) or ref.determinant_mod(original):
        return original.copy()
    for coefficient in range(1, ref.Q):
        candidate = (original + coefficient * S) % ref.Q
        if ref.determinant_mod(candidate):
            return candidate
    raise RuntimeError("archived reset-style matrix improvement failed")


def _archived_reset_abq(params: ref.Params) -> ref.ABQ:
    """Rebuild fixed ABQ with the archived reset-on-each-retry semantics."""
    correct = ref.reconstruct_fixed_abq(params)
    count = params.o * params.alpha * (
        params.r * params.r + params.r * params.l + 2 * params.l
    )
    raw = np.frombuffer(
        hashlib.shake_256(b"SNOVA_ABQ").digest(count), dtype=np.uint8
    ).astype(np.int64) % ref.Q
    position = 0
    size = params.o * params.alpha * params.r * params.r
    raw_A = raw[position:position + size].reshape(
        params.o, params.alpha, params.r, params.r
    )
    position += size
    size = params.o * params.alpha * params.r * params.l
    raw_B = raw[position:position + size].reshape(
        params.o, params.alpha, params.r, params.l
    )
    position += size + 2 * params.o * params.alpha * params.l
    if position != count:
        raise AssertionError("archived-reset ABQ expansion length mismatch")
    A = np.empty_like(raw_A)
    B = np.empty_like(raw_B)
    for output in range(params.o):
        for alpha in range(params.alpha):
            A[output, alpha] = _archived_reset_improve(
                raw_A[output, alpha], correct.S, params.l
            )
            B[output, alpha] = _archived_reset_improve(
                raw_B[output, alpha], correct.S, params.l
            )
    return ref.ABQ(A, B, correct.q1, correct.q2, correct.S)


def archived_reset_regression(
    kat: dict[str, str],
    params: ref.Params,
    public_matrices: np.ndarray,
    correct_abq: ref.ABQ,
) -> dict[str, object]:
    """Show that the archived reset semantics are not verifier-compatible."""
    public_key = bytes.fromhex(kat["pk"])
    message = bytes.fromhex(kat["msg"])
    signed_message = bytes.fromhex(kat["sm"])
    serialized = signed_message[:-len(message)]
    signature = ref.decode_base19(
        serialized[:-16], params.variables * params.r
    ).reshape(params.variables, params.r)
    salt = serialized[-16:]
    target = ref.message_target(public_key[:16], message, salt, params.outputs)
    reset_abq = _archived_reset_abq(params)
    reset_output = ref.literal_public_output(
        params, reset_abq, public_matrices, signature
    )
    mismatches = np.flatnonzero(reset_output != target).tolist()
    matching = params.outputs - len(mismatches)
    rho = np.asarray([1, 8, 9, 14], dtype=np.int64)
    correct_quotient, _ = ref.build_quotient(params, correct_abq, rho)
    reset_quotient, _ = ref.build_quotient(params, reset_abq, rho)
    quotient_differences = int(np.count_nonzero(correct_quotient != reset_quotient))
    abq_differences = int(
        np.count_nonzero(correct_abq.A != reset_abq.A)
        + np.count_nonzero(correct_abq.B != reset_abq.B)
    )
    if matching != 65 or quotient_differences != 144 or abq_differences != 12:
        raise AssertionError("archived-reset regression no longer matches the pinned KAT")
    return {
        "status": "archived reset-on-each-retry helper is verifier-incompatible",
        "kat_coordinates_matching": matching,
        "kat_coordinates_total": params.outputs,
        "kat_mismatch_coordinates": mismatches,
        "fixed_abq_scalar_differences": abq_differences,
        "quotient_scalar_differences": quotient_differences,
        "wrong_semantics_rejected": True,
    }


def known_answer_verification(
    kat: dict[str, str], params: ref.Params, public_matrices, abq
) -> dict[str, object]:
    public_key = bytes.fromhex(kat["pk"])
    message = bytes.fromhex(kat["msg"])
    signed_message = bytes.fromhex(kat["sm"])
    if signed_message[-len(message):] != message:
        raise AssertionError("KAT signed-message suffix does not equal msg")
    signature_bytes = signed_message[:-len(message)]
    expected = 16 + ref.encoded_length(params.variables * params.r)
    if len(signature_bytes) != expected:
        raise AssertionError("unexpected KAT signature length")
    packed_signature = signature_bytes[:-16]
    salt = signature_bytes[-16:]
    digits = ref.decode_base19(packed_signature, params.variables * params.r)
    if ref.encode_base19(digits) != packed_signature:
        raise AssertionError("signature serialization round trip failed")
    signature = digits.reshape(params.variables, params.r)
    output = ref.literal_public_output(params, abq, public_matrices, signature)
    target = ref.message_target(public_key[:16], message, salt, params.outputs)
    symmetric_blocks = ref.rejection_count(signature, params)
    threshold = params.n // 4 if params.l == 2 else 0
    accepted = symmetric_blocks <= threshold and np.array_equal(output, target)
    altered_target = ref.message_target(
        public_key[:16], message + b"\x00", salt, params.outputs
    )
    altered_rejected = not np.array_equal(output, altered_target)
    if not accepted or not altered_rejected:
        raise AssertionError("known-answer positive/negative test failed")
    return {
        "status": "known-answer verification only; supplied signer output",
        "accepted": True,
        "hash_matches_public_map": True,
        "signature_round_trip": True,
        "symmetric_blocks": int(symmetric_blocks),
        "rejection_threshold": int(threshold),
        "altered_message_rejected": True,
        "signature_bytes": len(signature_bytes),
        "message_bytes": len(message),
    }


def reduction_correspondence(
    public_seed: bytes,
    params: ref.Params,
    public_matrices: np.ndarray,
    abq: ref.ABQ,
    seed: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    rho = np.asarray([1, 8, 9, 14], dtype=np.int64)
    quotient, labels = ref.build_quotient(params, abq, rho)
    quotient_rank = ref.rank_mod(quotient)
    left_kernel = ref.nullspace_mod(quotient.T)
    offset, format_check = ref.safe_offsets(params, rho, rng)
    constant = ref.attack_public_output(params, abq, public_matrices, offset)
    linear = np.zeros((params.outputs, params.variables), dtype=np.int64)
    for index in range(params.variables):
        basis = np.zeros(params.variables, dtype=np.int64)
        basis[index] = 1
        related = ref.common_column_signature(basis, rho, offset)
        linear[:, index] = (
            ref.attack_public_output(params, abq, public_matrices, related)
            - constant
            - quotient
            @ ref.quotient_coordinates(params, public_matrices, basis, labels)
        ) % ref.Q

    for _ in range(32):
        x = rng.integers(0, ref.Q, size=params.variables, dtype=np.int64)
        related = ref.common_column_signature(x, rho, offset)
        attack_output = ref.attack_public_output(
            params, abq, public_matrices, related
        )
        literal_output = ref.literal_public_output(
            params, abq, public_matrices, related
        )
        decomposed = (
            constant
            + linear @ x
            + quotient
            @ ref.quotient_coordinates(params, public_matrices, x, labels)
        ) % ref.Q
        if not np.array_equal(attack_output, literal_output):
            raise AssertionError("attack evaluator and literal verifier disagree")
        if not np.array_equal(attack_output, decomposed):
            raise AssertionError("substituted-verifier decomposition failed")

    constraints = left_kernel @ linear % ref.Q
    constraint_rank = ref.rank_mod(constraints)
    target = ref.message_target(
        public_seed,
        b"SNOVA public-only reduction correspondence test",
        bytes.fromhex("00112233445566778899aabbccddeeff"),
        params.outputs,
    )
    base, kernel, _free, pivots = ref.affine_solve(
        constraints, left_kernel @ (target - constant) % ref.Q
    )
    if not np.array_equal(
        constraints @ base % ref.Q,
        left_kernel @ (target - constant) % ref.Q,
    ):
        raise AssertionError("affine target solve failed")
    return {
        "status": "official-key reduction correspondence; residual MQ not solved; no forgery",
        "public_inputs_only": True,
        "relation_rho": rho.tolist(),
        "quotient_shape": list(quotient.shape),
        "quotient_rank": int(quotient_rank),
        "left_kernel_dimension": int(left_kernel.shape[0]),
        "affine_constraint_rank": int(constraint_rank),
        "affine_kernel_dimension": int(kernel.shape[1]),
        "residual_quadratics": int(quotient_rank),
        "residual_variables": int(kernel.shape[1]),
        "substituted_verifier_random_checks": 32,
        "literal_verifier_cross_checks": 32,
        "format_offset_check": format_check,
        "affine_pivot_count": len(pivots),
        "residual_system_solved": False,
        "forgery_output": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--kat", type=Path, default=DEFAULT_KAT)
    parser.add_argument("--record", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260731)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--check",
        action="store_true",
        help="compare with committed JSON without writing files",
    )
    modes.add_argument(
        "--write",
        action="store_true",
        help="regenerate the committed JSON",
    )
    parser.add_argument("--expected", type=Path, default=DEFAULT_EXPECTED)
    args = parser.parse_args()

    params = ref.PUBLISHED_PARAMS[0]
    kat = ref.parse_public_kat_record(args.kat, args.record)
    public_key = bytes.fromhex(kat["pk"])
    public_matrices = ref.expand_public_key(public_key, params)
    abq = ref.reconstruct_fixed_abq(params)

    report = {
        "scope": (
            "public-key/verifier/reduction correspondence only; KAT signature is not a forgery; "
            "attack residual system is not solved"
        ),
        "parameter": asdict(params),
        "kat_file": args.kat.name,
        "kat_record": args.record,
        "kat_public_key_sha256": hashlib.sha256(public_key).hexdigest(),
        "public_key_bytes": len(public_key),
        "public_expansion_uses_secret": False,
        "abq_uses_pinned_reference_cumulative_improve_semantics": True,
        "archived_reset_semantics_regression": archived_reset_regression(
            kat, params, public_matrices, abq
        ),
        "known_answer_verification": known_answer_verification(
            kat, params, public_matrices, abq
        ),
        "attack_reduction_correspondence": reduction_correspondence(
            public_key[:16], params, public_matrices, abq, args.seed
        ),
    }
    output = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        expected = args.expected.read_text()
        if output != expected:
            difference = "".join(
                difflib.unified_diff(
                    expected.splitlines(keepends=True),
                    output.splitlines(keepends=True),
                    fromfile=args.expected.name,
                    tofile="computed",
                )
            )
            raise SystemExit(
                "correspondence output differs from committed JSON:\n" + difference
            )
        print(f"OK: {args.expected.name}")
    elif args.write:
        args.expected.parent.mkdir(parents=True, exist_ok=True)
        args.expected.write_text(output)
        print(f"WROTE: {args.expected.name}")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
