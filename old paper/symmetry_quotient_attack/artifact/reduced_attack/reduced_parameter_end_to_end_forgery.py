#!/usr/bin/env python3
"""Zero-offset toy forgery for an unofficial reduced SNOVA shape.

The shape is ``(v,o,q,ell,r)=(2,1,19,2,2)``.  All scheme-facing operations
(public expansion, hidden-UOV key generation, fixed ABQ data, target hashing,
packing, rejection, and literal verification) use the KAT-anchored Version 2.3
helpers in :mod:`snova_v23_reference`.  Only the public attack logic is local:
it interpolates the complete zero-offset restriction ``X_i=[u_i|0]``, extracts
its rank-three output image, filters salts by the remaining consistency
coordinate, and exhaustively searches public three-dimensional slices.

This is composition evidence at a deliberately tiny, unofficial shape.  It is
not a production-size solver, a cost measurement, or an invocation of the
upstream C verifier.
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

import snova_v23_reference as ref


PARAMS = ref.Params("reduced-v2.3-zero-offset-instance", 2, 1, 2, 2)
RHO = np.asarray([1, 0], dtype=np.int64)
DEFAULT_EXPECTED = Path(__file__).with_name(
    "reduced_parameter_end_to_end_forgery_output.json"
)
DEFAULT_KEY_LABEL = b"SNOVA reduced attack key seed v2"
DEFAULT_MESSAGE = b"chosen message for reduced SNOVA forgery"
DEFAULT_RNG_SEED = 0x534E4F5641
BATCH_TRIALS = 8


def reduced_keygen(label: bytes) -> tuple[bytes, dict]:
    """Create a valid reduced public key outside the attack boundary."""
    seed = hashlib.shake_256(label).digest(48)
    public_key, generated = ref.derive_public_key_from_seeds(
        seed[:16], seed[16:], PARAMS
    )
    parsed = ref.expand_public_key(public_key, PARAMS)
    if not np.array_equal(generated, parsed):
        raise AssertionError("compressed public-key round trip failed")
    return public_key, {
        "key_seed_sha256": hashlib.sha256(seed).hexdigest(),
        "pk_seed_hex": seed[:16].hex(),
        "public_key_bytes": len(public_key),
    }


def signature_from_u(u: np.ndarray) -> np.ndarray:
    """Return the zero-offset common-column signature ``X_i=[u_i|0]``."""
    signature = np.zeros((PARAMS.variables, PARAMS.r), dtype=np.int64)
    signature[:, 0] = np.asarray(u, dtype=np.int64).reshape(-1) % ref.Q
    return signature


def monomial_pairs(number_of_variables: int) -> list[tuple[int, int]]:
    return [
        (left, right)
        for left in range(number_of_variables)
        for right in range(left, number_of_variables)
    ]


def monomial_vector(point: np.ndarray, pairs) -> np.ndarray:
    point = np.asarray(point, dtype=np.int64).reshape(-1) % ref.Q
    return np.asarray(
        [int(point[left]) * int(point[right]) % ref.Q for left, right in pairs],
        dtype=np.int64,
    )


def monomial_matrix(points: np.ndarray, pairs) -> np.ndarray:
    points = np.asarray(points, dtype=np.int64) % ref.Q
    return np.stack(
        [(points[:, left] * points[:, right]) % ref.Q for left, right in pairs],
        axis=1,
    )


def row_echelon_transform(matrix: np.ndarray):
    """Reduced row echelon form together with the applied output transform."""
    reduced = np.asarray(matrix, dtype=np.int64).copy() % ref.Q
    rows, columns = reduced.shape
    transform = np.eye(rows, dtype=np.int64)
    row = 0
    pivots: list[int] = []
    for column in range(columns):
        pivot = next(
            (
                candidate
                for candidate in range(row, rows)
                if reduced[candidate, column] % ref.Q
            ),
            None,
        )
        if pivot is None:
            continue
        if pivot != row:
            reduced[[row, pivot]] = reduced[[pivot, row]]
            transform[[row, pivot]] = transform[[pivot, row]]
        scale = ref.inv_mod(int(reduced[row, column]))
        reduced[row] = reduced[row] * scale % ref.Q
        transform[row] = transform[row] * scale % ref.Q
        for other in range(rows):
            if other != row and reduced[other, column] % ref.Q:
                factor = int(reduced[other, column])
                reduced[other] = (reduced[other] - factor * reduced[row]) % ref.Q
                transform[other] = (
                    transform[other] - factor * transform[row]
                ) % ref.Q
        pivots.append(column)
        row += 1
        if row == rows:
            break
    return reduced, transform % ref.Q, row, pivots


def interpolate_restricted_map(
    abq: ref.ABQ, public_matrices: np.ndarray
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Interpolate every coefficient of the homogeneous restricted verifier."""
    pairs = monomial_pairs(PARAMS.variables)
    coefficients = np.zeros((PARAMS.outputs, len(pairs)), dtype=np.int64)
    basis_values = []
    for variable in range(PARAMS.variables):
        basis = np.zeros(PARAMS.variables, dtype=np.int64)
        basis[variable] = 1
        basis_values.append(
            ref.attack_public_output(
                PARAMS, abq, public_matrices, signature_from_u(basis)
            )
        )
    pair_to_column = {pair: column for column, pair in enumerate(pairs)}
    for variable in range(PARAMS.variables):
        coefficients[:, pair_to_column[(variable, variable)]] = basis_values[
            variable
        ]
    for left in range(PARAMS.variables):
        for right in range(left + 1, PARAMS.variables):
            point = np.zeros(PARAMS.variables, dtype=np.int64)
            point[left] = point[right] = 1
            value = ref.attack_public_output(
                PARAMS, abq, public_matrices, signature_from_u(point)
            )
            coefficients[:, pair_to_column[(left, right)]] = (
                value - basis_values[left] - basis_values[right]
            ) % ref.Q
    return coefficients, pairs


def literal_verify(public_key: bytes, message: bytes, signature: bytes) -> bool:
    """Serialized verifier through the distinct literal-evaluator path."""
    expected = ref.encoded_length(PARAMS.variables * PARAMS.r) + 16
    if len(signature) != expected:
        return False
    packed, salt = signature[:-16], signature[-16:]
    try:
        digits = ref.decode_base19(packed, PARAMS.variables * PARAMS.r)
        public_matrices = ref.expand_public_key(public_key, PARAMS)
    except ValueError:
        return False
    candidate = digits.reshape(PARAMS.variables, PARAMS.r)
    if ref.rejection_count(candidate, PARAMS) > PARAMS.n // 4:
        return False
    output = ref.literal_public_output(
        PARAMS,
        ref.reconstruct_fixed_abq(PARAMS),
        public_matrices,
        candidate,
    )
    target = ref.message_target(
        public_key[:16], message, salt, PARAMS.outputs
    )
    return bool(np.array_equal(output, target))


def forge(
    public_key: bytes,
    message: bytes,
    rng: np.random.Generator,
    *,
    maximum_salts: int = 10_000,
    maximum_slices: int = 2_000,
) -> tuple[bytes, dict]:
    """Forge using only a serialized public key, message, and public randomness."""
    public_matrices = ref.expand_public_key(public_key, PARAMS)
    abq = ref.reconstruct_fixed_abq(PARAMS)

    # The staged attack evaluator and literal verifier evaluator must agree on
    # unrestricted inputs before the restricted system is trusted.
    for _ in range(32):
        candidate = rng.integers(
            0,
            ref.Q,
            size=(PARAMS.variables, PARAMS.r),
            dtype=np.int64,
        )
        staged = ref.attack_public_output(PARAMS, abq, public_matrices, candidate)
        literal = ref.literal_public_output(PARAMS, abq, public_matrices, candidate)
        if not np.array_equal(staged, literal):
            raise AssertionError("staged and literal public evaluators disagree")

    coefficients, pairs = interpolate_restricted_map(abq, public_matrices)
    reduced, output_transform, coefficient_rank, _pivots = row_echelon_transform(
        coefficients
    )
    quotient, quotient_labels = ref.build_quotient(PARAMS, abq, RHO)
    quotient_rank = ref.rank_mod(quotient)
    if coefficient_rank != PARAMS.unordered or quotient_rank != PARAMS.unordered:
        raise RuntimeError(
            f"rank preflight failed: coefficients={coefficient_rank}, "
            f"quotient={quotient_rank}, expected={PARAMS.unordered}"
        )
    if ref.rank_mod(np.column_stack([coefficients, quotient])) != quotient_rank:
        raise AssertionError("interpolated and explicit quotient images differ")
    if np.any(reduced[coefficient_rank:] % ref.Q):
        raise AssertionError("row elimination failed")

    for _ in range(32):
        point = rng.integers(0, ref.Q, size=PARAMS.variables, dtype=np.int64)
        interpolated = coefficients @ monomial_vector(point, pairs) % ref.Q
        candidate = signature_from_u(point)
        staged = ref.attack_public_output(PARAMS, abq, public_matrices, candidate)
        literal = ref.literal_public_output(PARAMS, abq, public_matrices, candidate)
        if not np.array_equal(interpolated, staged):
            raise AssertionError("restricted interpolation mismatch")
        if not np.array_equal(interpolated, literal):
            raise AssertionError("restricted interpolation/literal mismatch")

    salt = None
    target = None
    transformed_target = None
    salt_counter = None
    for counter in range(maximum_salts):
        candidate_salt = counter.to_bytes(16, "little")
        candidate_target = ref.message_target(
            public_key[:16], message, candidate_salt, PARAMS.outputs
        )
        candidate_transformed = output_transform @ candidate_target % ref.Q
        if np.all(candidate_transformed[coefficient_rank:] == 0):
            salt = candidate_salt
            target = candidate_target
            transformed_target = candidate_transformed
            salt_counter = counter
            break
    if salt is None:
        raise RuntimeError("no target satisfying the output-image constraint")

    grid = np.asarray(
        list(itertools.product(range(ref.Q), repeat=coefficient_rank)),
        dtype=np.int64,
    )
    solution = None
    chosen_slice = None
    chosen_point = None
    slices_tested = 0
    for _ in range(maximum_slices):
        slice_matrix = rng.integers(
            0,
            ref.Q,
            size=(PARAMS.variables, coefficient_rank),
            dtype=np.int64,
        )
        if ref.rank_mod(slice_matrix) != coefficient_rank:
            continue
        slices_tested += 1
        points = grid @ slice_matrix.T % ref.Q
        features = monomial_matrix(points, pairs)
        values = features @ coefficients.T % ref.Q
        matches = np.all(values == target[None, :], axis=1)
        # X_i=[[u_2i,0],[u_2i+1,0]] is symmetric iff u_2i+1=0.
        # Here floor(n/4)=0, so every block must be nonsymmetric.
        matches &= np.all(points[:, 1::2] != 0, axis=1)
        indices = np.flatnonzero(matches)
        if indices.size:
            index = int(indices[0])
            solution = points[index]
            chosen_slice = slice_matrix
            chosen_point = grid[index]
            break
    if solution is None:
        raise RuntimeError("no root found on the searched public slices")

    candidate = signature_from_u(solution)
    signature = ref.encode_base19(candidate.reshape(-1)) + salt
    if not literal_verify(public_key, message, signature):
        raise AssertionError("literal verifier rejected the reconstructed signature")
    output = ref.literal_public_output(PARAMS, abq, public_matrices, candidate)
    if not np.array_equal(output, target):
        raise AssertionError("final public-map output differs from the target")

    return signature, {
        "common_column_vector_u": solution.tolist(),
        "explicit_quotient_labels": [list(label) for label in quotient_labels],
        "explicit_quotient_rank": quotient_rank,
        "monomial_pairs": [list(pair) for pair in pairs],
        "output_row_transform": output_transform.tolist(),
        "quotient_image_matches_interpolation": True,
        "reduced_quadratic_coefficients": coefficients.tolist(),
        "restricted_coefficient_rank": coefficient_rank,
        "salt_counter": salt_counter,
        "salt_hex": salt.hex(),
        "signature_matrices": candidate.reshape(
            PARAMS.n, PARAMS.l, PARAMS.r
        ).tolist(),
        "slice_matrix": chosen_slice.tolist(),
        "slice_point": chosen_point.tolist(),
        "slices_tested": slices_tested,
        "target": target.tolist(),
        "target_after_output_transform": transformed_target.tolist(),
        "target_consistency_dimension": PARAMS.outputs - coefficient_rank,
        "verifier_output": output.tolist(),
    }


def run_trial(key_label: bytes, message: bytes, rng_seed: int) -> dict:
    public_key, key_metadata = reduced_keygen(key_label)
    signature, solve = forge(
        public_key,
        message,
        np.random.default_rng(rng_seed),
    )
    altered_signature = bytearray(signature)
    altered_signature[0] ^= 1
    altered_salt = bytearray(signature)
    altered_salt[-1] ^= 1
    altered_public_key = bytearray(public_key)
    altered_public_key[0] ^= 1
    negative_tests = {
        "altered_message_rejected": not literal_verify(
            public_key, message + b"!", signature
        ),
        "altered_public_key_rejected": not literal_verify(
            bytes(altered_public_key), message, signature
        ),
        "altered_signature_rejected": not literal_verify(
            public_key, message, bytes(altered_signature)
        ),
        "altered_salt_rejected": not literal_verify(
            public_key, message, bytes(altered_salt)
        ),
    }
    if not literal_verify(public_key, message, signature):
        raise AssertionError("final positive verifier check failed")
    if not all(negative_tests.values()):
        raise AssertionError("a negative verifier control unexpectedly accepted")
    return {
        "digest_hex": hashlib.shake_256(message).digest(64).hex(),
        "key_generation": key_metadata,
        "message_hex": message.hex(),
        "negative_tests": negative_tests,
        "public_key_hex": public_key.hex(),
        "public_key_sha256": hashlib.sha256(public_key).hexdigest(),
        "signature_hex": signature.hex(),
        "signature_sha256": hashlib.sha256(signature).hexdigest(),
        "solve": solve,
        "verifier_accepts": True,
    }


def build_report() -> dict:
    main_trial = run_trial(DEFAULT_KEY_LABEL, DEFAULT_MESSAGE, DEFAULT_RNG_SEED)
    batch = []
    for trial in range(BATCH_TRIALS):
        record = run_trial(
            f"SNOVA reduced batch key {trial}".encode(),
            f"reduced SNOVA batch message {trial}".encode(),
            0x534E4F560000 + trial,
        )
        batch.append(
            {
                "public_key_sha256": record["public_key_sha256"],
                "restricted_coefficient_rank": record["solve"][
                    "restricted_coefficient_rank"
                ],
                "salt_counter": record["solve"]["salt_counter"],
                "signature_sha256": record["signature_sha256"],
                "slices_tested": record["solve"]["slices_tested"],
                "trial": trial,
                "verifier_accepts": record["verifier_accepts"],
            }
        )
    return {
        "attack_boundary": (
            "forge(public_key, message, public_rng) receives no secret key or "
            "planted witness"
        ),
        "batch_regression": {
            "records": batch,
            "successes": sum(record["verifier_accepts"] for record in batch),
            "trials": len(batch),
        },
        "implementation_boundary": (
            "KAT-anchored Python Version 2.3 transcription with distinct staged "
            "and literal public-map evaluators; not the upstream C verifier"
        ),
        "main_trial": main_trial,
        "official_parameter_set": False,
        "parameters": {
            "K": PARAMS.unordered,
            "M": PARAMS.outputs,
            "ell": PARAMS.l,
            "m1": PARAMS.m1,
            "n": PARAMS.n,
            "o": PARAMS.o,
            "q": ref.Q,
            "r": PARAMS.r,
            "v": PARAMS.v,
        },
        "production_size_cost_evidence": False,
        "scope": (
            "non-planted zero-offset end-to-end forgery at a deliberately "
            "reduced, unofficial shape; not production-feasibility evidence"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--check",
        action="store_true",
        help="recompute all nine trials and compare with committed JSON",
    )
    modes.add_argument(
        "--write",
        action="store_true",
        help="regenerate the committed deterministic JSON",
    )
    parser.add_argument("--expected", type=Path, default=DEFAULT_EXPECTED)
    args = parser.parse_args()

    output = json.dumps(build_report(), indent=2, sort_keys=True) + "\n"
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
                "reduced-forgery output differs from committed JSON:\n" + difference
            )
        print(f"OK: {args.expected.name} (main trial plus {BATCH_TRIALS} fresh-key regressions)")
    elif args.write:
        args.expected.parent.mkdir(parents=True, exist_ok=True)
        args.expected.write_text(output)
        print(f"WROTE: {args.expected.name}")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
