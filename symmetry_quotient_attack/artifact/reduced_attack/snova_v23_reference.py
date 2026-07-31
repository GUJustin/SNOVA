#!/usr/bin/env python3
"""Small, self-contained SNOVA Version-2.3 q=19 reference helpers.

The module provides finite-field linear algebra, pinned public-key and ABQ
expansion conventions, public-map evaluation, q=19 serialization, and the
common-column quotient decomposition used by the correspondence tests.  It has
no test-vector main routine and contains no planted target or known solution.
"""
from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


Q = 19
QA, QB, QC = 1, 3, 15


@dataclass(frozen=True)
class Params:
    name: str
    v: int
    o: int
    l: int
    r: int
    target: int = 0

    @property
    def n(self) -> int:
        return self.v + self.o

    @property
    def m1(self) -> int:
        return (self.o * self.r + self.l - 1) // self.l

    @property
    def alpha(self) -> int:
        return self.l * self.r + self.r

    @property
    def outputs(self) -> int:
        return self.o * self.r * self.l

    @property
    def variables(self) -> int:
        return self.n * self.l

    @property
    def unordered(self) -> int:
        return self.m1 * self.l * (self.l + 1) // 2


PUBLISHED_PARAMS = (
    Params("I-square-l4", 28, 5, 4, 4, 143),
    Params("I-square-l2", 48, 16, 2, 2, 143),
    Params("I-rect-l4xr5", 28, 4, 4, 5, 143),
    Params("III-square-l4", 40, 7, 4, 4, 207),
    Params("III-square-l2", 72, 24, 2, 2, 207),
    Params("III-rect-l4xr5", 38, 5, 4, 5, 207),
    Params("V-square-l4", 50, 9, 4, 4, 272),
    Params("V-square-l2", 96, 32, 2, 2, 272),
    Params("V-rect-l4xr6", 52, 6, 4, 6, 272),
)


@dataclass
class ABQ:
    A: np.ndarray
    B: np.ndarray
    q1: np.ndarray
    q2: np.ndarray
    S: np.ndarray


def inv_mod(value: int, modulus: int = Q) -> int:
    value = int(value) % modulus
    if value == 0:
        raise ZeroDivisionError
    return pow(value, modulus - 2, modulus)


def rref_mod(matrix: np.ndarray, modulus: int = Q):
    work = np.array(matrix, dtype=np.int64, copy=True) % modulus
    rows, columns = work.shape
    pivots: list[int] = []
    row = 0
    for column in range(columns):
        pivot = next(
            (candidate for candidate in range(row, rows) if work[candidate, column] % modulus),
            None,
        )
        if pivot is None:
            continue
        if pivot != row:
            work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * inv_mod(work[row, column], modulus) % modulus
        for other in range(rows):
            if other != row and work[other, column] % modulus:
                work[other] = (
                    work[other] - int(work[other, column]) * work[row]
                ) % modulus
        pivots.append(column)
        row += 1
        if row == rows:
            break
    return work, pivots


def rank_mod(matrix: np.ndarray, modulus: int = Q) -> int:
    return len(rref_mod(matrix, modulus)[1])


def nullspace_mod(matrix: np.ndarray, modulus: int = Q) -> np.ndarray:
    reduced, pivots = rref_mod(matrix, modulus)
    columns = matrix.shape[1]
    free = [column for column in range(columns) if column not in pivots]
    basis = []
    for free_column in free:
        vector = np.zeros(columns, dtype=np.int64)
        vector[free_column] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row, free_column] % modulus
        basis.append(vector)
    return np.asarray(basis, dtype=np.int64).reshape(len(basis), columns)


def determinant_mod(matrix: np.ndarray, modulus: int = Q) -> int:
    work = np.array(matrix, dtype=np.int64, copy=True) % modulus
    size = work.shape[0]
    determinant = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row, column]), None
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[[column, pivot]] = work[[pivot, column]]
            determinant = -determinant
        value = int(work[column, column])
        determinant = determinant * value % modulus
        work[column] = work[column] * inv_mod(value, modulus) % modulus
        for row in range(column + 1, size):
            if work[row, column]:
                work[row] = (
                    work[row] - int(work[row, column]) * work[column]
                ) % modulus
    return determinant % modulus


def affine_solve(coefficients: np.ndarray, target: np.ndarray, modulus: int = Q):
    coefficients = np.asarray(coefficients, dtype=np.int64) % modulus
    target = np.asarray(target, dtype=np.int64).reshape(-1, 1) % modulus
    variables = coefficients.shape[1]
    reduced, augmented_pivots = rref_mod(
        np.concatenate([coefficients, target], axis=1), modulus
    )
    if variables in augmented_pivots:
        raise ValueError("inconsistent affine system")
    pivots = [pivot for pivot in augmented_pivots if pivot < variables]
    free = [column for column in range(variables) if column not in pivots]
    base = np.zeros(variables, dtype=np.int64)
    kernel = np.zeros((variables, len(free)), dtype=np.int64)
    for index, free_column in enumerate(free):
        kernel[free_column, index] = 1
    for row, pivot in enumerate(pivots):
        base[pivot] = reduced[row, variables]
        for index, free_column in enumerate(free):
            kernel[pivot, index] = -reduced[row, free_column] % modulus
    if not np.array_equal(coefficients @ base % modulus, target[:, 0]):
        raise AssertionError("affine base point verification failed")
    if np.any(coefficients @ kernel % modulus):
        raise AssertionError("affine kernel verification failed")
    return base, kernel, free, pivots


def official_S(ell: int) -> np.ndarray:
    matrix = np.empty((ell, ell), dtype=np.int64)
    for row in range(ell):
        for column in range(ell):
            matrix[row, column] = (QA + row + column) & QB
    matrix[-1, -1] = QC
    return matrix % Q


def powers(matrix: np.ndarray, count: int) -> list[np.ndarray]:
    result = [np.eye(matrix.shape[0], dtype=np.int64)]
    for _ in range(1, count):
        result.append(result[-1] @ matrix % Q)
    return result


def block_power(matrix: np.ndarray, blocks: int) -> np.ndarray:
    return np.kron(np.eye(blocks, dtype=np.int64), matrix) % Q


def snova_xof(seed: bytes, count: int) -> bytes:
    result = bytearray()
    blocks = (count + 167) // 168
    for block in range(blocks):
        result += hashlib.shake_128(
            seed + block.to_bytes(8, "little")
        ).digest(168)
    return bytes(result[:count])


def field_S_element(coefficients: Sequence[int], S_powers: Sequence[np.ndarray]) -> np.ndarray:
    values = [int(value) % Q for value in coefficients]
    if values[-1] == 0:
        values[-1] = Q - (values[0] if values[0] else 1)
    return sum(
        (values[index] * S_powers[index] for index in range(len(values))),
        start=np.zeros_like(S_powers[0]),
    ) % Q


def pinned_improve(matrix: np.ndarray, S: np.ndarray, ell: int) -> np.ndarray:
    """Pinned reference-C cumulative `be_invertible_by_add_aS` behavior."""
    result = np.array(matrix, dtype=np.int64, copy=True) % Q
    if result.shape != (ell, ell) or determinant_mod(result):
        return result
    for coefficient in range(1, Q):
        result = (result + coefficient * S) % Q
        if determinant_mod(result):
            return result
    raise RuntimeError("pinned cumulative matrix improvement failed")


def _repair_field_S_coefficients(values: np.ndarray) -> np.ndarray:
    result = np.array(values, dtype=np.int64, copy=True) % Q
    for row in result.reshape(-1, result.shape[-1]):
        if row[-1] == 0:
            row[-1] = Q - (int(row[0]) if row[0] else 1)
    return result


def reconstruct_fixed_abq(params: Params) -> ABQ:
    S = official_S(params.l)
    count = params.o * params.alpha * (
        params.r * params.r + params.r * params.l + 2 * params.l
    )
    raw = np.frombuffer(
        hashlib.shake_256(b"SNOVA_ABQ").digest(count), dtype=np.uint8
    ).astype(np.int64) % Q
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
    position += size
    size = params.o * params.alpha * params.l
    q1 = raw[position:position + size].reshape(params.o, params.alpha, params.l)
    position += size
    q2 = raw[position:position + size].reshape(params.o, params.alpha, params.l)
    position += size
    if position != count:
        raise AssertionError("fixed ABQ expansion length mismatch")
    A = np.empty_like(raw_A)
    B = np.empty_like(raw_B)
    for output in range(params.o):
        for alpha in range(params.alpha):
            A[output, alpha] = pinned_improve(raw_A[output, alpha], S, params.l)
            B[output, alpha] = pinned_improve(raw_B[output, alpha], S, params.l)
    return ABQ(
        A,
        B,
        _repair_field_S_coefficients(q1),
        _repair_field_S_coefficients(q2),
        S,
    )


def encoded_length(number_of_digits: int) -> int:
    return math.ceil(8 * number_of_digits / 15)


def decode_base19(
    data: bytes, number_of_digits: int, *, require_canonical: bool = True
) -> np.ndarray:
    expected = encoded_length(number_of_digits)
    if len(data) != expected:
        raise ValueError(f"expected {expected} bytes, received {len(data)}")
    output: list[int] = []
    offset = 0
    remaining = number_of_digits
    while remaining:
        take_digits = min(15, remaining)
        take_bytes = min(8, len(data) - offset)
        value = int.from_bytes(data[offset:offset + take_bytes], "little")
        if require_canonical and value >= Q**take_digits:
            raise ValueError("noncanonical q=19 byte encoding")
        output.extend((value // Q**index) % Q for index in range(take_digits))
        offset += take_bytes
        remaining -= take_digits
    if offset != len(data):
        raise AssertionError("decoder did not consume the complete input")
    return np.asarray(output, dtype=np.int64)


def encode_base19(digits: np.ndarray) -> bytes:
    values = [int(value) % Q for value in np.asarray(digits).reshape(-1)]
    output = bytearray()
    for offset in range(0, len(values), 15):
        chunk = values[offset:offset + 15]
        integer = sum(value * Q**index for index, value in enumerate(chunk))
        output.extend(integer.to_bytes(8, "little"))
    return bytes(output[:encoded_length(len(values))])


def public_random_matrices(public_seed: bytes, params: Params):
    ell = params.l
    count = params.m1 * (
        params.v * ell * (ell + 1) // 2
        + (params.v * (params.v - 1) // 2 + params.v * params.o) * ell * ell
    )
    data = np.frombuffer(
        snova_xof(public_seed, count), dtype=np.uint8
    ).astype(np.int64) % Q
    position = 0
    triples = []
    for _ in range(params.m1):
        p11 = np.zeros((params.v * ell, params.v * ell), dtype=np.int64)
        p12 = np.zeros((params.v * ell, params.o * ell), dtype=np.int64)
        for block_i in range(params.v):
            for row in range(ell):
                for column in range(row, ell):
                    value = int(data[position])
                    position += 1
                    p11[block_i * ell + row, block_i * ell + column] = value
                    p11[block_i * ell + column, block_i * ell + row] = value
            for block_j in range(block_i + 1, params.v):
                for row in range(ell):
                    for column in range(ell):
                        value = int(data[position])
                        position += 1
                        p11[block_i * ell + row, block_j * ell + column] = value
                        p11[block_j * ell + column, block_i * ell + row] = value
            for oil_j in range(params.o):
                for row in range(ell):
                    for column in range(ell):
                        p12[block_i * ell + row, oil_j * ell + column] = int(
                            data[position]
                        )
                        position += 1
        triples.append((p11, p12, p12.T.copy()))
    if position != count:
        raise AssertionError("public XOF expansion length mismatch")
    return triples


def unpack_p22(packed: bytes, params: Params) -> np.ndarray:
    scalar_dimension = params.o * params.l
    digit_count = params.m1 * scalar_dimension * (scalar_dimension + 1) // 2
    values = decode_base19(packed, digit_count)
    position = 0
    blocks = []
    for _ in range(params.m1):
        matrix = np.zeros((scalar_dimension, scalar_dimension), dtype=np.int64)
        for block_i in range(params.o):
            for row in range(params.l):
                for column in range(row, params.l):
                    value = int(values[position])
                    position += 1
                    left = block_i * params.l + row
                    right = block_i * params.l + column
                    matrix[left, right] = matrix[right, left] = value
                for block_j in range(block_i + 1, params.o):
                    for column in range(params.l):
                        value = int(values[position])
                        position += 1
                        left = block_i * params.l + row
                        right = block_j * params.l + column
                        matrix[left, right] = matrix[right, left] = value
        blocks.append(matrix)
    if position != digit_count:
        raise AssertionError("P22 parser did not consume every coefficient")
    return np.asarray(blocks, dtype=np.int64)


def pack_p22(blocks: np.ndarray, params: Params) -> bytes:
    values: list[int] = []
    for matrix in blocks:
        for block_i in range(params.o):
            for row in range(params.l):
                for column in range(row, params.l):
                    values.append(
                        int(matrix[block_i * params.l + row, block_i * params.l + column])
                    )
                for block_j in range(block_i + 1, params.o):
                    for column in range(params.l):
                        values.append(
                            int(matrix[block_i * params.l + row, block_j * params.l + column])
                        )
    return encode_base19(np.asarray(values, dtype=np.int64))


def expand_public_key(public_key: bytes, params: Params) -> np.ndarray:
    random_parts = public_random_matrices(public_key[:16], params)
    p22 = unpack_p22(public_key[16:], params)
    matrices = []
    for (p11, p12, p21), block22 in zip(random_parts, p22):
        matrices.append(np.block([[p11, p12], [p21, block22]]) % Q)
    result = np.asarray(matrices, dtype=np.int64)
    if any(not np.array_equal(matrix, matrix.T) for matrix in result):
        raise AssertionError("odd-characteristic public matrices must be symmetric")
    return result


def expand_t12(secret_seed: bytes, params: Params) -> np.ndarray:
    needed = params.o * params.v * params.l
    stream = hashlib.shake_256(secret_seed).digest(max(64, 4 * needed))
    accepted = [
        byte % Q for byte in stream if byte < (256 // Q) * Q
    ]
    if len(accepted) < needed:
        raise RuntimeError("increase deterministic private expansion buffer")
    coefficients = accepted[:needed]
    S_powers = powers(official_S(params.l), params.l)
    result = np.zeros(
        (params.v * params.l, params.o * params.l), dtype=np.int64
    )
    for vinegar in range(params.v):
        for oil in range(params.o):
            offset = (vinegar * params.o + oil) * params.l
            result[
                vinegar * params.l:(vinegar + 1) * params.l,
                oil * params.l:(oil + 1) * params.l,
            ] = field_S_element(
                coefficients[offset:offset + params.l], S_powers
            )
    return result


def derive_public_key_from_seeds(
    public_seed: bytes, secret_seed: bytes, params: Params
) -> tuple[bytes, np.ndarray]:
    T12 = expand_t12(secret_seed, params)
    random_parts = public_random_matrices(public_seed, params)
    p22_blocks = []
    matrices = []
    for p11, p12, p21 in random_parts:
        p22 = -(T12.T @ ((p11 @ T12 + p12) % Q) + p21 @ T12) % Q
        p22_blocks.append(p22)
        matrices.append(np.block([[p11, p12], [p21, p22]]) % Q)
    public_key = public_seed + pack_p22(
        np.asarray(p22_blocks, dtype=np.int64), params
    )
    return public_key, np.asarray(matrices, dtype=np.int64)


def build_quotient(params: Params, abq: ABQ, rho: Sequence[int]):
    rho = np.asarray(rho, dtype=np.int64) % Q
    rho_outer = np.outer(rho, rho) % Q
    triangular = params.l * (params.l + 1) // 2
    quotient = np.zeros((params.outputs, params.m1 * triangular), dtype=np.int64)
    labels = [
        (public_index, left_power, right_power)
        for public_index in range(params.m1)
        for left_power in range(params.l)
        for right_power in range(left_power, params.l)
    ]
    for output in range(params.o):
        rows = slice(
            output * params.r * params.l,
            (output + 1) * params.r * params.l,
        )
        for alpha in range(params.alpha):
            public_index = (output + alpha) % params.m1
            core = abq.A[output, alpha] @ rho_outer @ abq.B[output, alpha] % Q
            coordinate = 0
            for left_power in range(params.l):
                for right_power in range(left_power, params.l):
                    coefficient = (
                        int(abq.q1[output, alpha, left_power])
                        * int(abq.q2[output, alpha, right_power])
                    )
                    if left_power != right_power:
                        coefficient += (
                            int(abq.q1[output, alpha, right_power])
                            * int(abq.q2[output, alpha, left_power])
                        )
                    column = public_index * triangular + coordinate
                    quotient[rows, column] = (
                        quotient[rows, column]
                        + (coefficient % Q) * core.reshape(-1)
                    ) % Q
                    coordinate += 1
    return quotient, labels


def quotient_coordinates(
    params: Params, public_matrices: np.ndarray, x: np.ndarray, labels
) -> np.ndarray:
    S_powers = powers(official_S(params.l), params.l)
    whipped = [block_power(power, params.n) @ x % Q for power in S_powers]
    output = []
    for public_index, left_power, right_power in labels:
        output.append(
            int(
                whipped[left_power].T
                @ public_matrices[public_index]
                @ whipped[right_power]
            ) % Q
        )
    return np.asarray(output, dtype=np.int64)


def common_column_signature(
    x: np.ndarray, rho: np.ndarray, offset: np.ndarray
) -> np.ndarray:
    return (x[:, None] * rho[None, :] + offset) % Q


def attack_public_output(
    params: Params, abq: ABQ, public_matrices: np.ndarray, signature: np.ndarray
) -> np.ndarray:
    S_powers = powers(abq.S, params.l)
    whipped = [block_power(power, params.n) @ signature % Q for power in S_powers]
    native = np.empty(
        (params.m1, params.l, params.l, params.r, params.r),
        dtype=np.int64,
    )
    for public_index in range(params.m1):
        for left_power in range(params.l):
            left = whipped[left_power].T @ public_matrices[public_index] % Q
            for right_power in range(params.l):
                native[public_index, left_power, right_power] = (
                    left @ whipped[right_power] % Q
                )
    result = np.zeros((params.o, params.r, params.l), dtype=np.int64)
    for output in range(params.o):
        for alpha in range(params.alpha):
            public_index = (output + alpha) % params.m1
            middle = np.zeros((params.r, params.r), dtype=np.int64)
            for left_power in range(params.l):
                for right_power in range(params.l):
                    middle = (
                        middle
                        + int(abq.q1[output, alpha, left_power])
                        * int(abq.q2[output, alpha, right_power])
                        * native[public_index, left_power, right_power]
                    ) % Q
            result[output] = (
                result[output]
                + abq.A[output, alpha] @ middle @ abq.B[output, alpha]
            ) % Q
    return result.reshape(-1)


def literal_public_output(
    params: Params, abq: ABQ, public_matrices: np.ndarray, signature: np.ndarray
) -> np.ndarray:
    S_powers = powers(abq.S, params.l)
    whipped = [block_power(power, params.n) @ signature % Q for power in S_powers]
    result = np.zeros((params.o, params.r, params.l), dtype=np.int64)
    for output in range(params.o):
        for alpha in range(params.alpha):
            public_index = (output + alpha) % params.m1
            middle = np.zeros((params.r, params.r), dtype=np.int64)
            for left_power in range(params.l):
                for right_power in range(params.l):
                    coefficient = (
                        int(abq.q1[output, alpha, left_power])
                        * int(abq.q2[output, alpha, right_power])
                    ) % Q
                    native = (
                        whipped[left_power].T
                        @ public_matrices[public_index]
                        @ whipped[right_power]
                    ) % Q
                    middle = (middle + coefficient * native) % Q
            result[output] = (
                result[output]
                + abq.A[output, alpha] @ middle @ abq.B[output, alpha]
            ) % Q
    return result.reshape(-1)


def _skew_vector(matrix: np.ndarray) -> np.ndarray:
    ell = matrix.shape[0]
    return np.asarray(
        [
            (int(matrix[row, column]) - int(matrix[column, row])) % Q
            for row in range(ell)
            for column in range(row + 1, ell)
        ],
        dtype=np.int64,
    )


def safe_offsets(params: Params, rho: np.ndarray, rng):
    offset = rng.integers(
        0, Q, size=(params.variables, params.r), dtype=np.int64
    )
    if params.l != params.r or params.l <= 2:
        return offset, {"deterministic": False}
    image = []
    for index in range(params.l):
        basis = np.zeros(params.l, dtype=np.int64)
        basis[index] = 1
        image.append(_skew_vector(np.outer(basis, rho)))
    image = np.asarray(image, dtype=np.int64).T
    image_rank = rank_mod(image)
    outside = None
    for index in range(image.shape[0]):
        basis = np.zeros(image.shape[0], dtype=np.int64)
        basis[index] = 1
        if rank_mod(np.column_stack([image, basis])) > image_rank:
            outside = basis
            break
    if outside is None:
        raise RuntimeError("failed to find alternating offset outside relation image")
    block_offset = np.zeros((params.l, params.l), dtype=np.int64)
    position = 0
    for row in range(params.l):
        for column in range(row + 1, params.l):
            block_offset[row, column] = outside[position]
            position += 1
    for block in range(params.n):
        offset[block * params.l:(block + 1) * params.l] = (
            offset[block * params.l:(block + 1) * params.l] + block_offset
        ) % Q
    return offset, {
        "deterministic": True,
        "skew_image_rank": image_rank,
        "alternating_dimension": image.shape[0],
    }


def rejection_count(signature: np.ndarray, params: Params) -> int:
    if params.l != params.r:
        return 0
    blocks = signature.reshape(params.n, params.l, params.r)
    return sum(np.array_equal(block, block.T) for block in blocks)


def message_target(
    public_seed: bytes, message: bytes, salt: bytes, outputs: int
) -> np.ndarray:
    digest = hashlib.shake_256(message).digest(64)
    packed = hashlib.shake_256(public_seed + digest + salt).digest(
        encoded_length(outputs)
    )
    return decode_base19(packed, outputs, require_canonical=False)


def parse_public_kat_record(path: Path, record: int = 0) -> dict[str, str]:
    """Parse only public verification fields from a signer KAT record."""
    current: dict[str, str] = {}
    current_index = -1
    for line in path.read_text().splitlines():
        if " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        key = key.strip()
        if key == "count":
            current_index += 1
            if current_index > record:
                break
            current = {}
        if current_index == record and key != "sk":
            current[key] = value.strip()
    required = {"count", "msg", "pk", "sm", "mlen"}
    if not required.issubset(current):
        raise ValueError(f"KAT record {record} is incomplete in {path}")
    return current
