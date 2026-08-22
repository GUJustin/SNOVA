#!/usr/bin/env python3
"""End-to-end reduced-parameter SNOVA symmetry-quotient forgery.

This is a clean Python transcription of the odd-q, symmetric, l=r=2 path of
SNOVA's pinned reference implementation (commit
9da14981336ede257c41ef53cc069989051e8181), specialized to the reduced row
(v,o,q,l,r)=(2,1,19,2,2).

The executable demonstration preserves the reference implementation's:
  * q=19 field and public S matrix;
  * fixed-ABQ SHAKE expansion and invertibility adjustment;
  * symmetric public-seed expansion and hidden-UOV P22 key derivation;
  * public-key and signature base-19 byte serialization;
  * verifier public-map algebra, target hashing, and rejection rule.

Only the dimensions are reduced.  The attack routine receives only serialized
public-key bytes and a chosen digest:
  1. expand the public key from its reference-format byte string;
  2. impose the zero-offset common-column relation X_i=[u_i | 0];
  3. interpolate the complete restricted public verifier;
  4. compute its rank-3 symmetry quotient and one consistency equation;
  5. search salts until the official target satisfies consistency;
  6. search public 3-dimensional slices and enumerate 19^3 points;
  7. serialize the resulting signature and verify it with a separately written
     direct verifier that does not reuse the attack evaluator.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

Q = 19
L = 2
R = 2
V = 2
O = 1
N = V + O
M1 = (O * R + L - 1) // L
ALPHA = L * R + R
M = O * L * R
K_EXPECTED = M1 * (L * (L + 1) // 2)
PACK_GF = 15
PACK_BYTES = 8
SEED_LENGTH_PUBLIC = 16
SEED_LENGTH_PRIVATE = 32
BYTES_SALT = 16
S = np.array([[1, 2], [2, 15]], dtype=np.int64)
I2 = np.eye(2, dtype=np.int64)
S_POW = np.stack([I2, S]) % Q


def mod(a):
    return np.asarray(a, dtype=np.int64) % Q


def matmul(a, b):
    return (np.asarray(a, dtype=np.int64) @ np.asarray(b, dtype=np.int64)) % Q


def matdet2(a):
    a = np.asarray(a, dtype=np.int64)
    return int((a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0]) % Q)


def shake(data: bytes, n: int) -> bytes:
    return hashlib.shake_256(data).digest(n)


def bytes_gf(num: int) -> int:
    return (PACK_BYTES * num + PACK_GF - 1) // PACK_GF


def compress_gf(values: np.ndarray) -> bytes:
    """Exact q=19 translation of snova_ref.c:compress_gf."""
    values = mod(values).reshape(-1)
    out = bytearray(bytes_gf(len(values)))
    idx = 0
    out_idx = 0
    while idx < len(values):
        val = 0
        fact = 1
        count = 0
        while count < PACK_GF and idx < len(values):
            val += fact * int(values[idx])
            idx += 1
            count += 1
            fact *= Q
        written = 0
        while written < PACK_BYTES and out_idx < len(out):
            out[out_idx] = val & 0xFF
            out_idx += 1
            val >>= 8
            written += 1
    return bytes(out)


def expand_gf(data: bytes, num: int) -> tuple[np.ndarray, bool]:
    """Exact q=19 translation of snova_ref.c:expand_gf.

    Returns (digits, canonical), where canonical is false when unused high bits
    encode a nonzero value, as the reference parser would reject.
    """
    expected = bytes_gf(num)
    if len(data) != expected:
        return np.zeros(num, dtype=np.int64), False
    values = np.zeros(num, dtype=np.int64)
    byte_idx = 0
    out_idx = 0
    residue = 0
    while out_idx < num:
        val = 0
        count = 0
        while count < PACK_BYTES and byte_idx < expected:
            val ^= int(data[byte_idx]) << (8 * count)
            byte_idx += 1
            count += 1
        digit_count = 0
        while digit_count < PACK_GF and out_idx < num:
            values[out_idx] = val % Q
            val //= Q
            out_idx += 1
            digit_count += 1
        residue |= val
    return values, residue == 0


def gen_a_fqs(coeff):
    coeff = [int(x) % Q for x in coeff]
    if coeff[-1] == 0:
        coeff[-1] = (Q - (coeff[0] + (1 if coeff[0] == 0 else 0))) % Q
    return (coeff[0] * I2 + coeff[1] * S) % Q, np.array(coeff, dtype=np.int64)


def improve_invertible(orig):
    a = mod(orig).copy()
    if matdet2(a) == 0:
        for f in range(1, Q):
            a = (a + f * S) % Q
            if matdet2(a) != 0:
                break
    if matdet2(a) == 0:
        raise RuntimeError("failed to improve matrix to invertible")
    return a


def fixed_abq():
    total = O * ALPHA * (R * R + L * R + 2 * L)
    raw = np.frombuffer(shake(b"SNOVA_ABQ", total), dtype=np.uint8).astype(np.int64) % Q
    p = 0
    araw = raw[p : p + O * ALPHA * R * R].reshape(O, ALPHA, R, R); p += O * ALPHA * R * R
    braw = raw[p : p + O * ALPHA * L * R].reshape(O, ALPHA, R, L); p += O * ALPHA * L * R
    q1raw = raw[p : p + O * ALPHA * L].reshape(O, ALPHA, L); p += O * ALPHA * L
    q2raw = raw[p : p + O * ALPHA * L].reshape(O, ALPHA, L)
    am = np.empty_like(araw)
    bm = np.empty_like(braw)
    q1 = np.empty_like(q1raw)
    q2 = np.empty_like(q2raw)
    for i in range(O):
        for a in range(ALPHA):
            am[i, a] = improve_invertible(araw[i, a])
            bm[i, a] = improve_invertible(braw[i, a])
            _, q1[i, a] = gen_a_fqs(q1raw[i, a])
            _, q2[i, a] = gen_a_fqs(q2raw[i, a])
    return am, bm, q1, q2


@dataclass
class ExpandedPublicKey:
    pk_seed: bytes
    P: np.ndarray  # (M1,N,N,L,L)
    Am: np.ndarray # (O,ALPHA,R,R)
    Bm: np.ndarray # (O,ALPHA,R,L)
    q1: np.ndarray # (O,ALPHA,L)
    q2: np.ndarray # (O,ALPHA,L)


def expand_public_seed(pk_seed: bytes):
    num_gen_pub_gf = M1 * (V * (V + 1) // 2 + V * O) * L * L
    raw = np.frombuffer(shake(pk_seed, num_gen_pub_gf), dtype=np.uint8).astype(np.int64) % Q
    p11 = np.zeros((M1, V, V, L, L), dtype=np.int64)
    p12 = np.zeros((M1, V, O, L, L), dtype=np.int64)
    p21 = np.zeros((M1, O, V, L, L), dtype=np.int64)
    p = 0
    for mi in range(M1):
        for ni in range(V):
            for i in range(L):
                for j in range(i, L):
                    value = int(raw[p]); p += 1
                    p11[mi, ni, ni, i, j] = value
                    p11[mi, ni, ni, j, i] = value
            for nj in range(ni + 1, V):
                block = raw[p:p + L * L].reshape(L, L); p += L * L
                p11[mi, ni, nj] = block
                p11[mi, nj, ni] = block.T
            for oj in range(O):
                block = raw[p:p + L * L].reshape(L, L); p += L * L
                p12[mi, ni, oj] = block
                p21[mi, oj, ni] = block.T
    return p11 % Q, p12 % Q, p21 % Q, p, num_gen_pub_gf


def expand_t12(sk_seed: bytes):
    need = O * V * L
    coeff = []
    stream = shake(sk_seed, 4096)
    limit = (256 // Q) * Q
    for byte in stream:
        if byte < limit:
            coeff.append(byte % Q)
            if len(coeff) == need:
                break
    if len(coeff) != need:
        raise RuntimeError("insufficient rejection-sampled bytes")
    coeff = np.array(coeff, dtype=np.int64).reshape(V, O, L)
    t12 = np.zeros((V, O, L, L), dtype=np.int64)
    for i in range(V):
        for j in range(O):
            t12[i, j], coeff[i, j] = gen_a_fqs(coeff[i, j])
    return t12


def compress_p22(p22: np.ndarray) -> bytes:
    packed = []
    for mi in range(M1):
        for ni in range(O):
            for i in range(L):
                for j in range(i, L):
                    packed.append(int(p22[mi, ni, ni, i, j]))
            for nj in range(ni + 1, O):
                for i in range(L):
                    for j in range(L):
                        packed.append(int(p22[mi, ni, nj, i, j]))
    return compress_gf(np.array(packed, dtype=np.int64))


def expand_p22(data: bytes) -> np.ndarray:
    numgf_pk = M1 * O * L * (O * L + 1) // 2
    packed, canonical = expand_gf(data, numgf_pk)
    if not canonical:
        raise ValueError("noncanonical public-key encoding")
    p22 = np.zeros((M1, O, O, L, L), dtype=np.int64)
    pos = 0
    for mi in range(M1):
        for ni in range(O):
            for i in range(L):
                for j in range(i, L):
                    value = int(packed[pos]); pos += 1
                    p22[mi, ni, ni, i, j] = value
                    p22[mi, ni, ni, j, i] = value
            for nj in range(ni + 1, O):
                for i in range(L):
                    for j in range(L):
                        value = int(packed[pos]); pos += 1
                        p22[mi, ni, nj, i, j] = value
                        p22[mi, nj, ni, j, i] = value
    return p22


def keygen_bytes(seed: bytes) -> tuple[bytes, dict]:
    if len(seed) != SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE:
        raise ValueError("seed must be 48 bytes")
    pk_seed, sk_seed = seed[:SEED_LENGTH_PUBLIC], seed[SEED_LENGTH_PUBLIC:]
    p11, p12, p21, consumed, allocated = expand_public_seed(pk_seed)
    t12 = expand_t12(sk_seed)

    f12 = np.zeros((M1, V, O, L, L), dtype=np.int64)
    for mi in range(M1):
        for j in range(V):
            for k in range(O):
                for idx in range(V):
                    f12[mi, j, k] = (f12[mi, j, k] + matmul(p11[mi, j, idx], t12[idx, k])) % Q
                f12[mi, j, k] = (f12[mi, j, k] + p12[mi, j, k]) % Q

    p22 = np.zeros((M1, O, O, L, L), dtype=np.int64)
    for mi in range(M1):
        for j in range(O):
            for k in range(O):
                for idx in range(V):
                    p22[mi, j, k] = (p22[mi, j, k] + matmul(t12[idx, j], f12[mi, idx, k])) % Q
                    p22[mi, j, k] = (p22[mi, j, k] + matmul(p21[mi, j, idx], t12[idx, k])) % Q
                p22[mi, j, k] = (-p22[mi, j, k]) % Q

    pk_bytes = pk_seed + compress_p22(p22)
    metadata = {
        "key_seed_sha256": hashlib.sha256(seed).hexdigest(),
        "pk_seed_hex": pk_seed.hex(),
        "public_symbols_consumed": consumed,
        "public_symbols_allocated": allocated,
        "public_key_bytes": len(pk_bytes),
        "secret_material_passed_to_attack": False,
    }
    return pk_bytes, metadata


def expand_public_key_bytes(pk_bytes: bytes) -> ExpandedPublicKey:
    numgf_pk = M1 * O * L * (O * L + 1) // 2
    expected = SEED_LENGTH_PUBLIC + bytes_gf(numgf_pk)
    if len(pk_bytes) != expected:
        raise ValueError(f"public key has {len(pk_bytes)} bytes, expected {expected}")
    pk_seed = pk_bytes[:SEED_LENGTH_PUBLIC]
    p22 = expand_p22(pk_bytes[SEED_LENGTH_PUBLIC:])
    p11, p12, p21, _, _ = expand_public_seed(pk_seed)
    p = np.zeros((M1, N, N, L, L), dtype=np.int64)
    p[:, :V, :V] = p11
    p[:, :V, V:] = p12
    p[:, V:, :V] = p21
    p[:, V:, V:] = p22
    am, bm, q1, q2 = fixed_abq()
    return ExpandedPublicKey(pk_seed=pk_seed, P=p, Am=am, Bm=bm, q1=q1, q2=q2)


def attack_public_map(pk: ExpandedPublicKey, x: np.ndarray) -> np.ndarray:
    """Evaluator used by interpolation and search; follows the reference staging."""
    x = mod(x).reshape(N, L, R)
    whipped = np.zeros((L, N, L, R), dtype=np.int64)
    for a in range(L):
        for idx in range(N):
            whipped[a, idx] = matmul(S_POW[a], x[idx])

    sum_t1 = np.zeros((M1, L, L, R, R), dtype=np.int64)
    for mi in range(M1):
        sum_t0 = np.zeros((L, N, L, R), dtype=np.int64)
        for b in range(L):
            for ni in range(N):
                acc = np.zeros((L, R), dtype=np.int64)
                for nj in range(N):
                    acc = (acc + matmul(pk.P[mi, ni, nj], whipped[b, nj])) % Q
                sum_t0[b, ni] = acc
        for a in range(L):
            for b in range(L):
                acc = np.zeros((R, R), dtype=np.int64)
                for ni in range(N):
                    acc = (acc + matmul(whipped[a, ni].T, sum_t0[b, ni])) % Q
                sum_t1[mi, a, b] = acc

    out = np.zeros((O, R, L), dtype=np.int64)
    for mi in range(O):
        for alpha in range(ALPHA):
            mip = (mi + alpha) % M1
            temp1 = np.zeros((R, R), dtype=np.int64)
            for a in range(L):
                for b in range(L):
                    temp1 = (temp1 + int(pk.q1[mi, alpha, a]) * int(pk.q2[mi, alpha, b]) * sum_t1[mip, a, b]) % Q
            out[mi] = (out[mi] + matmul(pk.Am[mi, alpha], matmul(temp1, pk.Bm[mi, alpha]))) % Q
    return out.reshape(-1) % Q


def verifier_public_map_direct(pk: ExpandedPublicKey, x: np.ndarray) -> np.ndarray:
    """Independent direct evaluator for transcript verification.

    Unlike attack_public_map, this routine does not use the staged sum_t0/sum_t1
    decomposition.  It directly computes every X_i^T S_a^T P_ij S_b X_j term.
    """
    x = mod(x).reshape(N, L, R)
    out = np.zeros((O, R, L), dtype=np.int64)
    for mi in range(O):
        for alpha in range(ALPHA):
            mip = (mi + alpha) % M1
            temp1 = np.zeros((R, R), dtype=np.int64)
            for a in range(L):
                for b in range(L):
                    channel = np.zeros((R, R), dtype=np.int64)
                    for ni in range(N):
                        left = matmul(S_POW[a], x[ni])
                        for nj in range(N):
                            right = matmul(S_POW[b], x[nj])
                            channel = (channel + matmul(left.T, matmul(pk.P[mip, ni, nj], right))) % Q
                    temp1 = (temp1 + int(pk.q1[mi, alpha, a]) * int(pk.q2[mi, alpha, b]) * channel) % Q
            out[mi] = (out[mi] + matmul(pk.Am[mi, alpha], matmul(temp1, pk.Bm[mi, alpha]))) % Q
    return out.reshape(-1) % Q


def target_from_digest(pk_seed: bytes, digest: bytes, salt: bytes) -> np.ndarray:
    raw = shake(pk_seed + digest + salt, bytes_gf(M))
    target, _ = expand_gf(raw, M)  # Reference verifier ignores the leftover residue.
    return target


def signature_from_u(u: np.ndarray) -> np.ndarray:
    u = mod(u).reshape(N, L)
    x = np.zeros((N, L, R), dtype=np.int64)
    x[:, :, 0] = u
    return x


def serialize_signature(x: np.ndarray, salt: bytes) -> bytes:
    if len(salt) != BYTES_SALT:
        raise ValueError("salt must be 16 bytes")
    return compress_gf(mod(x).reshape(-1)) + salt


def parse_signature(sig_bytes: bytes) -> tuple[np.ndarray, bytes, bool]:
    encoded_len = bytes_gf(N * L * R)
    if len(sig_bytes) != encoded_len + BYTES_SALT:
        return np.zeros((N, L, R), dtype=np.int64), b"", False
    values, canonical = expand_gf(sig_bytes[:encoded_len], N * L * R)
    return values.reshape(N, L, R), sig_bytes[encoded_len:], canonical


def format_accepts(x: np.ndarray) -> bool:
    num_sym = 0
    for block in mod(x).reshape(N, L, R):
        num_sym += int(block[0, 1] == block[1, 0])
    return num_sym <= N // 4


def verify_serialized(pk_bytes: bytes, sig_bytes: bytes, digest: bytes) -> bool:
    """Reference-format verifier using only serialized public inputs."""
    try:
        pk = expand_public_key_bytes(pk_bytes)
    except ValueError:
        return False
    x, salt, canonical = parse_signature(sig_bytes)
    if not canonical or not format_accepts(x):
        return False
    target = target_from_digest(pk.pk_seed, digest, salt)
    return bool(np.array_equal(verifier_public_map_direct(pk, x), target))


def monomial_pairs(nvars: int):
    return [(i, j) for i in range(nvars) for j in range(i, nvars)]


def monomial_vector(x: np.ndarray, pairs) -> np.ndarray:
    x = mod(x).reshape(-1)
    return np.array([(int(x[i]) * int(x[j])) % Q for i, j in pairs], dtype=np.int64)


def interpolate_restricted_map(pk: ExpandedPublicKey):
    nvars = N * L
    pairs = monomial_pairs(nvars)
    coefficients = np.zeros((M, len(pairs)), dtype=np.int64)
    basis_values = []
    for i in range(nvars):
        e = np.zeros(nvars, dtype=np.int64); e[i] = 1
        basis_values.append(attack_public_map(pk, signature_from_u(e)))
    pair_to_col = {pair: col for col, pair in enumerate(pairs)}
    for i in range(nvars):
        coefficients[:, pair_to_col[(i, i)]] = basis_values[i]
    for i in range(nvars):
        for j in range(i + 1, nvars):
            e = np.zeros(nvars, dtype=np.int64); e[i] = 1; e[j] = 1
            value = attack_public_map(pk, signature_from_u(e))
            coefficients[:, pair_to_col[(i, j)]] = (value - basis_values[i] - basis_values[j]) % Q
    return coefficients, pairs


def inv_mod(a: int) -> int:
    a %= Q
    if a == 0:
        raise ZeroDivisionError
    return pow(a, Q - 2, Q)


def row_echelon_transform(a: np.ndarray):
    reduced = mod(a).copy()
    rows, cols = reduced.shape
    transform = np.eye(rows, dtype=np.int64)
    row = 0
    pivots = []
    for col in range(cols):
        pivot = next((r for r in range(row, rows) if reduced[r, col] % Q != 0), None)
        if pivot is None:
            continue
        if pivot != row:
            reduced[[row, pivot]] = reduced[[pivot, row]]
            transform[[row, pivot]] = transform[[pivot, row]]
        scale = inv_mod(int(reduced[row, col]))
        reduced[row] = (reduced[row] * scale) % Q
        transform[row] = (transform[row] * scale) % Q
        for r in range(rows):
            if r != row and reduced[r, col] % Q:
                factor = int(reduced[r, col])
                reduced[r] = (reduced[r] - factor * reduced[row]) % Q
                transform[r] = (transform[r] - factor * transform[row]) % Q
        pivots.append(col)
        row += 1
        if row == rows:
            break
    return reduced, transform % Q, row, pivots


def matrix_rank_mod(a: np.ndarray) -> int:
    return row_echelon_transform(a)[2]


def enumerate_points(dim: int) -> np.ndarray:
    return np.array(list(itertools.product(range(Q), repeat=dim)), dtype=np.int64)


def monomial_matrix(points: np.ndarray, pairs) -> np.ndarray:
    columns = [(points[:, i] * points[:, j]) % Q for i, j in pairs]
    return np.stack(columns, axis=1)


def attack(pk_bytes: bytes, digest: bytes, rng: np.random.Generator, max_salts=100000, max_slices=10000):
    """Forge from serialized public key bytes only."""
    pk = expand_public_key_bytes(pk_bytes)

    # Cross-check the attack evaluator against the independently organized
    # direct verifier on unrestricted signatures.
    for _ in range(64):
        x = rng.integers(0, Q, size=(N, L, R), dtype=np.int64)
        if not np.array_equal(attack_public_map(pk, x), verifier_public_map_direct(pk, x)):
            raise AssertionError("attack and direct verifier evaluators disagree")

    coefficients, pairs = interpolate_restricted_map(pk)
    reduced, output_transform, rank, pivots = row_echelon_transform(coefficients)
    if rank != K_EXPECTED:
        raise RuntimeError(f"quotient rank {rank}, expected {K_EXPECTED}")
    if np.any(reduced[rank:] % Q):
        raise AssertionError("row elimination failed")

    # Confirm the interpolated restricted map against both evaluators.
    for _ in range(64):
        u = rng.integers(0, Q, size=N * L, dtype=np.int64)
        interpolated = (coefficients @ monomial_vector(u, pairs)) % Q
        x = signature_from_u(u)
        if not np.array_equal(interpolated, attack_public_map(pk, x)):
            raise AssertionError("restricted interpolation mismatch")
        if not np.array_equal(interpolated, verifier_public_map_direct(pk, x)):
            raise AssertionError("restricted interpolation/direct verifier mismatch")

    salt = None
    target = None
    transformed_target = None
    salt_counter = None
    for counter in range(max_salts):
        candidate_salt = counter.to_bytes(BYTES_SALT, "little")
        candidate_target = target_from_digest(pk.pk_seed, digest, candidate_salt)
        candidate_transformed = (output_transform @ candidate_target) % Q
        if np.all(candidate_transformed[rank:] == 0):
            salt = candidate_salt
            target = candidate_target
            transformed_target = candidate_transformed
            salt_counter = counter
            break
    if salt is None:
        raise RuntimeError("no consistent target found")

    grid = enumerate_points(rank)
    solution = None
    chosen_slice = None
    chosen_point = None
    slices_tested = 0
    for _ in range(max_slices):
        slice_matrix = rng.integers(0, Q, size=(N * L, rank), dtype=np.int64)
        if matrix_rank_mod(slice_matrix) != rank:
            continue
        slices_tested += 1
        points = (grid @ slice_matrix.T) % Q
        features = monomial_matrix(points, pairs)
        values = (features @ coefficients.T) % Q
        mask = np.all(values == target[None, :], axis=1)
        # For X_i=[u_i|0], a block is symmetric iff the second row entry is 0.
        mask &= np.all(points[:, 1::2] != 0, axis=1)
        indices = np.flatnonzero(mask)
        if indices.size:
            index = int(indices[0])
            solution = points[index]
            chosen_slice = slice_matrix
            chosen_point = grid[index]
            break
    if solution is None:
        raise RuntimeError("no solution found on searched slices")

    x_signature = signature_from_u(solution)
    sig_bytes = serialize_signature(x_signature, salt)
    if not verify_serialized(pk_bytes, sig_bytes, digest):
        raise AssertionError("forged serialized signature failed independent verifier")

    output = verifier_public_map_direct(pk, x_signature)
    transformed_output = (output_transform @ output) % Q
    if not np.array_equal(output, target) or not np.array_equal(transformed_output, transformed_target):
        raise AssertionError("canonical output split mismatch")

    # Negative controls: changing a signature byte or salt must invalidate it.
    mutated_signature = bytearray(sig_bytes)
    mutated_signature[0] ^= 1
    if verify_serialized(pk_bytes, bytes(mutated_signature), digest):
        raise AssertionError("mutated signature unexpectedly verified")
    mutated_salt = bytearray(sig_bytes)
    mutated_salt[-1] ^= 1
    if verify_serialized(pk_bytes, bytes(mutated_salt), digest):
        raise AssertionError("mutated salt unexpectedly verified")

    return {
        "parameters": {"v": V, "o": O, "q": Q, "l": L, "r": R, "m1": M1, "n": N, "M": M, "K": rank},
        "digest_hex": digest.hex(),
        "public_key_hex": pk_bytes.hex(),
        "public_key_sha256": hashlib.sha256(pk_bytes).hexdigest(),
        "signature_hex": sig_bytes.hex(),
        "signature_sha256": hashlib.sha256(sig_bytes).hexdigest(),
        "salt_hex": salt.hex(),
        "salt_counter": salt_counter,
        "target": target.tolist(),
        "target_after_output_transform": transformed_target.tolist(),
        "quotient_rank": rank,
        "target_consistency_dimension": M - rank,
        "output_row_transform": output_transform.tolist(),
        "reduced_quadratic_coefficients": coefficients.tolist(),
        "monomial_pairs": pairs,
        "slices_tested": slices_tested,
        "slice_matrix": chosen_slice.tolist(),
        "slice_point": chosen_point.tolist(),
        "common_column_vector_u": solution.tolist(),
        "signature_matrices": x_signature.tolist(),
        "verifier_output": output.tolist(),
        "format_accepts": format_accepts(x_signature),
        "serialized_signature_is_canonical": parse_signature(sig_bytes)[2],
        "independent_verifier_accepts": True,
        "negative_control_signature_mutation_rejected": True,
        "negative_control_salt_mutation_rejected": True,
        "attack_input_was_serialized_public_key_only": True,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("reduced_attack_transcript.json"))
    parser.add_argument("--seed", default="SNOVA reduced attack key seed v2")
    parser.add_argument("--digest", default="chosen message for reduced SNOVA forgery")
    args = parser.parse_args()

    key_seed = shake(args.seed.encode(), SEED_LENGTH_PUBLIC + SEED_LENGTH_PRIVATE)
    pk_bytes, key_metadata = keygen_bytes(key_seed)
    # The attack is deliberately called with public bytes only.
    rng = np.random.default_rng(0x534E4F5641)
    transcript = attack(pk_bytes, args.digest.encode(), rng)
    transcript["key_generation"] = key_metadata
    args.out.write_text(json.dumps(transcript, indent=2) + "\n")
    print(json.dumps({
        "independent_verifier_accepts": transcript["independent_verifier_accepts"],
        "parameters": transcript["parameters"],
        "public_key_bytes": len(bytes.fromhex(transcript["public_key_hex"])),
        "signature_bytes": len(bytes.fromhex(transcript["signature_hex"])),
        "salt_counter": transcript["salt_counter"],
        "slices_tested": transcript["slices_tested"],
        "target": transcript["target"],
        "u": transcript["common_column_vector_u"],
        "signature": transcript["signature_matrices"],
        "transcript": str(args.out),
    }, indent=2))


if __name__ == "__main__":
    main()
