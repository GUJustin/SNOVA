#!/usr/bin/env python3
"""Separately coded deterministic consistency checks for the SNOVA artifact.

The checker deliberately imports none of the generator code.  It transcribes
the formulas used in the paper, recomputes every ledger row, and only then
compares those values with ``all_nine_ledger.json``.  It also checks the scalar
multiplier exhaustively, the field-tower identities, target sampling,
dimension inequalities, spectrum tails, separator optima, and repair floors.
It does not prove the cryptographic reduction, the idealized XOF-transcript
model, the inherited homotopy hypothesis H_hom, or the size of kappa_hom.
"""
from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
Q = 19
RHO = Fraction(14, 256)
A2 = Q**2
A4 = Q**4
G2 = 692
G4 = 2628
REF = {"I": 143.0, "III": 207.0, "V": 272.0}
ROWS = (
    ("I", (28, 5, 4, 4)),
    ("I", (48, 16, 2, 2)),
    ("I", (28, 4, 4, 5)),
    ("III", (40, 7, 4, 4)),
    ("III", (72, 24, 2, 2)),
    ("III", (38, 5, 4, 5)),
    ("V", (50, 9, 4, 4)),
    ("V", (96, 32, 2, 2)),
    ("V", (52, 6, 4, 6)),
)


def log2_fraction(value: Fraction | int) -> float:
    if isinstance(value, int):
        return math.log2(value)
    return math.log2(value.numerator) - math.log2(value.denominator)


def logadd2(left: float, right: float) -> float:
    high, low = max(left, right), min(left, right)
    return high + math.log2(1.0 + 2.0 ** (low - high))


def eta128(h: int, K: int) -> Fraction:
    return Fraction(Q**h - 1, Q**K) + Fraction(1, 2**128)


def spectrum_tail(d: int, v: int, h: int, K: int, zero_offset: bool) -> Fraction:
    t = d * (v - 1 if zero_offset else v - 2)
    return Fraction(Q**h - 1) * (1 - Fraction(1, Q**K)) * RHO**t * 2**128


def atom_product(exponents: range) -> Fraction:
    value = Fraction(1)
    for exponent in exponents:
        value *= 1 - RHO**exponent
    return value


def separator(base_size: int, B: int, degree: int, gate_cost: int) -> tuple[Fraction, int]:
    candidates = []
    for ext_degree in range(1, 300):
        field_size = base_size**ext_degree
        if field_size < degree or field_size <= B * B:
            continue
        if 2 * ext_degree - 1 > base_size:
            continue
        success = 1 - Fraction(B * B, field_size)
        factor = Fraction((2 * ext_degree - 1) * gate_cost, 1) / success
        candidates.append((factor, ext_degree))
    return min(candidates)


def l4_acceptance_lower_bound(n: int, public_columns: int) -> Fraction:
    return (1 - Fraction(1, Q**4)) ** n if public_columns == 4 else Fraction(1)


def l2_acceptance(n: int) -> Fraction:
    numerator = sum(math.comb(n, j) * 18 ** (n - j) for j in range(n // 4 + 1))
    return Fraction(numerator, Q**n)


def accepted_root(alpha: Fraction, h: int, K: int, eta: Fraction) -> Fraction:
    a = alpha - eta / (Q - 1)
    return Fraction(1, Q ** (K - h)) * a * a / (a + eta)


def dense_quadratic_slp(nvars: int, neqs: int) -> int:
    monomials = math.comb(nvars + 1, 2)
    return monomials + 2 * neqs * (monomials + nvars + 1)


def l4_direct(level: str, params: tuple[int, int, int, int]) -> dict[str, object]:
    v, o, d, public_columns = params
    m1 = math.ceil(o * public_columns / d)
    K = m1 * math.comb(d + 1, 2)
    eta = eta128(K, K)
    root = accepted_root(l4_acceptance_lower_bound(v + o, public_columns), K, K, eta)
    B = 2**K
    Bplus = B * (K + 2) // 2
    slp = dense_quadratic_slp(K, K)
    H = slp + 2 * K + K * K
    factor, ext_degree = separator(A2, B, 2 * K, G2)
    work = factor * B * Bplus * H * K / root
    tail = spectrum_tail(d, v, K, K, False)
    normalized = work / (1 - tail)
    return {
        "direct": log2_fraction(normalized),
        "per_good": log2_fraction(work),
        "tail_log2": log2_fraction(tail),
        "extension": ext_degree,
    }


def l4_fast(level: str, params: tuple[int, int, int, int]) -> dict[str, object]:
    v, o, d, public_columns = params
    m1 = math.ceil(o * public_columns / d)
    K = m1 * math.comb(d + 1, 2)
    max_block_dimension = v + o - d * m1
    alpha = l4_acceptance_lower_bound(v + o, public_columns)
    candidates = []
    for channel_2 in range(m1 + 1):
        for channel_20 in range(m1 + 1):
            s = channel_2 + channel_20
            if not s or s > max_block_dimension or 4 * s > K:
                continue
            h = 4 * s
            eta = eta128(h, K)
            good = alpha - Fraction(channel_2 + Q * channel_20, A4)
            root = Fraction(1, Q ** (K - h)) * good * good / (good + eta)
            B = 2**channel_2 * 20**channel_20
            Bplus = Fraction(B) * (
                1 + Fraction(channel_2, 2) + Fraction(channel_20, 20)
            )
            assert Bplus.denominator == 1
            monomials = math.comb(s + 1, 2)
            slp_2 = (
                0
                if channel_2 == 0
                else monomials + 2 * channel_2 * (monomials + s + 1)
            )
            slp_20 = (
                0
                if channel_20 == 0
                else 6 * s + s * s + 2 * channel_20 * (s * s + 2 * s + 1)
            )
            degree = 2 * channel_2 + 20 * channel_20
            H = slp_2 + slp_20 + degree + s * s
            factor, ext_degree = separator(A4, B, degree, G4)
            work = factor * B * int(Bplus) * H * s / root
            projective_lines = Fraction(A4**s - 1, A4 - 1)
            jacobian_failure = projective_lines * RHO ** (4 * s)
            tail = spectrum_tail(d, v, h, K, False)
            normalized = work / (1 - jacobian_failure - tail)
            candidates.append(
                (
                    log2_fraction(normalized),
                    (s, channel_2, channel_20),
                    log2_fraction(work),
                    float(jacobian_failure),
                    float(tail),
                    ext_degree,
                )
            )
    _, profile, per_good, jacobian_failure, tail, ext_degree = min(candidates)
    return {
        "profile": profile,
        "per_good": per_good,
        "jacobian_failure": jacobian_failure,
        "tail": tail,
        "extension": ext_degree,
    }


def l4_adaptive(direct: dict[str, object], fast: dict[str, object]) -> float:
    tail = min(0.5, float(fast["tail"]) + 2.0 ** float(direct["tail_log2"]))
    fallback_probability = min(1.0, float(fast["jacobian_failure"]) / (1 - tail))
    fast_work = 2.0 ** float(fast["per_good"])
    fallback_work = 2.0 ** float(direct["per_good"])
    mixed = fast_work + fallback_probability * max(0.0, fallback_work - fast_work)
    return math.log2(mixed / (1 - tail))


def l2_complete(level: str, params: tuple[int, int, int, int]) -> dict[str, object]:
    v, m, d, _public_columns = params
    M = 4 * m
    K = 3 * m
    s = K // 2
    eta = eta128(K, K)
    root = accepted_root(l2_acceptance(v + m), K, K, eta)
    B = 2 ** (2 * m) * math.comb(m, m // 2)
    Bplus = Fraction(B) * (1 + 2 * m + Fraction(m * m, m + 2))
    assert Bplus.denominator == 1
    monomials = math.comb(s + 1, 2)
    slp = (
        2 * monomials
        + s * s
        + 4 * m * (monomials + s + 1)
        + 2 * m * (s + 1) ** 2
    )
    H = slp + 6 * m + K * K
    factor, ext_degree = separator(A2, B, 2 * s, G2)
    work = factor * B * int(Bplus) * H * K / root
    tail = spectrum_tail(d, v, K, K, True)
    solve = work / (1 - tail)
    target_log = (M - K) * math.log2(Q) + 32
    return {
        "direct": logadd2(log2_fraction(solve), target_log),
        "per_good": log2_fraction(work),
        "tail_log2": log2_fraction(tail),
        "target": target_log,
        "extension": ext_degree,
    }


def l2_diagonal(params: tuple[int, int, int, int]) -> dict[str, float]:
    v, m, d, _public_columns = params
    M = 4 * m
    K = 3 * m
    s = m
    h = 2 * s
    alpha = l2_acceptance(v + m)
    eta = eta128(h, K)
    good = alpha - Fraction(s, A2)
    root = Fraction(1, Q ** (K - h)) * max(good - eta, good * good / (good + eta))
    B = 2**s
    Bplus = B + s * 2 ** (s - 1)
    H = dense_quadratic_slp(s, s) + 2 * s + s * s
    factor, _ext_degree = separator(A2, B, 2 * s, G2)
    work = factor * B * Bplus * H * s / root
    jacobian_success = atom_product(range(2, 2 * s + 1, 2))
    tail = spectrum_tail(d, v, h, K, True)
    solve = work / (jacobian_success - tail)
    target_log = (M - K) * math.log2(Q) + 32
    return {
        "per_good": log2_fraction(work),
        "total": logadd2(log2_fraction(solve), target_log),
        "jacobian_failure": float(1 - jacobian_success),
        "tail": float(tail),
        "target": target_log,
    }


def l2_adaptive(
    level: str, complete: dict[str, object], diagonal: dict[str, float]
) -> float:
    if level != "I":
        return float(complete["direct"])
    tail = min(0.5, diagonal["tail"] + 2.0 ** float(complete["tail_log2"]))
    fallback_probability = min(1.0, diagonal["jacobian_failure"] / (1 - tail))
    fast_work = 2.0 ** diagonal["per_good"]
    fallback_work = 2.0 ** float(complete["per_good"])
    mixed = fast_work + fallback_probability * max(0.0, fallback_work - fast_work)
    return logadd2(math.log2(mixed / (1 - tail)), diagonal["target"])


def eval_netlist(net: dict, left: int, right: int) -> int:
    values = [bool((left >> bit) & 1) for bit in range(5)]
    values += [bool((right >> bit) & 1) for bit in range(5)]
    for operation, x, y in net["gates"]:
        xv = values[x[1]] if x[0] == "w" else bool(x[1])
        yv = values[y[1]] if y[0] == "w" else bool(y[1])
        if operation == "AND":
            values.append(xv and yv)
        elif operation == "XOR":
            values.append(xv ^ yv)
        elif operation == "XNOR":
            values.append(not (xv ^ yv))
        else:
            raise ValueError(operation)
    return sum(1 << bit for bit, wire in enumerate(net["outputs"]) if values[wire])


def poly_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % Q for x, y in zip(left, right))


def poly_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    raw = [0] * 7
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            raw[i + j] = (raw[i + j] + x * y) % Q
    for degree in range(6, 3, -1):
        coefficient = raw[degree] % Q
        raw[degree] = 0
        raw[degree - 3] = (raw[degree - 3] + coefficient) % Q
        raw[degree - 4] = (raw[degree - 4] + coefficient) % Q
    return tuple(raw[:4])


def determinant_mod_19(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column] % Q)
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            determinant = -determinant
        value = work[column][column] % Q
        determinant = determinant * value % Q
        inverse = pow(value, -1, Q)
        for j in range(column, len(work)):
            work[column][j] = work[column][j] * inverse % Q
        for row in range(column + 1, len(work)):
            coefficient = work[row][column] % Q
            for j in range(column, len(work)):
                work[row][j] = (work[row][j] - coefficient * work[column][j]) % Q
    return determinant % Q


def quadratic_factor_exists() -> bool:
    polynomial = (Q - 1, Q - 1, 0, 0, 1)
    for linear in range(Q):
        for constant in range(Q):
            remainder = list(polynomial)
            for degree in (4, 3, 2):
                coefficient = remainder[degree] % Q
                remainder[degree] = 0
                remainder[degree - 1] = (
                    remainder[degree - 1] - coefficient * linear
                ) % Q
                remainder[degree - 2] = (
                    remainder[degree - 2] - coefficient * constant
                ) % Q
            if remainder[0] == remainder[1] == 0:
                return True
    return False


def circuit_checks() -> None:
    net = json.loads((HERE / "f19_multiplier_netlist.json").read_text())
    payload = json.dumps(
        {"ninputs": net["ninputs"], "gates": net["gates"], "outputs": net["outputs"]},
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    assert hashlib.sha256(payload).hexdigest() == net["payload_sha256"]
    assert len(net["gates"]) == net["gate_count"] == 150
    for left in range(Q):
        for right in range(Q):
            assert eval_netlist(net, left, right) == left * right % Q
    assert all((x**4 - x - 1) % Q for x in range(Q))
    assert not quadratic_factor_exists()
    one = (1, 0, 0, 0)
    u = (1, 2, 15, 5)
    v = (4, 13, 13, 1)
    uv = poly_mul(u, v)
    assert poly_mul(u, u) == (Q - 1, 0, 0, 0)
    assert poly_mul(v, v) == poly_add(one, u)
    assert uv == (0, 2, 13, 0)
    basis = [one, u, v, uv]
    matrix = [[basis[column][row] for column in range(4)] for row in range(4)]
    assert determinant_mod_19(matrix) == 16
    tower = json.loads((HERE / "field_tower_circuits.json").read_text())
    assert set(tower) == {"F19_2", "F19_4", "model", "scalar"}
    assert set(tower["scalar"]) == {
        "addition",
        "multiplication",
        "multiplication_payload_sha256",
        "multiplication_stages",
        "negation",
        "subtraction",
        "validated_mul_pairs",
    }
    assert set(tower["F19_2"]) == {
        "addition",
        "multiplication",
        "representation",
        "subtraction",
    }
    assert set(tower["F19_4"]) == {
        "addition",
        "multiplication",
        "paper_basis",
        "subtraction",
        "tower",
        "tower_basis_determinant_mod_19",
        "u_t_coordinates",
        "uv_t_coordinates",
        "v_t_coordinates",
    }
    assert tower["scalar"]["multiplication_payload_sha256"] == net["payload_sha256"]
    assert tower["scalar"]["multiplication"] == 150
    assert tower["scalar"]["addition"] == 42
    assert tower["scalar"]["negation"] == 16
    assert tower["scalar"]["subtraction"] == 58
    assert tower["scalar"]["validated_mul_pairs"] == Q * Q
    assert tower["F19_2"]["multiplication"] == 5 * 42 + 2 * 16 + 3 * 150
    assert tower["F19_4"]["multiplication"] == (
        2 * 84 + 3 * 692 + (42 + 58) + 84 + (84 + 116)
    )
    assert tower["F19_4"]["tower_basis_determinant_mod_19"] == 16


def target_sampling_check() -> None:
    probabilities = []
    for _level, (_v, o, ell, public_columns) in ROWS:
        digits = o * ell * public_columns
        byte_length = math.ceil(8 * digits / 15)
        chunks, remainder_digits = divmod(digits, 15)
        bound = (2**64 // Q**15) * Q**15
        probability = Fraction(bound, 2**64) ** chunks
        if remainder_digits:
            remaining_bytes = byte_length - 8 * chunks
            bits = 8 * remaining_bytes
            partial_bound = (2**bits // Q**remainder_digits) * Q**remainder_digits
            probability *= Fraction(partial_bound, 2**bits)
        probabilities.append(probability)
    minimum = min(probabilities)
    assert abs(float(minimum) - 0.1524622985242944) < 1e-16


def structural_coefficient_bound_check() -> None:
    """Recompute the conditional fresh-coefficient union bound.

    This check assumes, but does not verify, the exact deterministic outer-map
    and alternating-source hypotheses stated in the paper.
    """
    failures = []
    for _level, (v, o, d, public_columns) in ROWS:
        if d != 4:
            continue
        m1 = math.ceil(o * public_columns / d)
        field_size = Q**d
        M = o * public_columns * d
        K = m1 * math.comb(d + 1, 2)
        delta_source = (
            Fraction(Q ** (M - K) - 1, Q - 1) * RHO ** (d * (v - 1))
        )
        projective_lines = Fraction(field_size ** (1 + o) - 1, field_size - 1) - 1
        delta_projection = (
            RHO ** (m1 * math.comb(d + 1, 2))
            + projective_lines * RHO ** (m1 * d * d)
        )
        failures.append(delta_source + delta_projection)
    assert max(failures) < Fraction(1, 2**209)


def dimension_and_repair_checks() -> None:
    combined_pairs, zero_offset_pairs = [], []
    for _level, (v, o, d, public_columns) in ROWS:
        x = o * public_columns
        M = d * x
        m1 = math.ceil(x / d)
        K = m1 * math.comb(d + 1, 2)
        assert K <= M
        if d == 4:
            assert d * (v - 1) - (M - K) >= K
        ordered_cap = min(M, m1 * d * d)
        symmetric_cap = math.comb(d + 1, 2)
        x_combined = max(
            v,
            math.ceil(ordered_cap / d),
            d * (math.ceil(ordered_cap / symmetric_cap) - 1) + 1,
        )
        combined_pairs.append((M, d * x_combined))
        m_square = d * v // symmetric_cap + 1
        x_zero = max(
            v + 1,
            d * (m_square - 1) + 1,
            math.ceil(ordered_cap / d),
            d * (math.ceil(ordered_cap / symmetric_cap) - 1) + 1,
        )
        zero_offset_pairs.append((M, d * x_zero))
    assert combined_pairs == [
        (80, 116),
        (64, 96),
        (80, 116),
        (112, 180),
        (96, 144),
        (100, 152),
        (144, 228),
        (128, 192),
        (144, 228),
    ]
    assert zero_offset_pairs == [
        (80, 180),
        (64, 130),
        (80, 180),
        (112, 260),
        (96, 194),
        (100, 244),
        (144, 324),
        (128, 258),
        (144, 324),
    ]
    combined = [Fraction(100 * (new - old), old) for old, new in combined_pairs]
    zero_offset = [
        Fraction(100 * (new - old), old) for old, new in zero_offset_pairs
    ]
    assert min(combined) == 45
    assert max(combined) == Fraction(425, 7)
    assert min(zero_offset) == Fraction(1625, 16)
    assert max(zero_offset) == 144


def ledger_checks() -> dict[str, dict[str, float]]:
    recomputed = {}
    for level, params in ROWS:
        key = str(params)
        if params[2] == 4:
            direct = l4_direct(level, params)
            fast = l4_fast(level, params)
            adaptive = l4_adaptive(direct, fast)
            recomputed[key] = {"direct": direct["direct"], "adaptive": adaptive}
        else:
            direct = l2_complete(level, params)
            diagonal = l2_diagonal(params)
            adaptive = l2_adaptive(level, direct, diagonal)
            recomputed[key] = {"direct": direct["direct"], "adaptive": adaptive}

    supplied = json.loads((HERE / "all_nine_ledger.json").read_text())
    for level, params in ROWS:
        key = str(params)
        row = supplied["rows"][key]
        assert "structural_probability" not in json.dumps(row)
        assert abs(
            float(recomputed[key]["direct"]) - row["direct_complete_square_exponent"]
        ) < 3e-9
        assert abs(
            float(recomputed[key]["adaptive"]) - row["adaptive"]["normalized_log2_AXN"]
        ) < 3e-9
    for route, field in (
        ("all_nine_direct_complete_square", "direct"),
        ("all_nine_adaptive", "adaptive"),
    ):
        for level in REF:
            maximum = max(
                float(recomputed[str(params)][field])
                for row_level, params in ROWS
                if row_level == level
            )
            assert abs(supplied[route][level]["exponent"] - maximum) < 3e-9
            assert abs(supplied[route][level]["headroom"] - (REF[level] - maximum)) < 3e-9
    return recomputed


def main() -> None:
    assert RHO == Fraction(14, 256)
    circuit_checks()
    target_sampling_check()
    structural_coefficient_bound_check()
    dimension_and_repair_checks()
    ledger_checks()
    print("Legacy ledger consistency checks passed")
    print("- direct/adaptive ledger rows recomputed without importing the generator")
    print("- no random structural-preflight probability appears in the ledger")
    print("- all 361 F19 multiplier pairs, payload SHA-256, field/tower identities, determinant")
    print("- target sampling, coefficient/spectrum bounds, dimensions, separators, repair floors")
    print("- stored ledger rows and separately recomputed levelwise maxima")
    print("- manuscript-table checks intentionally omitted because the legacy tables were removed")
    print("- excluded by design: verifier correspondence, reduction theorem, XOF model, H_hom, kappa_hom")


if __name__ == "__main__":
    main()
