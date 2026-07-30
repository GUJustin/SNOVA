#!/usr/bin/env python3
"""Exact arithmetic library for the SNOVA master-theorem upgrade.

The routines in this module provide the common structural probabilities,
seeded-spectrum tails, Bezout quantities, homotopy ledgers, and repair-floor
calculations used by the versioned wrapper scripts.  Version 4 uses the
zero-offset chosen-message construction for the missing Level-III orbit
slice; the older complete ordinary 70-by-70 routine is retained only as an
independent diagnostic fallback.

Exact integers and Fractions are used until displayed logarithms are taken.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

Q = 19
RHO = Fraction(14, 256)
G4_AXN = 6081

L2_ROWS = (
    ("I", (48, 16, 2, 2), 143),
    ("III", (72, 24, 2, 2), 207),
    ("V", (96, 32, 2, 2), 272),
)
L4_ROWS = (
    ("I-a", (28, 5, 4, 4), 143),
    ("I-b", (28, 4, 4, 5), 143),
    ("III-a", (40, 7, 4, 4), 207),
    ("III-b", (38, 5, 4, 5), 207),
    ("V-a", (50, 9, 4, 4), 272),
    ("V-b", (52, 6, 4, 6), 272),
)


def lg2(x: Fraction | int | float) -> float:
    if isinstance(x, Fraction):
        if x <= 0:
            raise ValueError("logarithm of nonpositive fraction")
        return math.log2(x.numerator) - math.log2(x.denominator)
    if x <= 0:
        raise ValueError("logarithm of nonpositive value")
    return math.log2(x)


def product(values: Iterable[Fraction]) -> Fraction:
    out = Fraction(1, 1)
    for x in values:
        out *= x
    return out


def p_atom_product(exponents: Iterable[int]) -> Fraction:
    return product(1 - RHO**e for e in exponents)


def ceil_extension_degree(base: int, threshold: int) -> int:
    r = 1
    while base**r < threshold:
        r += 1
    return r


def exact_l2_acceptance(n: int) -> Fraction:
    """Density with at most floor(n/4) affine-hyperplane block hits."""
    numerator = sum(
        math.comb(n, j) * (Q - 1) ** (n - j)
        for j in range(n // 4 + 1)
    )
    return Fraction(numerator, Q**n)


def l4_acceptance(n: int, r: int) -> Fraction:
    """Exact lower bound used for the public ell=4 rejection rule."""
    if r != 4:
        return Fraction(1, 1)
    return (1 - Fraction(1, Q**4)) ** n


def structural_preflight(d: int, v: int, o: int, m1: int, M: int) -> Fraction:
    """Theorem-12 random-XOF lower bound, retaining byte-to-field bias."""
    delta_src = Fraction(Q ** (m1 * math.comb(d, 2)) - 1, Q - 1) * RHO ** (d * (v - 1))
    A_size = Q**d
    other_lines = Fraction(A_size ** (1 + o) - 1, A_size - 1) - 1
    delta_proj = RHO ** (m1 * math.comb(d + 1, 2)) + other_lines * RHO ** (m1 * d * d)
    p_outer = p_atom_product(range(1, M + 1))
    return p_outer * (1 - delta_src) * (1 - delta_proj)


def spectrum_excess_tail(*, d: int, v: int, M: int, K: int, h: int, eta: Fraction) -> Fraction:
    """Theorem-27 excess-spectrum tail for a cross-data-selected slice."""
    t = d * (v - 2) - (M - K)
    if t < 0:
        raise ValueError("negative fresh-rank exponent")
    baseline = Fraction(Q**h - 1, Q**K)
    if eta <= baseline:
        raise ValueError("eta must exceed the full-rank baseline")
    return (
        Fraction(Q**h - 1, 1)
        * (1 - Fraction(1, Q**K))
        * RHO**t
        / (eta - baseline)
    )


def ordinary_square_root_probability(alpha: Fraction, K: int) -> tuple[Fraction, Fraction]:
    eta = Fraction(Q**K - 1, Q**K) + Fraction(1, Q)
    p = (alpha - eta / (Q - 1)) ** 2 / (1 + eta)
    return eta, p


def corrected_dense_quadratic_slp(nvars: int, neqs: int) -> int:
    """Conservative SLP for neqs affine quadratics in nvars variables."""
    quadratic_monomials = math.comb(nvars + 1, 2)
    terms_per_equation = quadratic_monomials + nvars + 1
    return quadratic_monomials + 2 * neqs * terms_per_equation


@dataclass
class L2Diagonal:
    level: str
    params: tuple[int, int, int, int]
    s: int
    K: int
    structural_probability: float
    jacobian_probability: float
    spectrum_failure_log2: float
    combined_key_probability: float
    root_probability_log2: float
    B_log2: float
    extension_degree_over_F19_4: int
    per_key_log2_AXN: float
    normalized_log2_AXN: float
    plus7_log2_AXN: float
    headroom: float
    output_MiB: float


def l2_diagonal(level: str, params: tuple[int, int, int, int], reference: int, s: int | None = None) -> L2Diagonal:
    v, o, d, rcols = params
    assert d == 2 and rcols == 2 and v == 3 * o
    m = o
    if s is None:
        s = m
    if not 1 <= s <= m:
        raise ValueError("invalid diagonal core size")
    n = v + o
    M = o * rcols * d
    K = 3 * m
    h = 2 * s
    alpha = exact_l2_acceptance(n)
    eta = Fraction(1, Q)
    delta = Fraction(s, Q**2)
    mu = Fraction(1, Q ** (K - h))
    root_constant = max(alpha - eta - delta, (alpha - delta) ** 2 / (1 + eta))
    p_root = mu * root_constant
    if p_root <= 0:
        raise ValueError("nonpositive root bound")

    B = 2**s
    Bplus = B + s * 2 ** (s - 1)
    slp = corrected_dense_quadratic_slp(s, s)
    r_ext = ceil_extension_degree(Q**4, max(2 * s, 8 * B * B))
    H = slp + 2 * s + s * s
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * s

    p_pre = structural_preflight(2, v, o, m, M)
    p_jac = p_atom_product(2 * c for c in range(1, s + 1))
    eps = spectrum_excess_tail(d=2, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (p_jac - eps)
    if p_key <= 0:
        raise ValueError("nonpositive combined key density")
    normalized = work / p_key

    output_bits = (s + 1) * B * (4 * r_ext) * math.log2(Q)
    return L2Diagonal(
        level=level,
        params=params,
        s=s,
        K=K,
        structural_probability=float(p_pre),
        jacobian_probability=float(p_jac),
        spectrum_failure_log2=lg2(eps),
        combined_key_probability=float(p_key),
        root_probability_log2=lg2(p_root),
        B_log2=lg2(B),
        extension_degree_over_F19_4=r_ext,
        per_key_log2_AXN=lg2(work),
        normalized_log2_AXN=lg2(normalized),
        plus7_log2_AXN=lg2(normalized) + 7,
        headroom=reference - lg2(normalized),
        output_MiB=output_bits / 8 / 2**20,
    )


@dataclass
class L2Complete:
    level: str
    params: tuple[int, int, int, int]
    s: int
    K: int
    structural_probability: float
    spectrum_failure_log2: float
    root_probability: float
    B_log2: float
    extension_degree_over_F19_4: int
    per_key_log2_AXN: float
    normalized_log2_AXN: float
    headroom: float


def l2_complete(level: str, params: tuple[int, int, int, int], reference: int) -> L2Complete:
    v, o, d, rcols = params
    assert d == 2 and rcols == 2 and o % 2 == 0
    m = o
    n = v + o
    M = 4 * m
    K = 3 * m
    s = 3 * m // 2
    h = K
    alpha = exact_l2_acceptance(n)
    eta, p_root = ordinary_square_root_probability(alpha, K)

    B = 2 ** (2 * m) * math.comb(m, m // 2)
    Bplus_f = Fraction(B, 1) * (1 + 2 * m + Fraction(m * m, m + 2))
    assert Bplus_f.denominator == 1
    Bplus = Bplus_f.numerator

    D = math.comb(s + 1, 2)
    diagonal_terms = D + s + 1
    cross_terms = (s + 1) ** 2
    slp = 2 * D + s * s + 4 * m * diagonal_terms + 2 * m * cross_terms
    r_ext = ceil_extension_degree(Q**4, max(2 * s, 8 * B * B))
    H = slp + 6 * m + K * K
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * K

    p_pre = structural_preflight(2, v, o, m, M)
    eps = spectrum_excess_tail(d=2, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - eps)
    normalized = work / p_key
    return L2Complete(
        level=level,
        params=params,
        s=s,
        K=K,
        structural_probability=float(p_pre),
        spectrum_failure_log2=lg2(eps),
        root_probability=float(p_root),
        B_log2=lg2(B),
        extension_degree_over_F19_4=r_ext,
        per_key_log2_AXN=lg2(work),
        normalized_log2_AXN=lg2(normalized),
        headroom=reference - lg2(normalized),
    )


@dataclass
class L2Ordinary:
    level: str
    params: tuple[int, int, int, int]
    K: int
    D0_lower: int
    spectrum_failure_log2: float
    root_probability: float
    normalized_log2_AXN: float
    headroom: float


def l2_ordinary(level: str, params: tuple[int, int, int, int], reference: int) -> L2Ordinary:
    v, o, d, rcols = params
    m = o
    M = 4 * m
    K = 3 * m
    D0 = 2 * (v - 1) - (M - K)
    if D0 < K:
        raise ValueError("ordinary square not guaranteed")
    alpha = exact_l2_acceptance(v + o)
    eta, p_root = ordinary_square_root_probability(alpha, K)
    B = 2**K
    Bplus = B + K * 2 ** (K - 1)
    slp = corrected_dense_quadratic_slp(K, K)
    r_ext = ceil_extension_degree(Q**4, max(2 * K, 8 * B * B))
    H = slp + 2 * K + K * K
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * K
    p_pre = structural_preflight(2, v, o, m, M)
    eps = spectrum_excess_tail(d=2, v=v, M=M, K=K, h=K, eta=eta)
    normalized = work / (p_pre * (1 - eps))
    return L2Ordinary(
        level=level,
        params=params,
        K=K,
        D0_lower=D0,
        spectrum_failure_log2=lg2(eps),
        root_probability=float(p_root),
        normalized_log2_AXN=lg2(normalized),
        headroom=reference - lg2(normalized),
    )


@dataclass
class L4Channel:
    name: str
    params: tuple[int, int, int, int]
    profile: tuple[int, int, int]
    structural_probability: float
    jacobian_failure_log2: float
    spectrum_failure_log2: float
    combined_key_probability: float
    per_key_log2_AXN: float
    normalized_log2_AXN: float
    plus7_log2_AXN: float
    plus7_headroom: float


L4_CHANNEL_PROFILES = {
    "I": (9, 5, 4),
    "III": (10, 7, 3),
    "V": (11, 9, 2),
}


def l4_channel(name: str, params: tuple[int, int, int, int], reference: int) -> L4Channel:
    v, o, d, rcols = params
    assert d == 4
    M = o * rcols * d
    m1 = (o * rcols + d - 1) // d
    K = 10 * m1
    level = name.split("-")[0]
    s, a, b = L4_CHANNEL_PROFILES[level]
    assert a + b == s and a <= m1 and b <= m1
    h = 4 * s
    alpha = l4_acceptance(v + o, rcols)
    eta = Fraction(1, 2)
    delta = Fraction(a + Q * b, Q**4)
    mu = Fraction(1, Q ** (K - h))
    p_root = mu * (alpha - delta) ** 2 / (1 + eta)

    B = 2**a * (Q + 1) ** b
    Bplus_f = Fraction(B, 1) * (1 + Fraction(a, 2) + Fraction(b, Q + 1))
    assert Bplus_f.denominator == 1
    Bplus = Bplus_f.numerator

    D = math.comb(s + 1, 2)
    slp0 = 0 if a == 0 else D + 2 * a * (D + s + 1)
    slp1 = 0 if b == 0 else 6 * s + s * s + 2 * b * (s * s + 2 * s + 1)
    slp = slp0 + slp1
    degree_sum = 2 * a + (Q + 1) * b
    r_ext = ceil_extension_degree(Q**4, max(degree_sum, 8 * B * B))
    H = slp + degree_sum + s * s
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * s

    p_pre = structural_preflight(4, v, o, m1, M)
    # Rank-deficient s-by-s A-Jacobian: union over projective left-kernel lines.
    A_size = Q**4
    projective_lines = Fraction(A_size**s - 1, A_size - 1)
    jac_fail = projective_lines * RHO ** (4 * s)
    eps = spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - jac_fail - eps)
    if p_key <= 0:
        raise ValueError("nonpositive l4 key-density lower bound")
    normalized = work / p_key
    return L4Channel(
        name=name,
        params=params,
        profile=(s, a, b),
        structural_probability=float(p_pre),
        jacobian_failure_log2=lg2(jac_fail),
        spectrum_failure_log2=lg2(eps),
        combined_key_probability=float(p_key),
        per_key_log2_AXN=lg2(work),
        normalized_log2_AXN=lg2(normalized),
        plus7_log2_AXN=lg2(normalized) + 7,
        plus7_headroom=reference - lg2(normalized) - 7,
    )


@dataclass
class L4OrdinaryLevelIII:
    params: tuple[int, int, int, int]
    D0_lower: int
    structural_probability: float
    spectrum_failure_log2: float
    per_key_log2_AXN: float
    normalized_log2_AXN: float
    headroom: float


def l4_level3_ordinary(params: tuple[int, int, int, int], reference: int = 207) -> L4OrdinaryLevelIII:
    v, o, d, rcols = params
    M = o * rcols * d
    m1 = (o * rcols + d - 1) // d
    K = 10 * m1
    assert K == 70
    D0 = d * (v - 1) - (M - K)
    if D0 < K:
        raise ValueError("complete ordinary Level-III square unavailable")
    alpha = l4_acceptance(v + o, rcols)
    eta, p_root = ordinary_square_root_probability(alpha, K)
    B = 2**K
    Bplus = B + K * 2 ** (K - 1)
    slp = corrected_dense_quadratic_slp(K, K)
    # The manuscript's finite-cardinality adaptation originally retained an
    # extra (K-1) in this threshold. We keep it to reproduce its conservative
    # complete-square ledger.
    r_ext = ceil_extension_degree(Q**4, max(2 * K, 8 * (K - 1) * B * B))
    H = slp + 2 * K + K * K
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * K
    p_pre = structural_preflight(4, v, o, m1, M)
    eps = spectrum_excess_tail(d=4, v=v, M=M, K=K, h=K, eta=eta)
    normalized = work / (p_pre * (1 - eps))
    return L4OrdinaryLevelIII(
        params=params,
        D0_lower=D0,
        structural_probability=float(p_pre),
        spectrum_failure_log2=lg2(eps),
        per_key_log2_AXN=lg2(work),
        normalized_log2_AXN=lg2(normalized),
        headroom=reference - lg2(normalized),
    )



def l4_channel_profile(name: str, params: tuple[int, int, int, int], reference: int,
                       s: int, a: int, b: int) -> L4Channel:
    """Evaluate an arbitrary whole-channel core (a degree-2, b degree-20)."""
    v, o, d, rcols = params
    assert d == 4 and a + b == s
    M = o * rcols * d
    m1 = (o * rcols + d - 1) // d
    K = 10 * m1
    Hdim = v + o - 4 * m1
    if not (1 <= s <= Hdim and 0 <= a <= m1 and 0 <= b <= m1 and 4 * s <= K):
        raise ValueError("inadmissible whole-channel profile")
    h = 4 * s
    alpha = l4_acceptance(v + o, rcols)
    eta = Fraction(1, 2)
    delta = Fraction(a + Q * b, Q**4)
    if alpha <= delta:
        raise ValueError("singular-density charge exhausts accepted density")
    mu = Fraction(1, Q ** (K - h))
    p_root = mu * (alpha - delta) ** 2 / (1 + eta)
    B = 2**a * (Q + 1) ** b
    Bplus_f = Fraction(B, 1) * (1 + Fraction(a, 2) + Fraction(b, Q + 1))
    if Bplus_f.denominator != 1:
        raise ValueError("nonintegral companion Bezout quantity")
    Bplus = Bplus_f.numerator
    D = math.comb(s + 1, 2)
    slp0 = 0 if a == 0 else D + 2 * a * (D + s + 1)
    slp1 = 0 if b == 0 else 6 * s + s * s + 2 * b * (s * s + 2 * s + 1)
    degree_sum = 2 * a + (Q + 1) * b
    r_ext = ceil_extension_degree(Q**4, max(degree_sum, 8 * B * B))
    H = slp0 + slp1 + degree_sum + s * s
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * G4_AXN * B * Bplus * H * s
    p_pre = structural_preflight(4, v, o, m1, M)
    A_size = Q**4
    projective_lines = Fraction(A_size**s - 1, A_size - 1)
    jac_fail = projective_lines * RHO ** (4 * s)
    eps = spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - jac_fail - eps)
    if p_key <= 0:
        raise ValueError("nonpositive l4 key density")
    normalized = work / p_key
    return L4Channel(
        name=name,
        params=params,
        profile=(s, a, b),
        structural_probability=float(p_pre),
        jacobian_failure_log2=lg2(jac_fail),
        spectrum_failure_log2=lg2(eps),
        combined_key_probability=float(p_key),
        per_key_log2_AXN=lg2(work),
        normalized_log2_AXN=lg2(normalized),
        plus7_log2_AXN=lg2(normalized) + 7,
        plus7_headroom=reference - lg2(normalized) - 7,
    )


def l4_channel_frontier(name: str, params: tuple[int, int, int, int], reference: int) -> dict[str, Any]:
    v, o, d, rcols = params
    m1 = (o * rcols + d - 1) // d
    Hdim = v + o - d * m1
    rows: list[tuple[L4Channel, float]] = []
    for a in range(m1 + 1):
        for b in range(m1 + 1):
            s = a + b
            if s == 0 or s > Hdim or 4 * s > 10 * m1:
                continue
            row = l4_channel_profile(name, params, reference, s, a, b)
            B = 2**a * (Q + 1) ** b
            r_ext = ceil_extension_degree(Q**4, max(2 * a + (Q + 1) * b, 8 * B * B))
            output_bits = (s + 1) * B * (4 * r_ext) * math.log2(Q)
            rows.append((row, output_bits / 8 / 2**30))
    fast, fast_output = min(rows, key=lambda z: z[0].normalized_log2_AXN)
    eligible = [z for z in rows if z[0].plus7_log2_AXN <= reference]
    low, low_output = min(eligible, key=lambda z: z[1])
    return {
        "fast": {**asdict(fast), "output_GiB": fast_output},
        "low_output_plus7": {**asdict(low), "output_GiB": low_output},
        "profiles_tested": len(rows),
    }


def fresh_square_floor(params: tuple[int, int, int, int]) -> dict[str, Any]:
    v, o, d, rcols = params
    M = o * rcols * d
    minimum = d * v
    return {
        "current_M": M,
        "fresh_square_floor_M": minimum,
        "fresh_square_increase_percent": 100 * (minimum / M - 1),
        "guarantee_persists_while": f"M' <= {d * (v - 1)}",
    }


def ordered_cap_floor(params: tuple[int, int, int, int]) -> int:
    """Proposition-49 floor under m1=ceil(or/d), for official Kord=M."""
    v, o, d, rcols = params
    M = o * rcols * d
    mprime = math.ceil(M / math.comb(d + 1, 2))
    xprime = max(math.ceil(M / d), d * (mprime - 1) + 1)
    return d * xprime


def combined_repair_floor(params: tuple[int, int, int, int]) -> dict[str, Any]:
    v, o, d, rcols = params
    M = o * rcols * d
    ordered = ordered_cap_floor(params)
    square = d * v
    combined = max(ordered, square)
    return {
        "current_M": M,
        "ordered_cap_floor_M": ordered,
        "fresh_square_floor_M": square,
        "combined_floor_M": combined,
        "combined_increase_percent": 100 * (combined / M - 1),
    }


def main() -> None:
    report: dict[str, Any] = {
        "constants": {"q": Q, "rho": [RHO.numerator, RHO.denominator], "g4_AXN": G4_AXN},
        "l2": {},
        "l4_channel": {},
        "l4_channel_frontier": {},
        "l4_level3_ordinary": {},
        "all_nine_fixed_core_free": {},
        "repair_floors": {},
        "current_pdf_regression": {},
    }

    for level, params, ref in L2_ROWS:
        diag = l2_diagonal(level, params, ref)
        complete = l2_complete(level, params, ref)
        ordinary = l2_ordinary(level, params, ref)
        report["l2"][level] = {
            "diagonal": asdict(diag),
            "complete_fixed_core_free": asdict(complete),
            "ordinary_square": asdict(ordinary),
        }
        report["repair_floors"][str(params)] = combined_repair_floor(params)

    for name, params, ref in L4_ROWS:
        ch = l4_channel(name, params, ref)
        report["l4_channel"][name] = asdict(ch)
        report["l4_channel_frontier"][name] = l4_channel_frontier(name, params, ref)
        report["repair_floors"][str(params)] = combined_repair_floor(params)
        if name.startswith("III"):
            report["l4_level3_ordinary"][name] = asdict(l4_level3_ordinary(params, ref))

    # Correct fixed-core-free frontier. Level-I/V use the current manuscript's
    # orbit-complete normalized values; Level III must use the complete ordinary
    # square for the shape whose guaranteed A-kernel dimension is only 15.
    fixed_l4 = {
        "I": {"normalized": 138.40973, "route": "16-orbit eigenblock sweep, h=48"},
        "III": {
            "normalized": max(x["normalized_log2_AXN"] for x in report["l4_level3_ordinary"].values()),
            "route": "complete ordinary 70-by-70 square system",
        },
        "V": {"normalized": 220.91810, "route": "16-orbit eigenblock sweep, h=88"},
    }
    for level, _, ref in L2_ROWS:
        l2_norm = report["l2"][level]["complete_fixed_core_free"]["normalized_log2_AXN"]
        l4_norm = fixed_l4[level]["normalized"]
        worst = max(l2_norm, l4_norm)
        report["all_nine_fixed_core_free"][level] = {
            "l2_normalized": l2_norm,
            "l4_normalized": l4_norm,
            "worst_normalized": worst,
            "headroom": ref - worst,
            "l4_route": fixed_l4[level]["route"],
        }

    # Explicitly expose the current-PDF dimension regression.
    report["current_pdf_regression"] = {
        "claimed_level_III_A_dimension": 17,
        "required_h": 68,
        "shape": (38, 5, 4, 5),
        "guaranteed_stable_A_dimension": 15,
        "maximum_guaranteed_h": 60,
        "valid_replacement": "complete ordinary K=70 square system",
    }

    out = Path(__file__).resolve().with_name("snova_master_upgrade_numbers.json")
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print("SNOVA master-theorem upgrade ledger")
    print("=" * 88)
    for level in ("I", "III", "V"):
        l2 = report["l2"][level]
        print(f"Level {level} ell=2:")
        print(
            "  one-eigenblock normalized "
            f"2^{l2['diagonal']['normalized_log2_AXN']:.6f}, +7 -> "
            f"2^{l2['diagonal']['plus7_log2_AXN']:.6f}"
        )
        print(
            "  fixed-core-free complete normalized "
            f"2^{l2['complete_fixed_core_free']['normalized_log2_AXN']:.6f}"
        )
    print()
    for name in report["l4_channel_frontier"]:
        fast = report["l4_channel_frontier"][name]["fast"]
        low = report["l4_channel_frontier"][name]["low_output_plus7"]
        print(
            f"{name} ell=4 channel fast: 2^{fast['normalized_log2_AXN']:.6f} "
            f"profile {tuple(fast['profile'])}, output {fast['output_GiB']:.3g} GiB"
        )
        print(
            f"  low-output +7: 2^{low['normalized_log2_AXN']:.6f}, "
            f"+7 -> 2^{low['plus7_log2_AXN']:.6f}, output {low['output_GiB']:.3g} GiB"
        )
    print()
    for level, x in report["all_nine_fixed_core_free"].items():
        print(
            f"All-nine fixed-core-free Level {level}: worst 2^{x['worst_normalized']:.6f}, "
            f"headroom {x['headroom']:.6f}"
        )
    print()
    floors = [x["combined_increase_percent"] for x in report["repair_floors"].values()]
    print(f"Combined repair-floor increases range from {min(floors):.2f}% to {max(floors):.2f}%")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
