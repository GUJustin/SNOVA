#!/usr/bin/env python3
"""Version-5 numerical ledger for the SNOVA symmetry-quotient paper.

New over v4:
  * direct F19^2 homotopy arithmetic for every ell=2 A-linear core;
  * a constructive F19^4 tower arithmetic bound for ell=4 cores;
  * exact separator-success accounting, optimizing
        (2r-1)/(1-B^2/|A|^r)
    instead of forcing |A|^r >= 8B^2 and then charging 8/7;
  * the corrected one-denominator zero-offset spectrum tail;
  * all tables regenerated from exact integers/Fractions.

The 16-orbit combinatorial ledger remains the v4 exhaustive artifact input.
It is re-costed only by the uniform F19^4 arithmetic improvement; no unpriced
per-orbit separator improvement is claimed.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


base = load_module("snova_base_v5", HERE / "snova_master_upgrade_numbers.py")
v4 = load_module("snova_v4_for_v5", HERE / "snova_master_upgrade_v4.py")

Q = 19
RHO = Fraction(14, 256)
G1_AXN = 441
_CIRCUITS = json.loads((HERE / "snova_field_tower_circuits_v5.json").read_text())
G2_AXN = int(_CIRCUITS["F19_2"]["multiplication"])
G4_AXN = int(_CIRCUITS["F19_4"]["multiplication"])
A2 = Q**2
A4 = Q**4
SPECTRUM_EXCESS_BITS = 32
SPECTRUM_EXCESS = Fraction(1, 2**SPECTRUM_EXCESS_BITS)

# v4 values include exact orbit enumeration, eta conversion, key density, and
# random within-family projection.  They used G4=6081.  The arithmetic schedule
# is linear in that primitive, so the re-costing shift is exact.
ORBIT_V4_COMPACT = {
    "I": v4.compact_orbit("I")["compact_corrected_normalized_log2_AXN"],
    "III": v4.compact_orbit("III")["compact_corrected_normalized_log2_AXN"],
    "V": v4.compact_orbit("V")["compact_corrected_normalized_log2_AXN"],
}
ORBIT_ARITH_SHIFT = math.log2(G4_AXN / base.G4_AXN)


def lg2(x: Fraction | int | float) -> float:
    return base.lg2(x)


def spectrum_threshold(h: int, K: int) -> Fraction:
    """Full-rank baseline plus an explicit 2^-32 excess allowance."""
    return Fraction(Q**h - 1, Q**K) + SPECTRUM_EXCESS


def accepted_nonsingular_probability(*, alpha: Fraction, h: int, K: int, eta: Fraction) -> Fraction:
    """Theorem-21 accepted full-Jacobian-root probability."""
    return Fraction(1, Q ** (K - h)) * (alpha - eta / (Q - 1)) ** 2 / (1 + eta)


def exact_homotopy_factor(
    *, base_size: int, B: int, degree_bound: int, gate_cost: int
) -> dict[str, Any]:
    """Optimize extension degree using exact separator success.

    The finite-cardinality proof needs enough base-extension elements for the
    deterministic start nodes and a random separating form.  At degree r, the
    latter fails with probability at most B^2/base_size^r.  Evaluation and
    interpolation use 2r-1 base-field products.
    """
    if B <= 0 or degree_bound <= 0:
        raise ValueError("invalid homotopy parameters")
    r0 = 1
    while base_size**r0 < degree_bound or base_size**r0 <= B * B:
        r0 += 1
    best: tuple[Fraction, int, Fraction] | None = None
    # The factor eventually grows linearly.  Sixty-four steps beyond the first
    # admissible r is ample and is independently checked in v5_checks.py.
    for r in range(r0, r0 + 65):
        field_size = base_size**r
        if 2 * r - 1 > base_size:
            raise ValueError("not enough base-field evaluation points")
        success = 1 - Fraction(B * B, field_size)
        if success <= 0:
            continue
        factor = Fraction((2 * r - 1) * gate_cost, 1) / success
        if best is None or factor < best[0]:
            best = (factor, r, success)
    if best is None:
        raise RuntimeError("no admissible extension")
    factor, r, success = best
    old_r = base.ceil_extension_degree(base_size, max(degree_bound, 8 * B * B))
    old_factor_same_gate = Fraction(8, 7) * (2 * old_r - 1) * gate_cost
    return {
        "extension_degree": r,
        "field_size": base_size**r,
        "separator_success": [success.numerator, success.denominator],
        "separator_success_float": float(success),
        "pointwise_products": 2 * r - 1,
        "gate_factor": [factor.numerator, factor.denominator],
        "gate_factor_log2": lg2(factor),
        "old_8_over_7_extension_degree": old_r,
        "gain_vs_8_over_7_same_gate_bits": lg2(old_factor_same_gate / factor),
    }


def factor_fraction(info: dict[str, Any]) -> Fraction:
    n, d = info["gate_factor"]
    return Fraction(n, d)


def l2_diagonal_v5(level: str, params: tuple[int, int, int, int], reference: int, s: int | None = None) -> dict[str, Any]:
    v, o, d, rcols = params
    if d != 2 or rcols != 2 or v != 3 * o:
        raise ValueError("unexpected ell=2 row")
    m = o
    if s is None:
        s = m
    if not 1 <= s <= m:
        raise ValueError("invalid diagonal core size")
    n = v + o
    M = 4 * m
    K = 3 * m
    h = 2 * s
    alpha = base.exact_l2_acceptance(n)
    eta = spectrum_threshold(h, K)
    delta = Fraction(s, A2)
    mu = Fraction(1, Q ** (K - h))
    root_constant = max(alpha - eta - delta, (alpha - delta) ** 2 / (1 + eta))
    p_root = mu * root_constant
    if p_root <= 0:
        raise ValueError("nonpositive root probability")

    B = 2**s
    Bplus = B + s * 2 ** (s - 1)
    slp = base.corrected_dense_quadratic_slp(s, s)
    H = slp + 2 * s + s * s
    hom = exact_homotopy_factor(base_size=A2, B=B, degree_bound=2 * s, gate_cost=G2_AXN)
    work = Fraction(1, 1) / p_root * factor_fraction(hom) * B * Bplus * H * s

    p_pre = base.structural_preflight(2, v, o, m, M)
    p_jac = base.p_atom_product(2 * c for c in range(1, s + 1))
    eps = base.spectrum_excess_tail(d=2, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (p_jac - eps)
    if p_key <= 0:
        raise ValueError("nonpositive key density")
    normalized = work / p_key
    r_ext = hom["extension_degree"]
    output_bits = (s + 1) * B * (2 * r_ext) * math.log2(Q)
    old = base.l2_diagonal(level, params, reference, s)
    return {
        "level": level,
        "parameters": list(params),
        "route": "one-eigenblock diagonal core over F19^2",
        "s": s,
        "h": h,
        "K": K,
        "eta": [eta.numerator, eta.denominator],
        "structural_probability": float(p_pre),
        "jacobian_probability": float(p_jac),
        "spectrum_failure_log2": lg2(eps),
        "combined_key_probability": float(p_key),
        "root_probability_log2": lg2(p_root),
        "B_log2": lg2(B),
        "Bplus_log2": lg2(Bplus),
        "homotopy": hom,
        "per_key_log2_AXN": lg2(work),
        "normalized_log2_AXN": lg2(normalized),
        "headroom": reference - lg2(normalized),
        "plus7_headroom": reference - lg2(normalized) - 7,
        "output_MiB": output_bits / 8 / 2**20,
        "v4_normalized_log2_AXN": old.normalized_log2_AXN,
        "v5_gain_bits": old.normalized_log2_AXN - lg2(normalized),
    }


def l2_complete_v5(level: str, params: tuple[int, int, int, int], reference: int) -> dict[str, Any]:
    v, o, d, rcols = params
    if d != 2 or rcols != 2 or o % 2:
        raise ValueError("unexpected ell=2 row")
    m = o
    n = v + o
    M = 4 * m
    K = 3 * m
    s = K // 2
    h = K
    alpha = base.exact_l2_acceptance(n)
    eta = spectrum_threshold(h, K)
    p_root = accepted_nonsingular_probability(alpha=alpha, h=h, K=K, eta=eta)

    B = 2 ** (2 * m) * math.comb(m, m // 2)
    Bplus_f = Fraction(B, 1) * (1 + 2 * m + Fraction(m * m, m + 2))
    if Bplus_f.denominator != 1:
        raise ValueError("nonintegral companion Bezout number")
    Bplus = Bplus_f.numerator
    D = math.comb(s + 1, 2)
    diagonal_terms = D + s + 1
    cross_terms = (s + 1) ** 2
    slp = 2 * D + s * s + 4 * m * diagonal_terms + 2 * m * cross_terms
    H = slp + 6 * m + K * K
    hom = exact_homotopy_factor(base_size=A2, B=B, degree_bound=2 * s, gate_cost=G2_AXN)
    work = Fraction(1, 1) / p_root * factor_fraction(hom) * B * Bplus * H * K

    p_pre = base.structural_preflight(2, v, o, m, M)
    eps = base.spectrum_excess_tail(d=2, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - eps)
    normalized = work / p_key
    r_ext = hom["extension_degree"]
    output_bits = (K + 1) * B * (2 * r_ext) * math.log2(Q)
    old = base.l2_complete(level, params, reference)
    return {
        "level": level,
        "parameters": list(params),
        "route": "complete two-block square system over F19^2",
        "s": s,
        "h": h,
        "K": K,
        "eta": [eta.numerator, eta.denominator],
        "structural_probability": float(p_pre),
        "spectrum_failure_log2": lg2(eps),
        "root_probability": float(p_root),
        "B_log2": lg2(B),
        "Bplus_log2": lg2(Bplus),
        "homotopy": hom,
        "per_key_log2_AXN": lg2(work),
        "normalized_log2_AXN": lg2(normalized),
        "headroom": reference - lg2(normalized),
        "plus7_headroom": reference - lg2(normalized) - 7,
        "output_GiB": output_bits / 8 / 2**30,
        "v4_normalized_log2_AXN": old.normalized_log2_AXN,
        "v5_gain_bits": old.normalized_log2_AXN - lg2(normalized),
    }


def l4_channel_profile_v5(
    name: str,
    params: tuple[int, int, int, int],
    reference: int,
    s: int,
    a: int,
    b: int,
    eta: Fraction,
) -> dict[str, Any]:
    v, o, d, rcols = params
    if d != 4 or a + b != s:
        raise ValueError("invalid whole-channel profile")
    M = o * rcols * d
    m1 = (o * rcols + d - 1) // d
    K = 10 * m1
    Hdim = v + o - d * m1
    if not (1 <= s <= Hdim and 0 <= a <= m1 and 0 <= b <= m1 and 4 * s <= K):
        raise ValueError("inadmissible whole-channel profile")
    h = 4 * s
    alpha = base.l4_acceptance(v + o, rcols)
    delta = Fraction(a + Q * b, A4)
    if alpha <= delta:
        raise ValueError("singular-density loss too large")
    mu = Fraction(1, Q ** (K - h))
    p_root = mu * (alpha - delta) ** 2 / (1 + eta)

    B = 2**a * (Q + 1) ** b
    Bplus_f = Fraction(B, 1) * (1 + Fraction(a, 2) + Fraction(b, Q + 1))
    if Bplus_f.denominator != 1:
        raise ValueError("nonintegral companion Bezout number")
    Bplus = Bplus_f.numerator
    D = math.comb(s + 1, 2)
    slp0 = 0 if a == 0 else D + 2 * a * (D + s + 1)
    slp1 = 0 if b == 0 else 6 * s + s * s + 2 * b * (s * s + 2 * s + 1)
    degree_bound = 2 * a + (Q + 1) * b
    H = slp0 + slp1 + degree_bound + s * s
    hom = exact_homotopy_factor(base_size=A4, B=B, degree_bound=degree_bound, gate_cost=G4_AXN)
    work = Fraction(1, 1) / p_root * factor_fraction(hom) * B * Bplus * H * s

    p_pre = base.structural_preflight(4, v, o, m1, M)
    projective_lines = Fraction(A4**s - 1, A4 - 1)
    jac_fail = projective_lines * RHO ** (4 * s)
    eps = base.spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - jac_fail - eps)
    if p_key <= 0:
        raise ValueError("nonpositive key density")
    normalized = work / p_key
    r_ext = hom["extension_degree"]
    output_bits = (s + 1) * B * (4 * r_ext) * math.log2(Q)
    old = v4.l4_channel_profile_eta(name, params, reference, s, a, b, eta)
    return {
        "name": name,
        "parameters": list(params),
        "profile": [s, a, b],
        "eta": [eta.numerator, eta.denominator],
        "structural_probability": float(p_pre),
        "jacobian_failure_log2": lg2(jac_fail),
        "spectrum_failure_log2": lg2(eps),
        "combined_key_probability": float(p_key),
        "root_probability_log2": lg2(p_root),
        "B_log2": lg2(B),
        "Bplus_log2": lg2(Bplus),
        "homotopy": hom,
        "per_key_log2_AXN": lg2(work),
        "normalized_log2_AXN": lg2(normalized),
        "plus7_log2_AXN": lg2(normalized) + 7,
        "plus7_headroom": reference - lg2(normalized) - 7,
        "output_GiB": output_bits / 8 / 2**30,
        "v4_same_profile_normalized_log2_AXN": old["normalized_log2_AXN"],
        "v5_same_profile_gain_bits": old["normalized_log2_AXN"] - lg2(normalized),
    }


def l4_channel_frontier_v5(
    name: str, params: tuple[int, int, int, int], reference: int
) -> dict[str, Any]:
    v, o, d, rcols = params
    m1 = (o * rcols + d - 1) // d
    Hdim = v + o - d * m1
    rows: list[dict[str, Any]] = []
    for a in range(m1 + 1):
        for b in range(m1 + 1):
            s = a + b
            if s == 0 or s > Hdim or 4 * s > 10 * m1:
                continue
            K = 10 * m1
            eta = spectrum_threshold(4 * s, K)
            rows.append(l4_channel_profile_v5(name, params, reference, s, a, b, eta))
    fast = min(rows, key=lambda z: z["normalized_log2_AXN"])
    eligible = [z for z in rows if z["plus7_log2_AXN"] <= reference]
    low = min(eligible, key=lambda z: z["output_GiB"])
    return {"fast": fast, "low_output_plus7": low, "profiles_tested": len(rows)}


def orbit_v5(level: str) -> dict[str, Any]:
    """Re-cost the exact 16-orbit schedule at a near-baseline spectrum threshold."""
    reference = {"I": 143, "III": 207, "V": 272}[level]
    h = {"I": 48, "III": 68, "V": 88}[level]
    K = {"I": 50, "III": 70, "V": 90}[level]
    v, M = {"I": (28, 80), "III": (40, 112), "V": (50, 144)}[level]
    alpha = v4._orbit_acceptance(level)
    eta_old = Fraction(1, Q)
    eta_new = spectrum_threshold(h, K)
    p_old = accepted_nonsingular_probability(alpha=alpha, h=h, K=K, eta=eta_old)
    p_new = accepted_nonsingular_probability(alpha=alpha, h=h, K=K, eta=eta_new)
    eps_new = base.spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta_new)
    old = ORBIT_V4_COMPACT[level]
    root_shift = lg2(p_old / p_new)
    spectrum_charge = -lg2(1 - eps_new)
    new = old + ORBIT_ARITH_SHIFT + root_shift + spectrum_charge
    return {
        "level": level,
        "h": h,
        "K": K,
        "eta": [eta_new.numerator, eta_new.denominator],
        "spectrum_failure_log2": lg2(eps_new),
        "v4_compact_projection_complete_log2_AXN": old,
        "arithmetic_shift_bits": ORBIT_ARITH_SHIFT,
        "root_probability_shift_bits": root_shift,
        "spectrum_charge_bits": spectrum_charge,
        "v5_projection_complete_log2_AXN": new,
        "headroom": reference - new,
        "scope": (
            f"exact re-costing of the v4 16-orbit schedule by {G4_AXN}/{base.G4_AXN}; "
            "the per-orbit 8/7 separator reserve is retained"
        ),
    }


def zero_offset_level3b_orbit_v5() -> dict[str, Any]:
    """Re-cost the repaired (38,5,4,5) zero-offset orbit route sharply."""
    old = v4.zero_offset_level3b_orbit()
    h, K, v = 68, 70, 38
    alpha = v4._orbit_acceptance("III")
    eta_old = Fraction(1, Q)
    eta_new = spectrum_threshold(h, K)
    p_old = accepted_nonsingular_probability(alpha=alpha, h=h, K=K, eta=eta_old)
    p_new = accepted_nonsingular_probability(alpha=alpha, h=h, K=K, eta=eta_new)
    eps_old = v4.zero_offset_spectrum_tail(d=4, v=v, K=K, h=h, eta=eta_old)
    eps_new = v4.zero_offset_spectrum_tail(d=4, v=v, K=K, h=h, eta=eta_new)
    root_shift = lg2(p_old / p_new)
    spectrum_shift = lg2((1 - eps_old) / (1 - eps_new))
    solve = old["solve_normalized_log2_AXN"] + ORBIT_ARITH_SHIFT + root_shift + spectrum_shift
    target = old["target_generation_envelope_log2_AXN"]
    combined = max(solve, target) + math.log2(1 + 2 ** (-abs(solve - target)))
    return {
        **old,
        "eta": [eta_new.numerator, eta_new.denominator],
        "v4_solve_normalized_log2_AXN": old["solve_normalized_log2_AXN"],
        "arithmetic_shift_bits": ORBIT_ARITH_SHIFT,
        "root_probability_shift_bits": root_shift,
        "spectrum_shift_bits": spectrum_shift,
        "spectrum_failure_log2": lg2(eps_new),
        "solve_normalized_log2_AXN": solve,
        "combined_log2_AXN": combined,
        "headroom": 207 - combined,
    }


def random_cut_row_v5(name: str, params: tuple[int, int, int, int], reference: int) -> dict[str, Any]:
    v, o, d, rcols = params
    if d != 4 or rcols != 4:
        raise ValueError("random-cut row requires square ell=4 shape")
    m = (o * rcols + d - 1) // d
    M, K = o * rcols * d, 10 * m
    h, s = K + 2, (K + 2) // 4
    if 4 * s != h or s > v + o - 4 * m:
        raise ValueError("required stable slice unavailable")
    B, Bplus = v4.random_cut_bezout(m, s)
    alpha = base.l4_acceptance(v + o, rcols)
    theta = spectrum_threshold(h, h)
    p_root = (alpha - theta / (Q - 1)) ** 2 / (1 + theta)
    D = math.comb(s + 1, 2)
    slp = 4 * D + 6 * s * s + 2 * (4 * m * (D + s + 1) + 6 * m * (s + 1) ** 2) + 4 * (4 * s + 1)
    degree_bound = 5 * m + 2
    H = slp + (20 * m + 8) + h * h
    hom = exact_homotopy_factor(base_size=A4, B=B, degree_bound=degree_bound, gate_cost=G4_AXN)
    work = Fraction(1, 1) / p_root * factor_fraction(hom) * B * Bplus * H * h
    p_pre = base.structural_preflight(4, v, o, m, M)
    eta_bar = Q**2 * theta
    eps = base.spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta_bar)
    normalized = work / (p_pre * (1 - eps))
    old = v4.random_cut_row(name, params, reference)
    return {
        "name": name,
        "parameters": list(params),
        "h": h,
        "K": K,
        "s": s,
        "B_log2": lg2(B),
        "Bplus_log2": lg2(Bplus),
        "root_probability": float(p_root),
        "spectrum_failure_log2": lg2(eps),
        "homotopy": hom,
        "normalized_log2_AXN": lg2(normalized),
        "headroom": reference - lg2(normalized),
        "v4_normalized_log2_AXN": old["normalized_log2_AXN"],
        "v5_gain_bits": old["normalized_log2_AXN"] - lg2(normalized),
    }


def main() -> None:
    references = {"I": 143, "III": 207, "V": 272}
    report: dict[str, Any] = {
        "constants": {
            "q": Q,
            "rho": [RHO.numerator, RHO.denominator],
            "g1_AXN": G1_AXN,
            "g2_pointwise_multiplication_AXN": G2_AXN,
            "g4_pointwise_multiplication_AXN": G4_AXN,
            "previous_g4_AXN": base.G4_AXN,
            "orbit_arithmetic_shift_bits": ORBIT_ARITH_SHIFT,
            "spectrum_excess_bits": SPECTRUM_EXCESS_BITS,
        },
        "homotopy_upgrade": {
            "success_at_extension_degree_r": "1-B^2/|A|^r",
            "optimized_multiplier": "(2r-1)/(1-B^2/|A|^r)",
            "old_multiplier": "(8/7)(2r_old-1), with |A|^r_old >= 8B^2",
        },
        "field_circuits": _CIRCUITS,
        "rows": {},
        "orbit_complete": {level: orbit_v5(level) for level in ("I", "III", "V")},
        "zero_offset_level3b": zero_offset_level3b_orbit_v5(),
        "all_nine_fastest_theorem_backed": {},
        "all_nine_core_complete": {},
        "combined_repair_floors": {},
        "zero_offset_square_floors": {},
        "random_cut_one_core": {},
    }

    l2_by_level: dict[str, dict[str, Any]] = {}
    for level, params, reference in base.L2_ROWS:
        diag = l2_diagonal_v5(level, params, reference)
        complete = l2_complete_v5(level, params, reference)
        fastest = min(diag["normalized_log2_AXN"], complete["normalized_log2_AXN"])
        row = {
            "level": level,
            "parameters": list(params),
            "diagonal": diag,
            "complete_two_block": complete,
            "fastest_route": diag["route"] if diag["normalized_log2_AXN"] <= complete["normalized_log2_AXN"] else complete["route"],
            "fastest_normalized_log2_AXN": fastest,
            "fastest_headroom": reference - fastest,
        }
        report["rows"][str(params)] = row
        l2_by_level[level] = row
        report["combined_repair_floors"][str(params)] = base.combined_repair_floor(params)
        report["zero_offset_square_floors"][str(params)] = v4.zero_offset_square_floor(params)

    l4_fast: dict[str, list[float]] = {"I": [], "III": [], "V": []}
    l4_core: dict[str, list[float]] = {"I": [], "III": [], "V": []}
    for name, params, reference in base.L4_ROWS:
        level = name.split("-")[0]
        compact = l4_channel_frontier_v5(name, params, reference)
        if name == "III-b":
            core = report["zero_offset_level3b"]["combined_log2_AXN"]
            core_route = "zero-offset filtered, projection-complete 16-orbit sweep"
        else:
            core = report["orbit_complete"][level]["v5_projection_complete_log2_AXN"]
            core_route = "projection-complete 16-orbit sweep"
        row = {
            "name": name,
            "level": level,
            "parameters": list(params),
            "fast_whole_channel_sharp": compact["fast"],
            "low_output_plus7_whole_channel_sharp": compact["low_output_plus7"],
            "whole_channel_sharp_profiles_tested": compact["profiles_tested"],
            "core_complete_route": core_route,
            "core_complete_normalized_log2_AXN": core,
            "core_complete_headroom": reference - core,
        }
        report["rows"][str(params)] = row
        l4_fast[level].append(compact["fast"]["normalized_log2_AXN"])
        l4_core[level].append(core)
        report["combined_repair_floors"][str(params)] = base.combined_repair_floor(params)
        report["zero_offset_square_floors"][str(params)] = v4.zero_offset_square_floor(params)
        if params[3] == 4:
            report["random_cut_one_core"][name] = random_cut_row_v5(name, params, reference)

    for level in ("I", "III", "V"):
        fast = max(l2_by_level[level]["fastest_normalized_log2_AXN"], max(l4_fast[level]))
        core = max(l2_by_level[level]["complete_two_block"]["normalized_log2_AXN"], max(l4_core[level]))
        report["all_nine_fastest_theorem_backed"][level] = {
            "worst_normalized_log2_AXN": fast,
            "headroom": references[level] - fast,
            "plus7_headroom": references[level] - fast - 7,
        }
        report["all_nine_core_complete"][level] = {
            "worst_normalized_log2_AXN": core,
            "headroom": references[level] - core,
            "plus7_headroom": references[level] - core - 7,
        }

    # Structural and numerical sanity checks.  Exact expected values are
    # recorded after the first independently checked run below.
    assert G2_AXN == report["field_circuits"]["F19_2"]["multiplication"]
    assert G4_AXN == report["field_circuits"]["F19_4"]["multiplication"]
    for level in ("I", "III", "V"):
        assert report["all_nine_fastest_theorem_backed"][level]["plus7_headroom"] > 0
    floor_incs = [x["combined_increase_percent"] for x in report["combined_repair_floors"].values()]
    zero_incs = [x["increase_percent"] for x in report["zero_offset_square_floors"].values()]
    assert abs(min(floor_incs) - 45.0) < 1e-12
    assert abs(max(floor_incs) - 60.71428571428572) < 1e-12
    assert abs(min(zero_incs) - 96.875) < 1e-12
    assert abs(max(zero_incs) - 128.0) < 1e-12

    out = HERE / "snova_master_upgrade_v5.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print("SNOVA master-upgrade v5 ledger")
    print("=" * 100)
    print(f"Arithmetic: pointwise G2={G2_AXN}, G4={G4_AXN} (old G4={base.G4_AXN})")
    print("\nell=2 direct-field rows:")
    for level in ("I", "III", "V"):
        row = l2_by_level[level]
        d = row["diagonal"]
        c = row["complete_two_block"]
        print(
            f"  Level {level}: diagonal 2^{d['normalized_log2_AXN']:.6f} "
            f"(gain {d['v5_gain_bits']:.6f}); complete 2^{c['normalized_log2_AXN']:.6f} "
            f"(gain {c['v5_gain_bits']:.6f})"
        )
    print("\nell=4 fast whole-channel rows (near-baseline spectrum):")
    for name, params, _ in base.L4_ROWS:
        row = report["rows"][str(params)]["fast_whole_channel_sharp"]
        print(
            f"  {name}: profile {tuple(row['profile'])}, 2^{row['normalized_log2_AXN']:.6f}, "
            f"same-profile gain {row['v5_same_profile_gain_bits']:.6f}"
        )
    print("\nAll-nine frontiers:")
    for level in ("I", "III", "V"):
        f = report["all_nine_fastest_theorem_backed"][level]
        c = report["all_nine_core_complete"][level]
        print(
            f"  Level {level}: fastest 2^{f['worst_normalized_log2_AXN']:.6f} "
            f"(headroom {f['headroom']:.6f}); core-complete "
            f"2^{c['worst_normalized_log2_AXN']:.6f} (headroom {c['headroom']:.6f})"
        )
    print("\nRandom-cut one-core rows:")
    for name, row in report["random_cut_one_core"].items():
        print(f"  {name}: 2^{row['normalized_log2_AXN']:.6f} (gain {row['v5_gain_bits']:.6f})")
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
