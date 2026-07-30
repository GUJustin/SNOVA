#!/usr/bin/env python3
"""Version-4 numerical ledger for the SNOVA master-theorem upgrade.

This wrapper reuses the exact arithmetic routines in
``snova_master_upgrade_numbers.py`` and adds the corrections/strengthenings
that the consolidated manuscript insertion needs:

1. the random within-family projection factor that turns a family-count
   pattern into a concrete nonsingular core;
2. the Level-III dimension correction: the (38,5,4,5) row cannot use an
   h=68 stable slice, but a zero-offset chosen-message consistency filter
   restores the full A-linear domain and makes the h=68 orbit sweep valid;
3. compact eta=1/19 whole-channel and orbit-complete frontiers for ell=4;
4. a chosen-message square-persistence theorem and its stronger dimension-only
   repair floor; and
5. both the fastest theorem-backed and the core-complete all-nine frontiers.

All newly derived probabilities are exact Fractions until logarithms are
printed.  The three pre-existing orbit-sweep exponents are pinned artifact
inputs from the manuscript's exhaustive 16-orbit certificate; this script
only applies the newly proved projection-success correction to them.
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
BASE_PATH = HERE / "snova_master_upgrade_numbers.py"
spec = importlib.util.spec_from_file_location("snova_base", BASE_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load {BASE_PATH}")
base = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = base
spec.loader.exec_module(base)

# Pinned outputs of the manuscript's exact 16-orbit enumeration before the
# random within-family projection correction.
ORBIT_ROBUST_NORMALIZED = {
    "I": 138.40973,
    "III": 179.84405,
    "V": 220.91810,
}
ORBIT_H = {"I": 48, "III": 68, "V": 88}


def projection_overhead_bits(h: int) -> float:
    chi = 1.0 - h / (base.Q ** 4)
    if chi <= 0:
        raise ValueError("projection field too small")
    return -math.log2(chi)


def corrected_orbit(level: str) -> dict[str, Any]:
    h = ORBIT_H[level]
    overhead = projection_overhead_bits(h)
    normalized = ORBIT_ROBUST_NORMALIZED[level] + overhead
    reference = {"I": 143, "III": 207, "V": 272}[level]
    return {
        "h": h,
        "projection_field_size": base.Q ** 4,
        "projection_success": 1.0 - h / (base.Q ** 4),
        "projection_overhead_bits": overhead,
        "robust_base_normalized_log2_AXN": ORBIT_ROBUST_NORMALIZED[level],
        "corrected_normalized_log2_AXN": normalized,
        "headroom": reference - normalized,
    }


def _orbit_acceptance(level: str) -> Fraction:
    # The pinned per-level orbit ledger conservatively charges the square-row
    # acceptance density for both public shapes at that level.
    n = {"I": 33, "III": 47, "V": 59}[level]
    return (1 - Fraction(1, base.Q**4)) ** n


def compact_orbit(level: str) -> dict[str, Any]:
    """Convert the pinned eta=1/2 orbit ledger to eta=1/q exactly."""
    h = ORBIT_H[level]
    alpha = _orbit_acceptance(level)
    eta_robust = Fraction(1, 2)
    eta_compact = Fraction(1, base.Q)
    p_robust = (alpha - eta_robust / (base.Q - 1)) ** 2 / (1 + eta_robust)
    p_compact = (alpha - eta_compact / (base.Q - 1)) ** 2 / (1 + eta_compact)
    work_shift = base.lg2(p_robust / p_compact)
    overhead = projection_overhead_bits(h)
    normalized = ORBIT_ROBUST_NORMALIZED[level] + work_shift + overhead
    reference = {"I": 143, "III": 207, "V": 272}[level]
    return {
        "h": h,
        "eta": [1, base.Q],
        "projection_field_size": base.Q**4,
        "projection_success": 1.0 - h / (base.Q**4),
        "projection_overhead_bits": overhead,
        "compact_work_shift_bits": work_shift,
        "robust_base_normalized_log2_AXN": ORBIT_ROBUST_NORMALIZED[level],
        "compact_corrected_normalized_log2_AXN": normalized,
        "robust_corrected_normalized_log2_AXN": ORBIT_ROBUST_NORMALIZED[level] + overhead,
        "headroom": reference - normalized,
    }


def l4_channel_profile_eta(
    name: str,
    params: tuple[int, int, int, int],
    reference: int,
    s: int,
    a: int,
    b: int,
    eta: Fraction,
) -> dict[str, Any]:
    """Exact whole-channel ledger at a caller-specified spectrum threshold."""
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
    delta = Fraction(a + base.Q * b, base.Q**4)
    if alpha <= delta:
        raise ValueError("singular-density charge exhausts accepted density")
    mu = Fraction(1, base.Q ** (K - h))
    p_root = mu * (alpha - delta) ** 2 / (1 + eta)

    B = 2**a * (base.Q + 1) ** b
    Bplus_f = Fraction(B, 1) * (1 + Fraction(a, 2) + Fraction(b, base.Q + 1))
    if Bplus_f.denominator != 1:
        raise ValueError("nonintegral companion Bezout quantity")
    Bplus = Bplus_f.numerator
    D = math.comb(s + 1, 2)
    slp0 = 0 if a == 0 else D + 2 * a * (D + s + 1)
    slp1 = 0 if b == 0 else 6 * s + s * s + 2 * b * (s * s + 2 * s + 1)
    degree_sum = 2 * a + (base.Q + 1) * b
    r_ext = base.ceil_extension_degree(base.Q**4, max(degree_sum, 8 * B * B))
    H = slp0 + slp1 + degree_sum + s * s
    work = (
        Fraction(8, 7)
        / p_root
        * (2 * r_ext - 1)
        * base.G4_AXN
        * B
        * Bplus
        * H
        * s
    )

    p_pre = base.structural_preflight(4, v, o, m1, M)
    A_size = base.Q**4
    projective_lines = Fraction(A_size**s - 1, A_size - 1)
    jac_fail = projective_lines * base.RHO ** (4 * s)
    eps = base.spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta)
    p_key = p_pre * (1 - jac_fail - eps)
    if p_key <= 0:
        raise ValueError("nonpositive l4 key density")
    normalized = work / p_key
    output_bits = (s + 1) * B * (4 * r_ext) * math.log2(base.Q)
    return {
        "name": name,
        "parameters": list(params),
        "profile": [s, a, b],
        "eta": [eta.numerator, eta.denominator],
        "structural_probability": float(p_pre),
        "jacobian_failure_log2": base.lg2(jac_fail),
        "spectrum_failure_log2": base.lg2(eps),
        "combined_key_probability": float(p_key),
        "root_probability_log2": base.lg2(p_root),
        "B_log2": base.lg2(B),
        "Bplus_log2": base.lg2(Bplus),
        "extension_degree_over_F19_4": r_ext,
        "per_key_log2_AXN": base.lg2(work),
        "normalized_log2_AXN": base.lg2(normalized),
        "plus7_log2_AXN": base.lg2(normalized) + 7,
        "plus7_headroom": reference - base.lg2(normalized) - 7,
        "output_GiB": output_bits / 8 / 2**30,
    }


def l4_channel_frontier_eta(
    name: str,
    params: tuple[int, int, int, int],
    reference: int,
    eta: Fraction,
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
            rows.append(l4_channel_profile_eta(name, params, reference, s, a, b, eta))
    fast = min(rows, key=lambda z: z["normalized_log2_AXN"])
    eligible = [z for z in rows if z["plus7_log2_AXN"] <= reference]
    low = min(eligible, key=lambda z: z["output_GiB"])
    return {"fast": fast, "low_output_plus7": low, "profiles_tested": len(rows)}


def zero_offset_spectrum_tail(*, d: int, v: int, K: int, h: int, eta: Fraction) -> Fraction:
    """Excess-spectrum tail with no affine constraint: D0=V'."""
    t = d * (v - 2)
    baseline = Fraction(base.Q**h - 1, base.Q**K)
    if eta <= baseline:
        raise ValueError("eta must exceed the full-rank baseline")
    return (
        Fraction(base.Q**h - 1, 1)
        * (1 - Fraction(1, base.Q**K))
        * base.RHO**t
        / (eta - baseline)
    )


def quotient_rank_probability(K: int) -> Fraction:
    """Biased-product lower bound for an M-by-K outer image to have rank K."""
    p = Fraction(1, 1)
    for c in range(1, K + 1):
        p *= 1 - base.RHO**c
    return p


def zero_offset_level3b_orbit() -> dict[str, Any]:
    """Target-filtered h=68 orbit route for (38,5,4,5)."""
    params = (38, 5, 4, 5)
    v, o, d, rcols = params
    M, K, h = o * rcols * d, 70, 68
    eta = Fraction(1, base.Q)
    eps = zero_offset_spectrum_tail(d=d, v=v, K=K, h=h, eta=eta)
    p_rank = quotient_rank_probability(K)
    projection = projection_overhead_bits(h)
    # The exhaustive artifact's per-certified-key eta=1/2 ledger.  It
    # conservatively charges the square-row acceptance density at Level III;
    # retaining that acceptance charge for the rectangular row is safe.
    per_key_robust = 179.75835
    alpha = _orbit_acceptance("III")
    eta_robust = Fraction(1, 2)
    p_robust = (alpha - eta_robust / (base.Q - 1)) ** 2 / (1 + eta_robust)
    p_compact = (alpha - eta / (base.Q - 1)) ** 2 / (1 + eta)
    compact_shift = base.lg2(p_robust / p_compact)
    normalized = (
        per_key_robust
        + compact_shift
        + projection
        - base.lg2(p_rank * (1 - eps))
    )
    target_filter_exponent = (M - K) * math.log2(base.Q)
    target_generation_exponent = target_filter_exponent + 32.0
    combined = max(normalized, target_generation_exponent) + math.log2(
        1 + 2 ** (-abs(normalized - target_generation_exponent))
    )
    return {
        "parameters": list(params),
        "M": M,
        "K": K,
        "h": h,
        "A_dimension": h // d,
        "eta": [eta.numerator, eta.denominator],
        "quotient_rank_probability": float(p_rank),
        "spectrum_failure_log2": base.lg2(eps),
        "projection_overhead_bits": projection,
        "compact_work_shift_bits": compact_shift,
        "solve_normalized_log2_AXN": normalized,
        "target_filter_log2_trials": target_filter_exponent,
        "target_generation_envelope_log2_AXN": target_generation_exponent,
        "combined_log2_AXN": combined,
        "headroom": 207 - combined,
    }


def zero_offset_square_floor(params: tuple[int, int, int, int]) -> dict[str, Any]:
    """Necessary coupling floor to eliminate the chosen-message square by dimensions."""
    v, o, d, rcols = params
    M = o * rcols * d
    C = math.comb(d + 1, 2)
    m_required = math.floor(d * (v - 1) / C) + 1
    x_required = d * (m_required - 1) + 1
    M_required = d * x_required
    return {
        "current_M": M,
        "required_base_forms": m_required,
        "required_x": x_required,
        "dimension_only_floor_M": M_required,
        "increase_percent": 100 * (M_required / M - 1),
    }



def _mul_theta_state(p0, p1, terms, cap):
    from collections import defaultdict
    n0, n1 = defaultdict(int), defaultdict(int)
    for exp, coeff in p0.items():
        for delta, dcoeff in terms:
            new = tuple(x + y for x, y in zip(exp, delta))
            if all(x <= cap for x in new):
                n0[new] += coeff * dcoeff
    for exp, coeff in p1.items():
        for delta, dcoeff in terms:
            new = tuple(x + y for x, y in zip(exp, delta))
            if all(x <= cap for x in new):
                n1[new] += coeff * dcoeff
    for exp, coeff in p0.items():
        n1[exp] += coeff
    return dict(n0), dict(n1)


def random_cut_bezout(m: int, s: int) -> tuple[int, int]:
    factors = []
    for a in range(4):
        delta = tuple(1 if i == a else 0 for i in range(4))
        factors.extend([[(delta, 2)]] * m)
    for a in range(4):
        for c in range(a + 1, 4):
            da = tuple(1 if i == a else 0 for i in range(4))
            dc = tuple(1 if i == c else 0 for i in range(4))
            factors.extend([[(da, 1), (dc, 1)]] * m)
    cut = [(tuple(1 if i == a else 0 for i in range(4)), 1) for a in range(4)]
    factors.extend([cut, cut])
    p0, p1 = {(0, 0, 0, 0): 1}, {}
    for factor in factors:
        p0, p1 = _mul_theta_state(p0, p1, factor, s)
    target = (s, s, s, s)
    B = p0[target]
    return B, B + sum(p1.values())


def random_cut_row(name: str, params: tuple[int, int, int, int], reference: int) -> dict[str, Any]:
    from fractions import Fraction
    import math
    v, o, d, rcols = params
    if d != 4 or rcols != 4:
        raise ValueError('random-cut row requires square ell=4 shape')
    m = (o * rcols + d - 1) // d
    M, K = o * rcols * d, 10 * m
    h, s = K + 2, (K + 2) // 4
    if 4 * s != h or s > v + o - 4 * m:
        raise ValueError('required stable slice unavailable')
    B, Bplus = random_cut_bezout(m, s)
    alpha = base.l4_acceptance(v + o, rcols)
    theta = Fraction(base.Q**h - 1, base.Q**h) + Fraction(1, base.Q)
    p_root = (alpha - theta / (base.Q - 1)) ** 2 / (1 + theta)
    D = math.comb(s + 1, 2)
    slp = 4 * D + 6 * s * s + 2 * (4 * m * (D + s + 1) + 6 * m * (s + 1) ** 2) + 4 * (4 * s + 1)
    degree_sum = 20 * m + 8
    r_ext = base.ceil_extension_degree(base.Q**4, max(5 * m + 2, 8 * B * B))
    H = slp + degree_sum + h * h
    work = Fraction(8, 7) / p_root * (2 * r_ext - 1) * base.G4_AXN * B * Bplus * H * h
    p_pre = base.structural_preflight(4, v, o, m, M)
    eta_bar = base.Q**2 * theta
    eps = base.spectrum_excess_tail(d=4, v=v, M=M, K=K, h=h, eta=eta_bar)
    normalized = work / (p_pre * (1 - eps))
    return {
        'name': name,
        'parameters': list(params),
        'h': h,
        'K': K,
        's': s,
        'B_log2': base.lg2(B),
        'Bplus_log2': base.lg2(Bplus),
        'root_probability': float(p_root),
        'spectrum_failure_log2': base.lg2(eps),
        'normalized_log2_AXN': base.lg2(normalized),
        'headroom': reference - base.lg2(normalized),
    }

def main() -> None:
    references = {"I": 143, "III": 207, "V": 272}
    report: dict[str, Any] = {
        "constants": {
            "q": base.Q,
            "rho": [base.RHO.numerator, base.RHO.denominator],
            "g4_AXN": base.G4_AXN,
        },
        "orbit_projection_robust": {
            level: corrected_orbit(level) for level in ("I", "III", "V")
        },
        "orbit_compact": {
            level: compact_orbit(level) for level in ("I", "III", "V")
        },
        "zero_offset_level3b": zero_offset_level3b_orbit(),
        "rows": {},
        "all_nine_fastest_theorem_backed": {},
        "all_nine_core_complete": {},
        "combined_repair_floors": {},
        "zero_offset_square_floors": {},
        "random_cut_one_core": {},
        "dimension_correction": {
            "row": [38, 5, 4, 5],
            "claimed_stable_orbit_h": 68,
            "required_stable_A_dimension": 17,
            "guaranteed_stable_A_dimension": 15,
            "valid_replacement": (
                "zero-offset chosen-message consistency filter, followed by "
                "an h=68 A-linear slice in the full public-vinegar domain"
            ),
            "ordinary_square_fallback": "complete ordinary 70-by-70 square system",
        },
    }

    # ell=2 rows.
    l2_by_level: dict[str, dict[str, Any]] = {}
    for level, params, reference in base.L2_ROWS:
        diag = base.l2_diagonal(level, params, reference)
        complete = base.l2_complete(level, params, reference)
        ordinary = base.l2_ordinary(level, params, reference)
        fastest = min(diag.normalized_log2_AXN, complete.normalized_log2_AXN)
        fastest_route = (
            "one-eigenblock diagonal core"
            if diag.normalized_log2_AXN <= complete.normalized_log2_AXN
            else "complete two-block system"
        )
        row = {
            "parameters": list(params),
            "level": level,
            "fastest_route": fastest_route,
            "fastest_normalized_log2_AXN": fastest,
            "fastest_headroom": reference - fastest,
            "diagonal": base.asdict(diag),
            "complete_two_block": base.asdict(complete),
            "ordinary_square": base.asdict(ordinary),
        }
        key = str(params)
        report["rows"][key] = row
        l2_by_level[level] = row
        report["combined_repair_floors"][key] = base.combined_repair_floor(params)
        report["zero_offset_square_floors"][key] = zero_offset_square_floor(params)

    # ell=4 rows.  Headline whole-channel rows use eta=1/q; eta=1/2 is
    # retained as a joint-stress sensitivity ledger.  The core-complete
    # III-b row uses the zero-offset chosen-message filter.
    l4_fast_by_level: dict[str, list[float]] = {"I": [], "III": [], "V": []}
    l4_core_by_level: dict[str, list[float]] = {"I": [], "III": [], "V": []}
    zero_IIIb = report["zero_offset_level3b"]
    for name, params, reference in base.L4_ROWS:
        level = name.split("-")[0]
        compact_frontier = l4_channel_frontier_eta(
            name, params, reference, Fraction(1, base.Q)
        )
        robust_frontier = l4_channel_frontier_eta(
            name, params, reference, Fraction(1, 2)
        )
        fast = compact_frontier["fast"]
        low = compact_frontier["low_output_plus7"]

        if name == "III-b":
            core_route = "zero-offset filtered projection-complete 16-orbit sweep"
            core_norm = zero_IIIb["combined_log2_AXN"]
            ordinary = base.l4_level3_ordinary(params, reference)
            core_extra = {
                "zero_offset_orbit": zero_IIIb,
                "ordinary_square_fallback": base.asdict(ordinary),
            }
        else:
            orbit = compact_orbit(level)
            core_route = "compact projection-complete 16-orbit sweep"
            core_norm = orbit["compact_corrected_normalized_log2_AXN"]
            core_extra = {"orbit": orbit}

        row = {
            "name": name,
            "parameters": list(params),
            "level": level,
            "fast_whole_channel_compact": fast,
            "low_output_plus7_whole_channel_compact": low,
            "whole_channel_compact_profiles_tested": compact_frontier["profiles_tested"],
            "whole_channel_robust": robust_frontier,
            "core_complete_route": core_route,
            "core_complete_normalized_log2_AXN": core_norm,
            "core_complete_headroom": reference - core_norm,
            **core_extra,
        }
        key = str(params)
        report["rows"][key] = row
        l4_fast_by_level[level].append(fast["normalized_log2_AXN"])
        l4_core_by_level[level].append(core_norm)
        report["combined_repair_floors"][key] = base.combined_repair_floor(params)
        report["zero_offset_square_floors"][key] = zero_offset_square_floor(params)

    for name, params, reference in base.L4_ROWS:
        if params[3] == 4:
            report["random_cut_one_core"][name] = random_cut_row(name, params, reference)

    for level in ("I", "III", "V"):
        fastest = max(
            l2_by_level[level]["fastest_normalized_log2_AXN"],
            max(l4_fast_by_level[level]),
        )
        report["all_nine_fastest_theorem_backed"][level] = {
            "worst_normalized_log2_AXN": fastest,
            "headroom": references[level] - fastest,
            "plus7_headroom": references[level] - fastest - 7,
        }

        core = max(
            l2_by_level[level]["complete_two_block"]["normalized_log2_AXN"],
            max(l4_core_by_level[level]),
        )
        report["all_nine_core_complete"][level] = {
            "worst_normalized_log2_AXN": core,
            "headroom": references[level] - core,
        }

    expected_fast = {
        "I": 134.83581037063414,
        "III": 186.31047399951058,
        "V": 236.37390592635347,
    }
    # Values below are filled from the exact run and asserted to guard the
    # paper tables.  Level III/V are dominated by the complete ell=2 system.
    expected_core = {
        "I": 137.82644602131035,
        "III": 186.31047399951058,
        "V": 236.37390592635347,
    }
    for level in ("I", "III", "V"):
        got_fast = report["all_nine_fastest_theorem_backed"][level]["worst_normalized_log2_AXN"]
        got_core = report["all_nine_core_complete"][level]["worst_normalized_log2_AXN"]
        assert abs(got_fast - expected_fast[level]) < 1e-9
        assert abs(got_core - expected_core[level]) < 1e-9
        assert expected_fast[level] + 7 < references[level]

    expected_cut = {
        "I-a": 138.59677640000518,
        "III-a": 180.17215997124003,
        "V-a": 221.187920052369,
    }
    for name, value in expected_cut.items():
        assert abs(report["random_cut_one_core"][name]["normalized_log2_AXN"] - value) < 1e-9

    increases = [
        row["combined_increase_percent"]
        for row in report["combined_repair_floors"].values()
    ]
    assert abs(min(increases) - 45.0) < 1e-12
    assert abs(max(increases) - 60.71428571428572) < 1e-12
    zero_increases = [
        row["increase_percent"] for row in report["zero_offset_square_floors"].values()
    ]
    assert abs(min(zero_increases) - 96.875) < 1e-12
    assert abs(max(zero_increases) - 128.0) < 1e-12
    assert zero_IIIb["target_generation_envelope_log2_AXN"] + 19 < zero_IIIb["solve_normalized_log2_AXN"]

    out = HERE / "snova_master_upgrade_v4.json"
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    print("SNOVA master-upgrade v4 ledger")
    print("=" * 100)
    print("Compact orbit-complete rows:")
    for level, row in report["orbit_compact"].items():
        print(
            f"  Level {level}: eta=1/19 shift {row['compact_work_shift_bits']:.6f} bits, "
            f"projection +{row['projection_overhead_bits']:.6f}, "
            f"normalized 2^{row['compact_corrected_normalized_log2_AXN']:.6f}"
        )
    z = report["zero_offset_level3b"]
    print(
        "\nIII-b zero-offset replacement: "
        f"solve 2^{z['solve_normalized_log2_AXN']:.6f}, "
        f"target envelope 2^{z['target_generation_envelope_log2_AXN']:.6f}, "
        f"combined 2^{z['combined_log2_AXN']:.6f}"
    )
    print("\nAll-nine frontiers:")
    for level in ("I", "III", "V"):
        fast = report["all_nine_fastest_theorem_backed"][level]
        core = report["all_nine_core_complete"][level]
        print(
            f"  Level {level}: fastest 2^{fast['worst_normalized_log2_AXN']:.6f} "
            f"(headroom {fast['headroom']:.6f}); core-complete "
            f"2^{core['worst_normalized_log2_AXN']:.6f} "
            f"(headroom {core['headroom']:.6f})"
        )
    print("\nRandom-cut one-core rows:")
    for name, row in report["random_cut_one_core"].items():
        print(f"  {name}: 2^{row['normalized_log2_AXN']:.6f}")
    print(
        "\nDimension-only zero-offset square floor range: "
        f"{min(zero_increases):.3f}% to {max(zero_increases):.3f}%"
    )
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
