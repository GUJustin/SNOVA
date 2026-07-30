#!/usr/bin/env python3
"""Exact direct-kernel bounds for recommended original-family core Jacobians.

At the fixed Frobenius-descent point used by the witness certificate, the first
row/column data of every public base form decompose, for each off-diagonal
A-coordinate pair, into an invertible F_q-linear image of 16 independent
mod-q expansion symbols.  For a nonzero left combination of the selected rows
belonging to one base form, the row-template map is injective, so at every one
of the s-1 off-diagonal coordinate positions at least one nonzero A-linear
constraint is imposed.  Each such constraint has F_q-rank four and therefore
probability at most rho^4.  Different coordinate positions and public forms
use disjoint source symbols.

Grouping a putative left-kernel vector by its active public-form blocks gives
  delta_core <= prod_i (1 + (|A|^t_i-1) rho^(4(s-1))) - 1,
where t_i is the number of selected core rows drawn from public form i.
"""
from __future__ import annotations

import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

sys.set_int_max_str_digits(0)
getcontext().prec = 140
ROOT = Path(__file__).resolve().parent
Q = 19
D = 4
A_SIZE = Q**D
RHO = Fraction(14, 256)


def log2_fraction(x: Fraction) -> Decimal:
    if x <= 0:
        raise ValueError("nonpositive logarithm")
    return (Decimal(x.numerator).ln() - Decimal(x.denominator).ln()) / Decimal(2).ln()


def frac_obj(x: Fraction) -> dict[str, object]:
    return {
        "numerator": x.numerator,
        "denominator": x.denominator,
        "decimal": float(Decimal(x.numerator) / Decimal(x.denominator)),
        "log2": float(log2_fraction(x)),
    }


def selected_profiles() -> dict[str, dict[str, object]]:
    data = json.loads((ROOT / "eigenblock_core_certificate.json").read_text())
    return {p["level"]: p for p in data["profiles"] if p["role"] == "recommended"}


def spectrum_by_level() -> dict[str, dict[str, object]]:
    data = json.loads((ROOT / "eigenblock_core_certificate.json").read_text())
    return {x["level"]: x for x in data["seeded_spectrum_bounds"]}


def structural_rows() -> list[dict[str, object]]:
    return json.loads((ROOT / "structural_preflight_probability.json").read_text())["results"]


def as_fraction(obj: dict[str, object]) -> Fraction:
    return Fraction(int(obj["numerator"]), int(obj["denominator"]))


def form_row_counts(p: dict[str, object]) -> list[int]:
    par = p["parameters"]
    m = int(par["m1"])
    counts = [0] * m
    for c in list(par["loop_counts"]) + list(par["edge_counts_lexicographic"]):
        for i in range(int(c)):
            counts[i] += 1
    if sum(counts) != int(par["variables_N"]):
        raise AssertionError((p["key"], counts, sum(counts), par["variables_N"]))
    return counts


def core_failure(p: dict[str, object]) -> tuple[list[int], Fraction]:
    par = p["parameters"]
    s = int(par["slice_dimension_s"])
    t = form_row_counts(p)
    per_active = RHO ** (D * (s - 1))
    product = Fraction(1, 1)
    for ti in t:
        product *= 1 + (A_SIZE**ti - 1) * per_active
    return t, product - 1


def main() -> None:
    prof = selected_profiles()
    spec = spectrum_by_level()
    level_core: dict[str, dict[str, object]] = {}
    for level in ("I", "III", "V"):
        p = prof[level]
        counts, delta = core_failure(p)
        success = 1 - delta
        level_core[level] = {
            "profile": p["key"],
            "s": p["parameters"]["slice_dimension_s"],
            "core_rows_by_form": counts,
            "active_form_constraint_exponent": D * (int(p["parameters"]["slice_dimension_s"]) - 1),
            "failure": frac_obj(delta),
            "success": frac_obj(success),
        }

    rows = []
    for st in structural_rows():
        level = str(st["level"])
        p = prof[level]
        delta = as_fraction(level_core[level]["failure"])
        eps_compact = as_fraction(spec[level]["compact"])
        eps_robust = as_fraction(spec[level]["robust"])
        pre = as_fraction(st["complete_cross_preflight_success"])
        row = {
            "shape": st["key"],
            "level": level,
            "profile": p["key"],
            "core_failure": frac_obj(delta),
            "cross_preflight": frac_obj(pre),
            "costs": {},
        }
        reference = int(p["parameters"]["reference_bits"])
        for mode, eps in (("compact", eps_compact), ("robust", eps_robust)):
            internal = 1 - delta - eps
            if internal <= 0:
                raise AssertionError((level, mode, internal))
            density = pre * internal
            penalty = -log2_fraction(density)
            per_key = Decimal(str(p["costs"][mode]["total_bits"]))
            normalized = per_key + penalty
            normalized_256 = per_key + Decimal(8) + penalty
            row["costs"][mode] = {
                "spectrum_failure": frac_obj(eps),
                "internal_success": frac_obj(internal),
                "vulnerable_key_density": frac_obj(density),
                "density_penalty_bits": float(penalty),
                "per_certified_key_bits": float(per_key),
                "normalized_bits": float(normalized),
                "normalized_margin_bits": float(Decimal(reference) - normalized),
                "normalized_kappa_256_bits": float(normalized_256),
                "normalized_kappa_256_margin_bits": float(Decimal(reference) - normalized_256),
            }
        rows.append(row)

    data = {
        "field": {"q": Q, "A_size": A_SIZE, "rho": "14/256"},
        "proof_ledger": {
            "per_active_form_constraints": "4(s-1) independent F_19 constraints",
            "bound": "prod_i(1+(|A|^t_i-1)*rho^(4(s-1)))-1",
            "independence": "off-diagonal A-coordinate pairs and public base forms use disjoint expansion symbols",
            "combination_with_spectrum": "union bound inside the internal V' x V' coefficient region",
            "combination_with_structural_preflight": "multiplication, because the coordinate regions are disjoint",
        },
        "levels": level_core,
        "results": rows,
    }
    (ROOT / "core_jacobian_probability.json").write_text(json.dumps(data, indent=2) + "\n")

    lines = [
        "# Direct core-Jacobian probability certificate",
        "",
        "For the recommended original-family cores, a nonzero left-kernel combination supported on one public base form forces at least `s-1` independent A-valued equations on disjoint off-diagonal coordinate-pair blocks.  Each A-valued equation has F_19-rank four.  With `rho=14/256`,",
        "",
        "```text",
        "delta_core <= product_i (1 + (|A|^t_i-1) rho^(4(s-1))) - 1.",
        "```",
        "",
        "Here `t_i` is the exact number of selected core rows drawn from public form `i`.  The row-template map is injective: loop multipliers appear in distinct diagonal blocks, adjacent-edge multipliers appear in an exposed first/last block, and the opposite-edge multiplier appears in its own H component.  Thus every active form contributes a nonzero coefficient map.",
        "",
        "| Level | s | rows by public form | log2 core failure | core success |",
        "|:--:|--:|:--|--:|--:|",
    ]
    for level in ("I", "III", "V"):
        z = level_core[level]
        lines.append(f"| {level} | {z['s']} | `{z['core_rows_by_form']}` | {z['failure']['log2']:.2f} | {z['success']['decimal']:.15f} |")
    lines += [
        "",
        "Combining this internal-core event with the exact aggregate-spectrum tail by a union bound, and multiplying by the disjoint structural-preflight probability, gives:",
        "",
        "| Shape | robust per certified key | vulnerable-key density | normalized robust | margin | normalized robust + 2^8 | margin |",
        "|:--:|--:|--:|--:|--:|--:|--:|",
    ]
    for r in rows:
        c = r["costs"]["robust"]
        lines.append(
            f"| {r['shape']} | {c['per_certified_key_bits']:.3f} | {c['vulnerable_key_density']['decimal']:.12f} | "
            f"{c['normalized_bits']:.3f} | {c['normalized_margin_bits']:.3f} | "
            f"{c['normalized_kappa_256_bits']:.3f} | {c['normalized_kappa_256_margin_bits']:.3f} |"
        )
    lines += [
        "",
        "The normalized ledger is a time-success ratio over a randomly generated public key in the pinned random-XOF idealization.  It does not mean that an adversary can resample a fixed target key.  The `+2^8` column is an exposed implementation-sensitivity point, not a proof that the suppressed homotopy factor is at most 256.",
        "",
    ]
    (ROOT / "CORE_JACOBIAN_PROBABILITY_CERTIFICATE.md").write_text("\n".join(lines))

    for level in ("I", "III", "V"):
        z = level_core[level]
        print(level, z["failure"]["log2"], z["success"]["decimal"])
    for r in rows[::2]:
        c = r["costs"]["robust"]
        print(r["level"], c["normalized_bits"], c["normalized_kappa_256_bits"])


if __name__ == "__main__":
    main()
