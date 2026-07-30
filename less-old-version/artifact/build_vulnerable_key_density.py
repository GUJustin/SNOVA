#!/usr/bin/env python3
"""Combine exact public-preflight, core-Jacobian, and spectrum-tail bounds.

This script certifies a lower bound on the random-XOF density of keys for
which the recommended original-family eigenblock-core attack is applicable.
All probability arithmetic is exact over Q.  Floating-point values are used
only for logarithms and presentation.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

sys.set_int_max_str_digits(1000000)

HERE = Path(__file__).resolve().parent


def frac(obj: dict[str, Any], prefix: str) -> Fraction:
    return Fraction(int(obj[f"{prefix}_numerator"]), int(obj[f"{prefix}_denominator"]))


def log2_fraction(x: Fraction) -> float:
    if x <= 0:
        raise ValueError("log2 requires a positive fraction")
    # Avoid converting very large integers directly to float.
    return math.log2(x.numerator) - math.log2(x.denominator)


def dec(x: Fraction, digits: int = 15) -> str:
    return f"{float(x):.{digits}g}"


def human_prob(x: Fraction) -> str:
    return f"{float(x):.12g}"


def load(name: str) -> dict[str, Any]:
    with (HERE / name).open("r", encoding="utf-8") as f:
        return json.load(f)


core = load("eigenblock_core_certificate.json")
pre = load("structural_preflight_probability.json")
anti = load("preflight_anticoncentration.json")

recommended_key = {"I": "L1-practical", "III": "L3-practical", "V": "L5-practical"}
profiles = {p["key"]: p for p in core["profiles"]}
anti_profiles = {p["level"]: p for p in anti["profiles"]}
spectrum = {r["level"]: r for r in core["seeded_spectrum_bounds"]}

results: list[dict[str, Any]] = []
for sh in pre["results"]:
    level = sh["level"]
    p = profiles[recommended_key[level]]
    ap = anti_profiles[level]
    sp = spectrum[level]

    cross = Fraction(
        int(sh["complete_cross_preflight_success"]["numerator"]),
        int(sh["complete_cross_preflight_success"]["denominator"]),
    )
    core_nonzero = Fraction(
        int(ap["nonzero_probability_lower_numerator"]),
        int(ap["nonzero_probability_lower_denominator"]),
    )
    eps = Fraction(int(sp["robust"]["numerator"]), int(sp["robust"]["denominator"]))
    internal = core_nonzero - eps
    if internal <= 0:
        raise AssertionError(f"nonpositive internal success lower bound for {sh['key']}")
    vulnerable = cross * internal

    robust_bits = float(p["costs"]["robust"]["total_bits"])
    density_penalty = -log2_fraction(vulnerable)
    normalized = robust_bits + density_penalty
    normalized_kappa256 = normalized + 8.0
    reference = int(p["parameters"]["reference_bits"])

    result = {
        "shape": sh["key"],
        "level": level,
        "M": sh["M"],
        "K": sh["K"],
        "profile": p["key"],
        "cross_preflight_probability": {
            "numerator": cross.numerator,
            "denominator": cross.denominator,
            "decimal": float(cross),
            "log2": log2_fraction(cross),
        },
        "core_jacobian_nonzero_probability": {
            "numerator": core_nonzero.numerator,
            "denominator": core_nonzero.denominator,
            "decimal": float(core_nonzero),
            "log2": log2_fraction(core_nonzero),
        },
        "robust_spectrum_failure_probability": {
            "numerator": eps.numerator,
            "denominator": eps.denominator,
            "decimal": float(eps),
            "log2": log2_fraction(eps),
        },
        "internal_core_and_spectrum_probability": {
            "numerator": internal.numerator,
            "denominator": internal.denominator,
            "decimal": float(internal),
            "log2": log2_fraction(internal),
        },
        "vulnerable_key_density_lower": {
            "numerator": vulnerable.numerator,
            "denominator": vulnerable.denominator,
            "decimal": float(vulnerable),
            "log2": log2_fraction(vulnerable),
            "density_penalty_bits": density_penalty,
        },
        "costs": {
            "per_vulnerable_key_robust_bits": robust_bits,
            "random_key_normalized_robust_bits": normalized,
            "random_key_normalized_robust_plus_kappa256_bits": normalized_kappa256,
            "reference_bits": reference,
            "normalized_margin_bits": reference - normalized,
            "normalized_plus_kappa256_margin_bits": reference - normalized_kappa256,
        },
    }
    results.append(result)

out = {
    "model": "random-XOF idealization with coordinatewise modulo-19 bytes",
    "rho": "14/256",
    "combination_rule": (
        "Cross-preflight byte regions are disjoint from the internal V' x V' region, "
        "so their lower bound multiplies.  The core-Jacobian event and spectrum event "
        "share the internal region, so they are combined by the union bound beta_core-epsilon."
    ),
    "results": results,
}

with (HERE / "vulnerable_key_density.json").open("w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
    f.write("\n")

lines = [
    "# Random-XOF vulnerable-key density certificate",
    "",
    "This certificate combines three exact ingredients for each official `ell=4` shape:",
    "",
    "1. the disjoint-coordinate public preflight bound (full quotient/source rank, outer-map ranks, and projection injectivity);",
    "2. the blockwise anti-concentration bound for a nonzero core-Jacobian witness polynomial; and",
    "3. the robust seeded-spectrum tail `Pr[bar_eta_L > 1/2]`.",
    "",
    "The cross-preflight coordinates are disjoint from the internal `V' x V'` coordinates, so those probabilities multiply.  The Jacobian and spectrum events share the internal coordinates, so the exact certified lower bound is",
    "",
    "```text",
    "p_vuln >= p_cross * (p_Jac - eps_spectrum).",
    "```",
    "",
    "All probability arithmetic in the JSON is exact rational arithmetic.  Logarithms below are presentation values.",
    "",
    "| Shape | p_cross | p_Jac | eps_spectrum | p_vuln lower | penalty (bits) | per-key robust | normalized robust | margin | normalized + 2^8 | margin |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
]
for r in results:
    c = r["costs"]
    lines.append(
        f"| {r['shape']} | {r['cross_preflight_probability']['decimal']:.9g} | "
        f"{r['core_jacobian_nonzero_probability']['decimal']:.9g} | "
        f"{r['robust_spectrum_failure_probability']['decimal']:.3e} | "
        f"{r['vulnerable_key_density_lower']['decimal']:.9g} | "
        f"{r['vulnerable_key_density_lower']['density_penalty_bits']:.3f} | "
        f"{c['per_vulnerable_key_robust_bits']:.3f} | "
        f"{c['random_key_normalized_robust_bits']:.3f} | "
        f"{c['normalized_margin_bits']:.3f} | "
        f"{c['random_key_normalized_robust_plus_kappa256_bits']:.3f} | "
        f"{c['normalized_plus_kappa256_margin_bits']:.3f} |"
    )

lines += [
    "",
    "## Interpretation",
    "",
    "The per-key column is the attack ledger conditional on the public certificate.  The normalized column additionally pays the inverse certified key density, equivalently measuring work per successful forgery over a random generated public key in the idealized distribution.  It is not a claim that an adversary can force a target installation to regenerate its key.",
    "",
    "The `+ 2^8` column exposes a factor-256 implementation sensitivity on the leading homotopy term.  It is a sensitivity point, not a theorem that all hidden constants are at most 256.",
]

with (HERE / "VULNERABLE_KEY_DENSITY_CERTIFICATE.md").open("w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")

print("wrote vulnerable_key_density.json and VULNERABLE_KEY_DENSITY_CERTIFICATE.md")
for r in results:
    c = r["costs"]
    print(
        r["shape"],
        "p=", f"{r['vulnerable_key_density_lower']['decimal']:.12g}",
        "penalty=", f"{r['vulnerable_key_density_lower']['density_penalty_bits']:.6f}",
        "normalized=", f"{c['random_key_normalized_robust_bits']:.6f}",
        "margin=", f"{c['normalized_margin_bits']:.6f}",
    )
