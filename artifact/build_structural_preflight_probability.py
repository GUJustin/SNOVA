#!/usr/bin/env python3
"""Exact random-XOF bounds for the structured affine/public preflights.

For q=19, ell=4, R_j=S^j, Gamma_0=1, Gamma_1=0:

* Modulo the 10-dimensional symmetric-square image of the raw (j,k)=(1,0)
  labels, the offset-linear map contains one copy of wedge^2(A) per base
  form.  For every nonzero dual alternating tuple, the coefficient map from
  the fresh W x V' block to (V')^* is surjective.  A projective union bound
  therefore certifies rank 6*m1 for this alternating affine source.

* Failure of projection injectivity is equivalent to the existence of a
  nonzero A-line Ax in J=W+O on which all base forms vanish against W.  The
  line W itself costs 10 independent symmetric coordinates per form; every
  other A-line costs 16 independent W x O coordinates per form.

* Conditional on an M-dimensional source subspace containing all K quotient
  columns, the independent outer map is invertible on that subspace with
  probability at least prod_{c=1}^M (1-rho^c).  This single event implies both
  rank Q=K and rank[Q Lambda]=M.

All probability arithmetic is exact over Q.  The only distributional input is
independence of the XOF bytes and maximum scalar atom rho=14/256.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

sys.set_int_max_str_digits(1000000)
ROOT = Path(__file__).resolve().parent
Q = 19
D = 4
A_SIZE = Q**D
RHO = Fraction(14, 256)
SHAPES = [
    {"key":"I-a","level":"I","v":28,"o":5,"r":4,"m1":5,"M":80,"K":50},
    {"key":"I-b","level":"I","v":28,"o":4,"r":5,"m1":5,"M":80,"K":50},
    {"key":"III-a","level":"III","v":40,"o":7,"r":4,"m1":7,"M":112,"K":70},
    {"key":"III-b","level":"III","v":38,"o":5,"r":5,"m1":7,"M":100,"K":70},
    {"key":"V-a","level":"V","v":50,"o":9,"r":4,"m1":9,"M":144,"K":90},
    {"key":"V-b","level":"V","v":52,"o":6,"r":6,"m1":9,"M":144,"K":90},
]


def log2f(x: Fraction) -> float:
    if x <= 0:
        raise ValueError("positive rational required")
    return math.log2(x.numerator) - math.log2(x.denominator)


def outer_success(M: int) -> Fraction:
    ans = Fraction(1, 1)
    for c in range(1, M + 1):
        ans *= 1 - RHO**c
    return ans


def certify(sh: dict[str, Any]) -> dict[str, Any]:
    m = sh["m1"]
    alternating_dim = math.comb(D, 2) * m
    vprime_dim = D * (sh["v"] - 1)
    source_projective_duals = (Q**alternating_dim - 1) // (Q - 1)
    source_failure = Fraction(source_projective_duals, 1) * RHO**vprime_dim

    # A-lines in J=W+O other than W itself.
    j_A_dim = 1 + sh["o"]
    non_W_A_lines = (A_SIZE**j_A_dim - 1) // (A_SIZE - 1) - 1
    projection_W_failure = RHO ** (math.comb(D + 1, 2) * m)
    projection_other_failure = Fraction(non_W_A_lines, 1) * RHO ** (D * D * m)
    projection_failure = projection_W_failure + projection_other_failure

    if source_failure >= 1 or projection_failure >= 1:
        raise AssertionError(sh["key"])
    source_success = 1 - source_failure
    projection_success = 1 - projection_failure
    outer = outer_success(sh["M"])
    cross = source_success * projection_success * outer

    return {
        **sh,
        "alternating_affine_rank_certified": alternating_dim,
        "affine_rank_needed": sh["M"] - sh["K"],
        "Vprime_base_dimension": vprime_dim,
        "source_projective_dual_count": source_projective_duals,
        "source_failure": {
            "numerator": source_failure.numerator,
            "denominator": source_failure.denominator,
            "decimal": float(source_failure),
            "log2": log2f(source_failure),
        },
        "source_success": {
            "numerator": source_success.numerator,
            "denominator": source_success.denominator,
            "decimal": float(source_success),
            "log2": log2f(source_success),
        },
        "projection_non_W_A_line_count": non_W_A_lines,
        "projection_W_failure": {
            "numerator": projection_W_failure.numerator,
            "denominator": projection_W_failure.denominator,
            "decimal": float(projection_W_failure),
            "log2": log2f(projection_W_failure),
        },
        "projection_other_failure": {
            "numerator": projection_other_failure.numerator,
            "denominator": projection_other_failure.denominator,
            "decimal": float(projection_other_failure),
            "log2": log2f(projection_other_failure),
        },
        "projection_failure": {
            "numerator": projection_failure.numerator,
            "denominator": projection_failure.denominator,
            "decimal": float(projection_failure),
            "log2": log2f(projection_failure),
        },
        "projection_success": {
            "numerator": projection_success.numerator,
            "denominator": projection_success.denominator,
            "decimal": float(projection_success),
            "log2": log2f(projection_success),
        },
        "outer_success": {
            "numerator": outer.numerator,
            "denominator": outer.denominator,
            "decimal": float(outer),
            "log2": log2f(outer),
        },
        "complete_cross_preflight_success": {
            "numerator": cross.numerator,
            "denominator": cross.denominator,
            "decimal": float(cross),
            "log2": log2f(cross),
            "penalty_bits": -log2f(cross),
        },
    }


def main() -> None:
    rows = [certify(x) for x in SHAPES]
    data = {
        "field": {"q": Q, "ell": D, "A_size": A_SIZE, "rho": "14/256"},
        "relation": "R_j=S^j",
        "offsets": "Gamma_0=1, Gamma_1=0; later Gamma_j arbitrary public choices",
        "coordinate_regions": {
            "source": "W x V'",
            "projection": "symmetric W x W and W x O",
            "outer": "independent public outer-map bytes",
            "internal_core": "V' x V' (not used in this certificate)",
        },
        "results": rows,
    }
    (ROOT / "structural_preflight_probability.json").write_text(json.dumps(data, indent=2) + "\n")

    md = [
        "# Structural public-preflight probability certificate",
        "",
        "This certificate replaces the earlier single-minor anti-concentration bound by direct kernel-union bounds.  The source, projection, outer-map, and internal-core byte regions are disjoint in the pinned expansion.",
        "",
        "For the affine source, the raw `(j,k)=(1,0)` labels quotient by the symmetric-square image to one `wedge^2(A)` per public base form.  Every nonzero dual alternating tuple imposes `4(v-1)` independent affine constraints on the fresh `W x V'` symbols.  Hence",
        "",
        "```text",
        "delta_src <= ((19^(6m1)-1)/18) * (14/256)^(4(v-1)).",
        "```",
        "",
        "For projection injectivity, the bad `A`-line `W` forces all `10m1` symmetric `W x W` symbols to vanish.  Any other `A`-line in `J=W+O` forces `16m1` independent `W x O` constraints.  Conditional on these source facts, one full-rank outer-map event implies both `rank Q=K` and `rank[Q Lambda]=M`.",
        "",
        "| Shape | affine rank proved / needed | log2 delta_src | log2 delta_proj | outer success | complete cross success | penalty bits |",
        "|:--:|:--:|--:|--:|--:|--:|--:|",
    ]
    for r in rows:
        md.append(
            f"| {r['key']} | {r['alternating_affine_rank_certified']} / {r['affine_rank_needed']} | "
            f"{r['source_failure']['log2']:.2f} | {r['projection_failure']['log2']:.2f} | "
            f"{r['outer_success']['decimal']:.9f} | {r['complete_cross_preflight_success']['decimal']:.9f} | "
            f"{r['complete_cross_preflight_success']['penalty_bits']:.4f} |"
        )
    md += [
        "",
        "All displayed probabilities are exact rationals in the JSON.  The decimal cross-preflight values differ only far beyond the displayed digits because the source and projection failures are already negligible; the visible `0.942322...` factor is the conservative biased-product full-rank bound for the outer map.",
    ]
    (ROOT / "STRUCTURAL_PREFLIGHT_PROBABILITY_CERTIFICATE.md").write_text("\n".join(md) + "\n")

    print("wrote structural_preflight_probability.json and STRUCTURAL_PREFLIGHT_PROBABILITY_CERTIFICATE.md")
    for r in rows:
        print(r["key"], "src", f"2^{r['source_failure']['log2']:.3f}", "proj", f"2^{r['projection_failure']['log2']:.3f}", "cross", f"{r['complete_cross_preflight_success']['decimal']:.12f}")


if __name__ == "__main__":
    main()
