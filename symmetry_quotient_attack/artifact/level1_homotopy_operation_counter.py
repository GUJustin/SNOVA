#!/usr/bin/env python3
"""Route-specific Level-I feasibility counter for the homotopy ledger.

The all-nine ledger charges the unit-leading term

    B * Bplus * H * N

from Proposition 5 of Safey El Din and Schost, plus the already explicit
separator and field-multiplication schedule.  It leaves the soft polynomial
arithmetic and the evaluation/interpolation transforms in ``kappa_hom``.

This script does not invent a value for that multiplier.  It exposes the exact
Level-I branch mixture and applies several explicit arithmetic profiles:

* an exact polynomial-basis realization of the small separator extension
  fields, using sparse irreducible moduli;
* the bilinear multiplication floor for the two polynomial products appearing
  in the cited one-parameter lifting bound;
* multiplication-only floors for padded Karatsuba and schoolbook polynomial
  multiplication.

All three polynomial rows are necessary-cost audits, not complete solver
bounds.  They omit polynomial additions, the outer soft constant, and the
later reconstruction, cleaning, and filtering work.  Therefore a row that
exceeds the NIST reference rules out that kernel, while a row that fits does
not prove the complete attack bound.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable

HERE = Path(__file__).resolve().parent
SOURCE_LEDGER = HERE / "all_nine_ledger.json"
OUTPUT_LEDGER = HERE / "level1_homotopy_operation_ledger.json"

REFERENCE = 143.0


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def logadd2(left: float, right: float) -> float:
    high, low = max(left, right), min(left, right)
    return high + math.log2(1.0 + 2.0 ** (low - high))


def adaptive_log2(
    fast_log2: float,
    fallback_log2: float,
    fallback_probability: float,
    common_failure: float,
    target_log2: float | None,
) -> float:
    if not 0.0 <= fallback_probability <= 1.0:
        raise ValueError("invalid fallback probability")
    if fallback_probability == 0.0:
        mixed = fast_log2
    elif fallback_probability == 1.0:
        mixed = fallback_log2
    else:
        mixed = logadd2(
            fast_log2 + math.log2(1.0 - fallback_probability),
            fallback_log2 + math.log2(fallback_probability),
        )
    normalized = mixed - math.log2(1.0 - common_failure)
    return normalized if target_log2 is None else logadd2(normalized, target_log2)


def additional_common_multiplier_bits(
    route: "Route", fast_log2: float, fallback_log2: float
) -> float:
    """Largest common extra factor on both solve branches below Level I."""

    low, high = 0.0, 256.0
    for _ in range(100):
        middle = (low + high) / 2
        total = adaptive_log2(
            fast_log2 + middle,
            fallback_log2 + middle,
            route.fallback_probability,
            route.common_failure,
            route.target_log2,
        )
        if total < REFERENCE:
            low = middle
        else:
            high = middle
    return low


def extension_multiplication_cost(
    extension_degree: int, base_multiplication_gates: int, base_addition_gates: int
) -> dict[str, int | float | str]:
    """Constructive polynomial-basis upper bound for one extension product."""

    points = 2 * extension_degree - 1
    multiplications = extension_degree * extension_degree
    convolution_additions = (extension_degree - 1) ** 2
    if base_multiplication_gates == 692 and extension_degree == 5:
        # E = F19^2[z]/(z^5+(2+u)).  Multiplication by -(2+u) maps
        # (a,b) to (b-2a,-a-2b).  Two doublings and the needed
        # additions/subtractions cost 200 scalar AXN gates.
        modulus = "z^5 + (2 + u)"
        reduction_additions = extension_degree - 1
        special_linear_gates = (extension_degree - 1) * 200
    elif base_multiplication_gates == 692 and extension_degree == 12:
        # E = F19^2[z]/(z^12+(1+u)).  Multiplication by -(1+u) maps
        # (a,b) to (b-a,-a-b), costing 116 scalar AXN gates.
        modulus = "z^12 + (1 + u)"
        reduction_additions = extension_degree - 1
        special_linear_gates = (extension_degree - 1) * 116
    elif base_multiplication_gates == 692 and extension_degree == 13:
        # E = F19^2[z]/(z^13+z+1).  Every high coefficient is negated once
        # and added into two low coefficients.
        modulus = "z^13 + z + 1"
        reduction_additions = 2 * (extension_degree - 1)
        special_linear_gates = (extension_degree - 1) * 32
    elif base_multiplication_gates == 2628 and extension_degree == 4:
        # E = F19^4[z]/(z^4+t).  In the paper basis t^4=t+1,
        # multiplication by -t costs four scalar negations and one scalar add.
        modulus = "z^4 + t"
        reduction_additions = extension_degree - 1
        special_linear_gates = (extension_degree - 1) * (4 * 16 + 42)
    else:
        raise ValueError(
            f"no sparse separator field for e={extension_degree}, "
            f"base multiplication={base_multiplication_gates}"
        )
    additions = convolution_additions + reduction_additions
    gates = (
        multiplications * base_multiplication_gates
        + additions * base_addition_gates
        + special_linear_gates
    )
    unit_ledger_gates = points * base_multiplication_gates
    return {
        "method": "schoolbook convolution with sparse irreducible reduction",
        "modulus": modulus,
        "extension_degree": extension_degree,
        "unit_ledger_pointwise_products": points,
        "base_multiplications": multiplications,
        "base_additions": additions,
        "special_linear_AXN_gates": special_linear_gates,
        "explicit_AXN_gates": gates,
        "unit_ledger_AXN_gates": unit_ledger_gates,
        "overhead_factor": gates / unit_ledger_gates,
        "overhead_bits": math.log2(gates / unit_ledger_gates),
    }


def schoolbook_multiplication_floor(length: int) -> int:
    """Coefficient multiplications used by the schoolbook kernel."""

    return length * length


def karatsuba_multiplication_floor(length: int) -> int:
    """Coefficient multiplications used by power-of-two padded Karatsuba."""

    padded = 1 << (length - 1).bit_length()
    multiplications = 1
    size = 1
    while size < padded:
        multiplications *= 3
        size *= 2
    return multiplications


def bilinear_floor_operations(length: int) -> int:
    """Multiplication-only rank floor for two length-``length`` polynomials."""

    return 2 * length - 1


@dataclass(frozen=True)
class Branch:
    name: str
    base_log2: float
    bezout: int
    companion_degree: int
    extension_degree: int
    base_multiplication_gates: int
    base_addition_gates: int


@dataclass(frozen=True)
class Route:
    name: str
    parameters: tuple[int, int, int, int]
    fast: Branch
    fallback: Branch
    fallback_probability: float
    common_failure: float
    target_log2: float | None


def l2_route(source: dict) -> Route:
    params = (48, 16, 2, 2)
    row = source["rows"][str(params)]
    m = params[1]
    fast_B = 2**m
    fast_Bplus = fast_B + m * 2 ** (m - 1)
    fallback_B = 2 ** (2 * m) * math.comb(m, m // 2)
    fallback_Bplus_fraction = Fraction(fallback_B) * (
        1 + 2 * m + Fraction(m * m, m + 2)
    )
    assert fallback_Bplus_fraction.denominator == 1
    return Route(
        name="ell2-level1-adaptive",
        parameters=params,
        fast=Branch(
            "one-eigenblock",
            row["diagonal"]["per_good_key_log2_AXN"],
            fast_B,
            fast_Bplus,
            row["diagonal"]["homotopy"]["extension_degree"],
            692,
            84,
        ),
        fallback=Branch(
            "complete-two-block",
            row["complete"]["per_good_key_log2_AXN"],
            fallback_B,
            fallback_Bplus_fraction.numerator,
            row["complete"]["homotopy"]["extension_degree"],
            692,
            84,
        ),
        fallback_probability=row["adaptive"]["fast_preflight_failure_upper"]
        / (1.0 - row["adaptive"]["common_spectrum_failure_upper"]),
        common_failure=row["adaptive"]["common_spectrum_failure_upper"],
        target_log2=row["diagonal"]["target_filter_log2_AXN"],
    )


def l4_route(source: dict, params: tuple[int, int, int, int]) -> Route:
    row = source["rows"][str(params)]
    fast_profile = row["fast"]["profile"]
    _s, channel_two, channel_twenty = fast_profile
    fast_B = 2**channel_two * 20**channel_twenty
    fast_Bplus_fraction = Fraction(fast_B) * (
        1 + Fraction(channel_two, 2) + Fraction(channel_twenty, 20)
    )
    assert fast_Bplus_fraction.denominator == 1
    K = row["complete_square"]["K"]
    fallback_B = 2**K
    fallback_Bplus = fallback_B * (K + 2) // 2
    return Route(
        name=f"ell4-level1-adaptive-{params}",
        parameters=params,
        fast=Branch(
            "five-quadratic-five-degree-twenty-channel",
            row["fast"]["per_good_key_log2_AXN"],
            fast_B,
            fast_Bplus_fraction.numerator,
            row["fast"]["homotopy"]["extension_degree"],
            2628,
            168,
        ),
        fallback=Branch(
            "complete-ordinary-square",
            row["complete_square"]["per_good_key_log2_AXN"],
            fallback_B,
            fallback_Bplus,
            row["complete_square"]["homotopy"]["extension_degree"],
            692,
            84,
        ),
        fallback_probability=row["adaptive"]["fast_preflight_failure_upper"]
        / (1.0 - row["adaptive"]["common_spectrum_failure_upper"]),
        common_failure=row["adaptive"]["common_spectrum_failure_upper"],
        target_log2=None,
    )


def polynomial_ratio(branch: Branch, operation_count: Callable[[int], int]) -> Fraction:
    first = operation_count(branch.bezout)
    second = operation_count(4 * branch.companion_degree)
    return Fraction(first * second, branch.bezout * branch.companion_degree)


def branch_factor(
    branch: Branch, operation_count: Callable[[int], int] | None
) -> tuple[Fraction, dict]:
    field = extension_multiplication_cost(
        branch.extension_degree,
        branch.base_multiplication_gates,
        branch.base_addition_gates,
    )
    field_ratio = Fraction(
        int(field["explicit_AXN_gates"]), int(field["unit_ledger_AXN_gates"])
    )
    if operation_count is None:
        poly_ratio = Fraction(1)
    else:
        poly_ratio = polynomial_ratio(branch, operation_count)
    return field_ratio * poly_ratio, {
        "B": branch.bezout,
        "Bplus": branch.companion_degree,
        "lifting_precision_4Bplus": 4 * branch.companion_degree,
        "extension_field_realization": field,
        "polynomial_operation_ratio": [poly_ratio.numerator, poly_ratio.denominator],
        "combined_factor": [
            (field_ratio * poly_ratio).numerator,
            (field_ratio * poly_ratio).denominator,
        ],
        "combined_overhead_bits": log2_fraction(field_ratio * poly_ratio),
    }


def count_route(route: Route) -> dict:
    profiles: list[tuple[str, Callable[[int], int] | None, str]] = [
        (
            "explicit-extension-field-only",
            None,
            "Constructive field realization; polynomial soft factors remain unit-leading.",
        ),
        (
            "bilinear-polynomial-floor",
            bilinear_floor_operations,
            "Feasibility lower bound; it omits all polynomial additions and all outer lifting work.",
        ),
        (
            "padded-karatsuba-multiplication-floor",
            karatsuba_multiplication_floor,
            "Necessary multiplication count; omits additions, outer lifting, and later solver work.",
        ),
        (
            "schoolbook-multiplication-floor",
            schoolbook_multiplication_floor,
            "Necessary multiplication count; omits additions, outer lifting, and later solver work.",
        ),
    ]
    base = adaptive_log2(
        route.fast.base_log2,
        route.fallback.base_log2,
        route.fallback_probability,
        route.common_failure,
        route.target_log2,
    )
    results = {}
    for name, operation_count, meaning in profiles:
        fast_factor, fast_detail = branch_factor(route.fast, operation_count)
        fallback_factor, fallback_detail = branch_factor(route.fallback, operation_count)
        total = adaptive_log2(
            route.fast.base_log2 + log2_fraction(fast_factor),
            route.fallback.base_log2 + log2_fraction(fallback_factor),
            route.fallback_probability,
            route.common_failure,
            route.target_log2,
        )
        extra_bits = additional_common_multiplier_bits(
            route,
            route.fast.base_log2 + log2_fraction(fast_factor),
            route.fallback.base_log2 + log2_fraction(fallback_factor),
        )
        results[name] = {
            "meaning": meaning,
            "fast": fast_detail,
            "fallback": fallback_detail,
            "adaptive_log2_AXN": total,
            "remaining_bits_below_143": REFERENCE - total,
            "maximum_additional_common_multiplier_bits": extra_bits,
            "maximum_additional_common_multiplier": 2.0**extra_bits,
            "fits_143": total < REFERENCE,
        }
    return {
        "parameters": list(route.parameters),
        "fast_branch": route.fast.name,
        "fallback_branch": route.fallback.name,
        "fallback_probability": route.fallback_probability,
        "common_failure": route.common_failure,
        "target_log2_AXN": route.target_log2,
        "unit_leading_adaptive_log2_AXN": base,
        "unit_leading_headroom_bits": REFERENCE - base,
        "profiles": results,
    }


def main() -> None:
    source = json.loads(SOURCE_LEDGER.read_text())
    routes = [
        l2_route(source),
        l4_route(source, (28, 5, 4, 4)),
        l4_route(source, (28, 4, 4, 5)),
    ]
    route_results = {route.name: count_route(route) for route in routes}
    profile_names = next(iter(route_results.values()))["profiles"].keys()
    worst = {}
    for profile in profile_names:
        route_name, result = max(
            (
                (name, route["profiles"][profile])
                for name, route in route_results.items()
            ),
            key=lambda item: item[1]["adaptive_log2_AXN"],
        )
        worst[profile] = {
            "route": route_name,
            "adaptive_log2_AXN": result["adaptive_log2_AXN"],
            "remaining_bits_below_143": result["remaining_bits_below_143"],
            "fits_143": result["fits_143"],
        }

    report = {
        "scope": (
            "Level-I route-specific necessary-cost audit. This is not a complete "
            "homotopy solver bound and does not assign a value to kappa_hom."
        ),
        "source_claims": {
            "unit_leading_ledger": "all_nine_ledger.json",
            "homotopy_bound": "Safey-El-Din--Schost Proposition 5",
            "polynomial_kernel_shape": (
                "Schost Theorem 2, one parameter: multiplication at B and at "
                "precision bounded by 4 Bplus"
            ),
        },
        "reference_log2_AXN": REFERENCE,
        "routes": route_results,
        "worst_level1_by_profile": worst,
        "interpretation": {
            "exceeds": (
                "A profile above 143 is ruled out even before the omitted "
                "lifting, reconstruction, cleaning, and filtering costs."
            ),
            "fits": (
                "A profile below 143 only survives this necessary-cost audit; "
                "it is not yet a proved attack bound."
            ),
        },
    }
    OUTPUT_LEDGER.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("SNOVA Level-I homotopy necessary-cost audit")
    for name, result in worst.items():
        print(
            f"{name}: {result['adaptive_log2_AXN']:.6f} bits; "
            f"remaining {result['remaining_bits_below_143']:.6f}; "
            f"fits={result['fits_143']}"
        )
    print(f"wrote {OUTPUT_LEDGER}")


if __name__ == "__main__":
    main()
