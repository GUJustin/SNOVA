#!/usr/bin/env python3
"""Independent checks for the Level-I homotopy necessary-cost audit."""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def log2_fraction(value: Fraction) -> float:
    return math.log2(value.numerator) - math.log2(value.denominator)


def logadd(left: float, right: float) -> float:
    high, low = max(left, right), min(left, right)
    return high + math.log2(1 + 2 ** (low - high))


def mix(
    fast: float,
    fallback: float,
    probability: float,
    failure: float,
    target: float | None,
) -> float:
    value = logadd(
        fast + math.log2(1 - probability),
        fallback + math.log2(probability),
    ) - math.log2(1 - failure)
    return value if target is None else logadd(value, target)


def field_ratio(e: int, mul: int, add: int) -> Fraction:
    points = 2 * e - 1
    multiplications = e * e
    convolution_additions = (e - 1) ** 2
    if mul == 692 and e == 5:
        reduction_additions = e - 1
        special = (e - 1) * 200
    elif mul == 692 and e == 12:
        reduction_additions = e - 1
        special = (e - 1) * 116
    elif mul == 692 and e == 13:
        reduction_additions = 2 * (e - 1)
        special = (e - 1) * 32
    elif mul == 2628 and e == 4:
        reduction_additions = e - 1
        special = (e - 1) * (4 * 16 + 42)
    else:
        raise AssertionError((e, mul, add))
    explicit = (
        multiplications * mul
        + (convolution_additions + reduction_additions) * add
        + special
    )
    return Fraction(explicit, points * mul)


def a2_add(x, y):
    return ((x[0] + y[0]) % 19, (x[1] + y[1]) % 19)


def a2_sub(x, y):
    return ((x[0] - y[0]) % 19, (x[1] - y[1]) % 19)


def a2_mul(x, y):
    return (
        (x[0] * y[0] - x[1] * y[1]) % 19,
        (x[0] * y[1] + x[1] * y[0]) % 19,
    )


def a4_add(x, y):
    return tuple((a + b) % 19 for a, b in zip(x, y))


def a4_sub(x, y):
    return tuple((a - b) % 19 for a, b in zip(x, y))


def a4_mul(x, y):
    raw = [0] * 7
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            raw[i + j] = (raw[i + j] + a * b) % 19
    for degree in range(6, 3, -1):
        coefficient = raw[degree]
        raw[degree] = 0
        raw[degree - 3] = (raw[degree - 3] + coefficient) % 19
        raw[degree - 4] = (raw[degree - 4] + coefficient) % 19
    return tuple(raw[:4])


def field_pow(value, exponent, mul, one):
    result = one
    while exponent:
        if exponent & 1:
            result = mul(result, value)
        value = mul(value, value)
        exponent //= 2
    return result


def poly_trim(poly, zero):
    result = list(poly)
    while len(result) > 1 and result[-1] == zero:
        result.pop()
    return result


def poly_sub(left, right, sub, zero):
    size = max(len(left), len(right))
    result = []
    for i in range(size):
        a = left[i] if i < len(left) else zero
        b = right[i] if i < len(right) else zero
        result.append(sub(a, b))
    return poly_trim(result, zero)


def poly_mul(left, right, add, mul, zero):
    result = [zero] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = add(result[i + j], mul(a, b))
    return poly_trim(result, zero)


def poly_divmod(numerator, denominator, sub, mul, inv, zero):
    remainder = poly_trim(numerator, zero)
    denominator = poly_trim(denominator, zero)
    if denominator == [zero]:
        raise ZeroDivisionError
    quotient = [zero] * max(1, len(remainder) - len(denominator) + 1)
    lead_inverse = inv(denominator[-1])
    while len(remainder) >= len(denominator) and remainder != [zero]:
        shift = len(remainder) - len(denominator)
        coefficient = mul(remainder[-1], lead_inverse)
        quotient[shift] = coefficient
        subtractor = [zero] * shift + [mul(coefficient, x) for x in denominator]
        remainder = poly_sub(remainder, subtractor, sub, zero)
    return poly_trim(quotient, zero), remainder


def poly_mod_pow(base, exponent, modulus, add, sub, mul, inv, zero, one):
    result = [one]
    base = poly_divmod(base, modulus, sub, mul, inv, zero)[1]
    while exponent:
        if exponent & 1:
            result = poly_divmod(
                poly_mul(result, base, add, mul, zero),
                modulus,
                sub,
                mul,
                inv,
                zero,
            )[1]
        base = poly_divmod(
            poly_mul(base, base, add, mul, zero),
            modulus,
            sub,
            mul,
            inv,
            zero,
        )[1]
        exponent //= 2
    return result


def poly_gcd(left, right, sub, mul, inv, zero):
    left, right = poly_trim(left, zero), poly_trim(right, zero)
    while right != [zero]:
        left, right = right, poly_divmod(left, right, sub, mul, inv, zero)[1]
    scale = inv(left[-1])
    return [mul(scale, coefficient) for coefficient in left]


def prime_divisors(value: int) -> set[int]:
    result = set()
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result.add(divisor)
            value //= divisor
        divisor += 1
    if value > 1:
        result.add(value)
    return result


def assert_irreducible(modulus, q, add, sub, mul, zero, one):
    degree = len(modulus) - 1
    inv = lambda x: field_pow(x, q - 2, mul, one)
    variable = [zero, one]
    assert poly_sub(
        poly_mod_pow(variable, q**degree, modulus, add, sub, mul, inv, zero, one),
        variable,
        sub,
        zero,
    ) == [zero]
    for prime in prime_divisors(degree):
        witness = poly_sub(
            poly_mod_pow(
                variable,
                q ** (degree // prime),
                modulus,
                add,
                sub,
                mul,
                inv,
                zero,
                one,
            ),
            variable,
            sub,
            zero,
        )
        assert len(poly_gcd(modulus, witness, sub, mul, inv, zero)) == 1


def bilinear(length: int) -> int:
    return 2 * length - 1


def schoolbook_multiplications(length: int) -> int:
    return length * length


def karatsuba_multiplications(length: int) -> int:
    padded = 1 << (length - 1).bit_length()
    muls, size = 1, 1
    while size < padded:
        size *= 2
        muls *= 3
    return muls


def polynomial_ratio(B: int, Bplus: int, operation_count) -> Fraction:
    return Fraction(
        operation_count(B) * operation_count(4 * Bplus), B * Bplus
    )


def check_close(actual: float, expected: float) -> None:
    assert abs(actual - expected) < 1e-9, (actual, expected)


report = json.loads((HERE / "level1_homotopy_operation_ledger.json").read_text())
source = json.loads((HERE / "all_nine_ledger.json").read_text())
assert report["reference_log2_AXN"] == 143.0
assert len(report["routes"]) == 3

# Rabin irreducibility checks for every sparse separator-field modulus.
a2_zero, a2_one = (0, 0), (1, 0)
for degree, constant in ((5, (2, 1)), (12, (1, 1))):
    assert_irreducible(
        [constant] + [a2_zero] * (degree - 1) + [a2_one],
        19**2,
        a2_add,
        a2_sub,
        a2_mul,
        a2_zero,
        a2_one,
    )
assert_irreducible(
    [a2_one, a2_one] + [a2_zero] * 11 + [a2_one],
    19**2,
    a2_add,
    a2_sub,
    a2_mul,
    a2_zero,
    a2_one,
)
a4_zero, a4_one, a4_t = (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0)
assert_irreducible(
    [a4_t, a4_zero, a4_zero, a4_zero, a4_one],
    19**4,
    a4_add,
    a4_sub,
    a4_mul,
    a4_zero,
    a4_one,
)

# Reconstruct the ell=2 route without importing the generator or counter.
params = (48, 16, 2, 2)
source_row = source["rows"][str(params)]
route = report["routes"]["ell2-level1-adaptive"]
m = 16
B_fast = 2**m
Bplus_fast = B_fast + m * 2 ** (m - 1)
B_fallback = 2 ** (2 * m) * math.comb(m, m // 2)
Bplus_fallback_fraction = Fraction(B_fallback) * (
    1 + 2 * m + Fraction(m * m, m + 2)
)
assert Bplus_fallback_fraction.denominator == 1
Bplus_fallback = Bplus_fallback_fraction.numerator
probability = source_row["adaptive"]["fast_preflight_failure_upper"] / (
    1 - source_row["adaptive"]["common_spectrum_failure_upper"]
)
failure = source_row["adaptive"]["common_spectrum_failure_upper"]
target = source_row["diagonal"]["target_filter_log2_AXN"]
fast_base = source_row["diagonal"]["per_good_key_log2_AXN"]
fallback_base = source_row["complete"]["per_good_key_log2_AXN"]
check_close(
    route["unit_leading_adaptive_log2_AXN"],
    mix(fast_base, fallback_base, probability, failure, target),
)

models = {
    "bilinear-polynomial-floor": bilinear,
    "padded-karatsuba-multiplication-floor": karatsuba_multiplications,
    "schoolbook-multiplication-floor": schoolbook_multiplications,
}
fast_field = field_ratio(5, 692, 84)
fallback_field = field_ratio(12, 692, 84)

field_only = route["profiles"]["explicit-extension-field-only"]
field_total = mix(
    fast_base + log2_fraction(fast_field),
    fallback_base + log2_fraction(fallback_field),
    probability,
    failure,
    target,
)
check_close(field_only["adaptive_log2_AXN"], field_total)

for name, operation_count in models.items():
    fast_factor = fast_field * polynomial_ratio(
        B_fast, Bplus_fast, operation_count
    )
    fallback_factor = fallback_field * polynomial_ratio(
        B_fallback, Bplus_fallback, operation_count
    )
    expected = mix(
        fast_base + log2_fraction(fast_factor),
        fallback_base + log2_fraction(fallback_factor),
        probability,
        failure,
        target,
    )
    check_close(route["profiles"][name]["adaptive_log2_AXN"], expected)

# Check the two ell=4 branch shapes and exact Level-I parameters.
for params in ((28, 5, 4, 4), (28, 4, 4, 5)):
    name = f"ell4-level1-adaptive-{params}"
    counted = report["routes"][name]
    source_row = source["rows"][str(params)]
    check_close(
        counted["unit_leading_adaptive_log2_AXN"],
        source_row["adaptive"]["normalized_log2_AXN"],
    )
    fast = counted["profiles"]["explicit-extension-field-only"]["fast"]
    fallback = counted["profiles"]["explicit-extension-field-only"]["fallback"]
    assert fast["B"] == 2**5 * 20**5
    assert fast["Bplus"] == (2**5 * 20**5) * 15 // 4
    assert fallback["B"] == 2**50
    assert fallback["Bplus"] == 26 * 2**50
    assert fast["extension_field_realization"]["modulus"] == "z^4 + t"
    assert fallback["extension_field_realization"]["modulus"] == "z^13 + z + 1"

worst = report["worst_level1_by_profile"]
check_close(worst["explicit-extension-field-only"]["adaptive_log2_AXN"], 133.86429722728064)
check_close(worst["bilinear-polynomial-floor"]["adaptive_log2_AXN"], 137.86428604306457)
assert worst["bilinear-polynomial-floor"]["remaining_bits_below_143"] < 5.136
ell2_floor = report["routes"]["ell2-level1-adaptive"]["profiles"][
    "bilinear-polynomial-floor"
]
assert ell2_floor["maximum_additional_common_multiplier"] < 35.157
assert worst["explicit-extension-field-only"]["fits_143"]
assert worst["bilinear-polynomial-floor"]["fits_143"]
check_close(
    worst["padded-karatsuba-multiplication-floor"]["adaptive_log2_AXN"],
    195.64271074890303,
)
check_close(
    worst["schoolbook-multiplication-floor"]["adaptive_log2_AXN"],
    239.4526026074797,
)
assert not worst["padded-karatsuba-multiplication-floor"]["fits_143"]
assert not worst["schoolbook-multiplication-floor"]["fits_143"]

print("Level-I homotopy operation-counter checks passed")
print("- exact adaptive branch mixtures reproduced for all three Level-I rows")
print("- all four sparse separator-field moduli pass Rabin irreducibility checks")
print("- sparse polynomial-basis extension-field counts reproduced")
print("- bilinear floor leaves less than 5.136 bits on the ell=2 route")
print("- all remaining solve work must fit a common factor below 35.157")
print("- multiplication-only padded Karatsuba and schoolbook floors exceed 143")
