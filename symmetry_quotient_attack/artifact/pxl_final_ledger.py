#!/usr/bin/env python3
"""Reproduce the final PXL route ledgers used by the revision.

The script deliberately separates three kinds of quantities:

* exact arithmetic: dimensions, Hilbert coefficients, finite-field
  probabilities, integer operation counts, target work, and storage;
* route premises: the affine-rank and coefficient-support assumptions needed
  before the PXL model applies; and
* estimates: the PXL/Macaulay model and the conversion of one F_19 operation
  to 150 binary gates.

Run with ``--write`` to regenerate ``pxl_final_ledger.json`` and with
``--check`` (the default) to compare a fresh computation with the checked-in
ledger.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable


HERE = Path(__file__).resolve().parent
LEDGER = HERE / "pxl_final_ledger.json"

Q = 19
GATES_PER_F19_OP = 150
BITS_PER_F19_ELEMENT = (Q - 1).bit_length()  # ceil(log2(Q))
TARGET_SEARCH_BITS = 32
OMEGAS = (Fraction(281, 100), Fraction(3, 1))
REFERENCE_EXPONENTS = {"I": 143, "III": 207, "V": 272}


ELL2_SHAPES = (
    {"level": "I", "v": 48, "o": 16, "ell": 2, "r": 2, "m1": 16, "M": 64},
    {"level": "III", "v": 72, "o": 24, "ell": 2, "r": 2, "m1": 24, "M": 96},
    {"level": "V", "v": 96, "o": 32, "ell": 2, "r": 2, "m1": 32, "M": 128},
)

ELL4_SHAPES = (
    {"level": "I", "v": 28, "o": 5, "ell": 4, "r": 4, "m1": 5, "M": 80},
    {"level": "I", "v": 28, "o": 4, "ell": 4, "r": 5, "m1": 5, "M": 80},
    {"level": "III", "v": 40, "o": 7, "ell": 4, "r": 4, "m1": 7, "M": 112},
    {"level": "III", "v": 38, "o": 5, "ell": 4, "r": 5, "m1": 7, "M": 100},
    {"level": "V", "v": 50, "o": 9, "ell": 4, "r": 4, "m1": 9, "M": 144},
    {"level": "V", "v": 52, "o": 6, "ell": 4, "r": 6, "m1": 9, "M": 144},
)


# These are regression checks, not inputs to the optimization.  Each number
# must be reproduced independently from the formulas below.
EXPECTED_ROUNDED = {
    ("ell2_supported", "I", 48, 16, "2.81"): (9, 12, 128.833447),
    ("ell2_supported", "III", 72, 24, "2.81"): (13, 17, 187.029907),
    ("ell2_supported", "V", 96, 32, "2.81"): (15, 23, 245.041754),
    ("ell2_supported", "I", 48, 16, "3"): (9, 12, 133.667445),
    ("ell2_supported", "III", 72, 24, "3"): (14, 16, 194.626213),
    ("ell2_supported", "V", 96, 32, "3"): (19, 20, 255.549586),
    ("ell4_supported", "I", 28, 5, "2.81"): (7, 7, 145.867752),
    ("ell4_supported", "I", 28, 4, "2.81"): (7, 7, 145.867386),
    ("ell4_supported", "III", 40, 7, "2.81"): (8, 9, 197.338396),
    ("ell4_supported", "III", 38, 5, "2.81"): (8, 8, 201.265361),
    ("ell4_supported", "V", 50, 9, "2.81"): (9, 10, 252.202579),
    ("ell4_supported", "V", 52, 6, "2.81"): (9, 10, 249.538960),
    ("ell4_supported", "I", 28, 5, "3"): (7, 7, 148.400390),
    ("ell4_supported", "I", 28, 4, "3"): (7, 7, 148.400025),
    ("ell4_supported", "III", 40, 7, "3"): (9, 9, 201.379285),
    ("ell4_supported", "III", 38, 5, "3"): (8, 8, 205.022213),
    ("ell4_supported", "V", 50, 9, "3"): (10, 10, 257.208548),
    ("ell4_supported", "V", 52, 6, "3"): (9, 10, 254.991164),
    ("ell4_H_struct", "I", 28, 5, "2.81"): (10, 12, 134.037961),
    ("ell4_H_struct", "I", 28, 4, "2.81"): (10, 12, 134.037595),
    ("ell4_H_struct", "III", 40, 7, "2.81"): (12, 17, 182.006674),
    ("ell4_H_struct", "III", 38, 5, "2.81"): (12, 17, 182.006152),
    ("ell4_H_struct", "V", 50, 9, "2.81"): (15, 21, 230.444951),
    ("ell4_H_struct", "V", 52, 6, "2.81"): (15, 21, 230.444296),
    ("ell4_H_struct", "I", 28, 5, "3"): (10, 12, 138.911992),
    ("ell4_H_struct", "I", 28, 4, "3"): (10, 12, 138.911626),
    ("ell4_H_struct", "III", 40, 7, "3"): (13, 16, 189.638021),
    ("ell4_H_struct", "III", 38, 5, "3"): (13, 16, 189.637499),
    ("ell4_H_struct", "V", 50, 9, "3"): (18, 19, 240.174424),
    ("ell4_H_struct", "V", 52, 6, "3"): (18, 19, 240.173769),
}


@dataclass(frozen=True)
class PXLParameters:
    K: int
    h: int
    k: int
    D: int
    A: int
    B: int


def omega_label(omega: Fraction) -> str:
    if omega.denominator == 1:
        return str(omega.numerator)
    return f"{float(omega):.2f}"


def log2_int(value: int) -> float:
    """Accurate log2 for arbitrarily large positive integers."""
    if value <= 0:
        raise ValueError("log2 requires a positive integer")
    shift = max(value.bit_length() - 53, 0)
    return math.log2(value >> shift) + shift


def log2_fraction(value: Fraction) -> float:
    if value <= 0:
        raise ValueError("log2 requires a positive fraction")
    return log2_int(value.numerator) - log2_int(value.denominator)


def log2_sum(log_values: Iterable[float]) -> float:
    values = tuple(log_values)
    top = max(values)
    return top + math.log2(sum(2.0 ** (value - top) for value in values))


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
    }


def polynomial_coefficients(K: int, exponent_minus: int) -> list[int]:
    """Coefficients of (1-t)^exponent_minus (1+t)^K."""
    left = [((-1) ** i) * math.comb(exponent_minus, i) for i in range(exponent_minus + 1)]
    right = [math.comb(K, i) for i in range(K + 1)]
    result = [0] * (len(left) + len(right) - 1)
    for i, lhs in enumerate(left):
        for j, rhs in enumerate(right):
            result[i + j] += lhs * rhs
    return result


@lru_cache(maxsize=None)
def degree_and_A(K: int, h: int, k: int) -> tuple[int, int]:
    degree_coefficients = polynomial_coefficients(K, K - h + k - 1)
    D = next(
        degree
        for degree in range(2, len(degree_coefficients))
        if degree_coefficients[degree] <= 1
    )
    hilbert_coefficients = polynomial_coefficients(K, K - h + k)
    A = sum(max(hilbert_coefficients[degree], 0) for degree in range(D + 1))
    return D, A


@lru_cache(maxsize=None)
def pxl_parameters(K: int, h: int, k: int) -> PXLParameters:
    D, A = degree_and_A(K, h, k)
    B = math.comb(h - k + D, D)
    return PXLParameters(K=K, h=h, k=k, D=D, A=A, B=B)


def component_ledger(parameters: PXLParameters, omega: Fraction) -> dict[str, Any]:
    K, h, k, D, A, B = (
        parameters.K,
        parameters.h,
        parameters.k,
        parameters.D,
        parameters.A,
        parameters.B,
    )
    exact_components = {
        "C2": k * k * (h - k) * B * B,
        "C3": k * k * A * B * math.comb(h + D, D),
        "Cfix": (Q**k) * A * A * math.comb(k + D, D),
    }
    logs = {
        "C1": float(omega) * log2_int(B),
        "C2": log2_int(exact_components["C2"]),
        "C3": log2_int(exact_components["C3"]),
        "Cfix": log2_int(exact_components["Cfix"]),
        "Clin": k * math.log2(Q) + float(omega) * log2_int(A),
    }

    exact_integer: dict[str, str | None] = {
        "C1": None,
        "C2": str(exact_components["C2"]),
        "C3": str(exact_components["C3"]),
        "Cfix": str(exact_components["Cfix"]),
        "Clin": None,
    }
    total_exact: int | None = None
    if omega.denominator == 1:
        power = omega.numerator
        c1 = B**power
        clin = (Q**k) * (A**power)
        exact_integer["C1"] = str(c1)
        exact_integer["Clin"] = str(clin)
        total_exact = c1 + sum(exact_components.values()) + clin

    total_log = log2_sum(logs.values())
    if total_exact is not None:
        # The log-sum expression and the exact integer sum must agree.
        assert abs(total_log - log2_int(total_exact)) < 1e-10

    return {
        "formula": {
            "C1": "B^omega",
            "C2": "k^2*(h-k)*B^2",
            "C3": "k^2*A*B*binom(h+D,D)",
            "Cfix": "19^k*A^2*binom(k+D,D)",
            "Clin": "19^k*A^omega",
        },
        "component_log2_F19_ops": logs,
        "component_exact_integer_F19_ops": exact_integer,
        "C_PXL_log2_F19_ops": total_log,
        "C_PXL_exact_integer_F19_ops": str(total_exact) if total_exact is not None else None,
        "arithmetic_status": (
            "exact integer five-component sum"
            if total_exact is not None
            else "real-exponent operation-count estimate; C2, C3, and Cfix are exact integers"
        ),
    }


def ell2_probability(v: int, o: int, h: int, K: int) -> dict[str, Any]:
    n = v + o
    alpha_numerator = sum(
        math.comb(n, j) * (18 ** (n - j)) for j in range(n // 4 + 1)
    )
    alpha = Fraction(alpha_numerator, Q**n)
    eta = Fraction(Q**h - 1, Q**K) + Fraction(1, 2**128)
    accepted = alpha - eta / 18
    regular = Fraction(Q**h, Q**K) * accepted * accepted / (accepted + eta)
    # Exact union-bound subtraction for non-transverse directions.
    direct = regular - Fraction(1, Q**4)
    if direct <= 0:
        direct = Fraction(0, 1)
    return {
        "alpha": alpha,
        "eta": eta,
        "accepted_lower_bound": accepted,
        "p_regular_lower_bound": regular,
        "p_trial_lower_bound": direct,
        "transversality_subtraction": Fraction(1, Q**4),
    }


def ell4_probability(v: int, o: int, ell: int, r: int, h: int, K: int) -> dict[str, Any]:
    if r == ell:
        alpha = Fraction((Q**ell - 1) ** (v + o), Q ** (ell * (v + o)))
        alpha_source = "exact product probability (1-19^-4)^(v+o)"
    else:
        alpha = Fraction(1, 1)
        alpha_source = "rectangular-signature branch; alpha=1"
    eta = Fraction(Q**h - 1, Q**K) + Fraction(1, 2**128)
    accepted = alpha - eta / 18
    regular = Fraction(Q**h, Q**K) * accepted * accepted / (accepted + eta)
    return {
        "alpha": alpha,
        "alpha_source": alpha_source,
        "eta": eta,
        "accepted_lower_bound": accepted,
        "p_regular_lower_bound": regular,
        "p_trial_lower_bound": regular,
        "transversality_subtraction": Fraction(0, 1),
    }


def serialize_probability(probability: dict[str, Any]) -> dict[str, Any]:
    p_trial = probability["p_trial_lower_bound"]
    output: dict[str, Any] = {}
    for key, value in probability.items():
        output[key] = fraction_record(value) if isinstance(value, Fraction) else value
    output["p_trial_log2"] = log2_fraction(p_trial)
    output["retry_log2_upper_bound"] = -log2_fraction(p_trial)
    output["retry_multiplier_upper_bound"] = fraction_record(1 / p_trial)
    output["status"] = (
        "exact rational lower bound on trial probability and exact rational upper bound on retries"
    )
    return output


def work_ledger(
    pxl: dict[str, Any],
    probability: Fraction,
    target_codimension: int,
    level: str,
    target_charge_basis: str,
) -> dict[str, Any]:
    target_work = (2**TARGET_SEARCH_BITS) * (Q**target_codimension)
    pxl_gate_log = math.log2(GATES_PER_F19_OP) + pxl["C_PXL_log2_F19_ops"]
    target_log = log2_int(target_work)
    numerator_log = log2_sum((pxl_gate_log, target_log))
    total_log = numerator_log - log2_fraction(probability)
    exact_work = None
    if pxl["C_PXL_exact_integer_F19_ops"] is not None:
        exact_numerator = (
            GATES_PER_F19_OP * int(pxl["C_PXL_exact_integer_F19_ops"])
            + target_work
        )
        exact_work = Fraction(exact_numerator, 1) / probability
        assert abs(total_log - log2_fraction(exact_work)) < 1e-10
    reference = REFERENCE_EXPONENTS[level]
    return {
        "normalization_binary_gates_per_F19_operation": GATES_PER_F19_OP,
        "solver_binary_gate_log2": pxl_gate_log,
        "target_consistency_codimension": target_codimension,
        "target_work_arithmetic_value": str(target_work),
        "target_work_log2_binary_gates": target_log,
        "target_charge_basis": target_charge_basis,
        "target_charge_status": (
            "stipulated unit-leading model convention: 2^32 binary-gate units per retained "
            "target, times 19^target_consistency_codimension; the displayed integer is an "
            "exact evaluation of that convention, not a measured or circuit-level bound"
        ),
        "reported_work_formula": "(150*C_PXL + C_target) / p_trial",
        "reported_work_log2_binary_gates": total_log,
        "reported_work_exact_model_arithmetic": (
            fraction_record(exact_work) if exact_work is not None else None
        ),
        "reference_exponent": reference,
        "headroom_bits_before_C_aux": reference - total_log,
        "augmented_work_formula": "(150*C_PXL + C_target + C_aux) / p_trial",
        "C_aux_status": "symbolic and unpriced; it is not included in reported_work_log2_binary_gates",
        "status": "gate-count estimate conditional on the route premises and PXL cost model",
    }


def storage_ledger(parameters: PXLParameters) -> dict[str, Any]:
    coefficient_count = (
        parameters.B
        * parameters.B
        * math.comb(parameters.k + parameters.D, parameters.D)
    )
    bits = BITS_PER_F19_ELEMENT * coefficient_count
    return {
        "formula": "ceil(log2(19))*B^2*binom(k+D,D) bits",
        "bits_per_F19_element": BITS_PER_F19_ELEMENT,
        "exact_F19_coefficient_slots": str(coefficient_count),
        "exact_bits": str(bits),
        "log2_bits": log2_int(bits),
        "status": (
            "exact arithmetic value of the stated five-bit fixed-width dense-array proxy; "
            "not an information-theoretic packing minimum, total state, traffic count, or "
            "measured peak-memory bound"
        ),
    }


def route_entry(
    *,
    family: str,
    shape: dict[str, int | str],
    omega: Fraction,
    parameters: PXLParameters,
    probability: dict[str, Any],
    target_codimension: int,
    target_charge_basis: str,
    premises: list[str],
    hypothesis_metadata: dict[str, Any],
) -> dict[str, Any]:
    pxl = component_ledger(parameters, omega)
    trial_probability = probability["p_trial_lower_bound"]
    entry = {
        "route_family": family,
        "shape": shape,
        "omega": omega_label(omega),
        "PXL_parameters": {
            "K": parameters.K,
            "h": parameters.h,
            "k": parameters.k,
            "main_variables_h_minus_k": parameters.h - parameters.k,
            "D": parameters.D,
            "A": str(parameters.A),
            "B": str(parameters.B),
        },
        "route_premises": premises,
        "hypothesis_metadata": hypothesis_metadata,
        "probability": serialize_probability(probability),
        "PXL_components": pxl,
        "work": work_ledger(
            pxl,
            trial_probability,
            target_codimension,
            str(shape["level"]),
            target_charge_basis,
        ),
        "storage": storage_ledger(parameters),
    }
    return entry


def best_entry(
    candidates: Iterable[tuple[PXLParameters, dict[str, Any], int]],
    *,
    family: str,
    shape: dict[str, int | str],
    omega: Fraction,
    premises: list[str],
    hypothesis_metadata: dict[str, Any],
    target_charge_basis: str,
) -> dict[str, Any]:
    entries = [
        route_entry(
            family=family,
            shape=shape,
            omega=omega,
            parameters=parameters,
            probability=probability,
            target_codimension=target_codimension,
            target_charge_basis=target_charge_basis,
            premises=premises,
            hypothesis_metadata=hypothesis_metadata,
        )
        for parameters, probability, target_codimension in candidates
    ]
    return min(entries, key=lambda entry: entry["work"]["reported_work_log2_binary_gates"])


def ell4_public_premise_metadata(
    shape: dict[str, int | str], K: int, h: int
) -> dict[str, Any]:
    """Record the exact rank gates and the separate spectrum premise."""
    M = int(shape["M"])
    return {
        "exact_public_preflights": {
            "quotient_rank": {
                "claim": "rank(Q_R)=K",
                "value": K,
                "status": "exact public preflight",
            },
            "affine_source_rank": {
                "claim": "rank(C_R,v)=M-K",
                "value": M - K,
                "status": "exact public preflight",
            },
            "full_affine_rank": {
                "claim": "rank([Q_R|C_R,v])=M",
                "value": M,
                "status": "exact consequence of the two passing public rank preflights",
            },
        },
        "fixed_key_spectrum_premise": {
            "claim": "bar_eta_L^U <= eta_128(h,K)",
            "h": h,
            "K": K,
            "threshold_formula": "(19^h-1)*19^(-K)+2^(-128)",
            "status": (
                "fixed-key premise or ideal-transcript event; not an exact public "
                "preflight unless certified by an exact projective rank histogram"
            ),
        },
    }


def ell2_entries() -> list[dict[str, Any]]:
    entries = []
    premises = [
        "ordinary coefficient-supported slice from the exact affine reduction",
        "simultaneous graph quotient injectivity and the stated spectrum event",
        "direct PXL/Macaulay cost model with the complete five components",
    ]
    for shape in ELL2_SHAPES:
        v, o = int(shape["v"]), int(shape["o"])
        K = 3 * int(shape["m1"])
        h = v - 2
        probability = ell2_probability(v, o, h, K)
        target_codimension = int(shape["M"]) - K
        for omega in OMEGAS:
            candidates = (
                (pxl_parameters(K, h, k), probability, target_codimension)
                for k in range(1, h)
            )
            entries.append(
                best_entry(
                    candidates,
                    family="ell2_supported",
                    shape=shape,
                    omega=omega,
                    premises=premises,
                    hypothesis_metadata={
                        "ordinary_affine_support": "exact route theorem",
                        "H_PXL": "ordinary MQ degree, Hilbert, and Macaulay-rank premise",
                        "H_route": "instancewise solver-rank and completeness premise",
                    },
                    target_charge_basis=(
                        "stipulated 2^32 charge per retained target multiplied by the exact "
                        "zero-offset target-consistency factor 19^(M-K)"
                    ),
                )
            )
    return entries


def ell4_entries() -> list[dict[str, Any]]:
    entries = []
    supported_premises = [
        "exact public preflights rank(Q_R)=K and rank(C_R,v)=M-K, giving full affine rank M",
        "H_off: an oil-supported offset has lambda=M-K and satisfies rank(C|U_vin)=lambda",
        "H_off: A*ker(C|U_vin)=U_vin, so the affine kernel contains a coefficient-supporting main graph",
        "fixed-map homogeneous coefficient support holds on that main graph with h-k <= v",
        "homogeneous support is not a claim of full affine coefficient support or of a uniform MQ distribution",
        "fixed-key spectrum premise bar_eta_L^U <= eta_128(h,K), or the corresponding ideal-transcript event",
        "direct PXL/Macaulay cost model with the complete five components",
    ]
    structured_premises = [
        "exact-public-preflight: rank(Q_R)=K and rank(C_R,v)=M-K, giving full affine rank M",
        "fixed-key-spectrum: bar_eta_L^U <= eta_128(h,K), or the corresponding ideal-transcript event",
        "H_struct-degree: the instancewise structured system has the displayed solving degree D",
        "H_struct-Hilbert: its residual dimensions equal the displayed Hilbert-model values A and B",
        "H_struct-rank: its preprocessing and linearization Macaulay matrices have the modeled ranks",
        "H_struct-completeness: the modeled route extracts and verifies an accepted root whenever one is present",
    ]
    for shape in ELL4_SHAPES:
        v, o, ell, r = (
            int(shape["v"]),
            int(shape["o"]),
            int(shape["ell"]),
            int(shape["r"]),
        )
        K = 10 * int(shape["m1"])
        target_codimension = 0
        for omega in OMEGAS:
            supported_candidates = []
            for h in range(2, K + 1):
                probability = ell4_probability(v, o, ell, r, h, K)
                for k in range(1, h):
                    if h - k <= v:
                        supported_candidates.append(
                            (pxl_parameters(K, h, k), probability, target_codimension)
                        )
            supported_entry = best_entry(
                supported_candidates,
                family="ell4_supported",
                shape=shape,
                omega=omega,
                premises=supported_premises,
                hypothesis_metadata={
                    "H_off": {
                        "status": "unproved instancewise offset premise",
                        "lambda": int(shape["M"]) - K,
                        "rank_C_restricted_to_public_vinegar": int(shape["M"]) - K,
                        "A_span_of_kernel_is_public_vinegar": True,
                    },
                    "fixed_map_support": (
                        "homogeneous coefficient support on the selected main graph only"
                    ),
                    "full_affine_coefficient_support": False,
                },
                target_charge_basis=(
                    "stipulated 2^32 charge with consistency codimension zero under the "
                    "H_off full-rank value lambda=M-K"
                ),
            )
            supported_entry["hypothesis_metadata"] = {
                **ell4_public_premise_metadata(
                    shape,
                    K,
                    supported_entry["PXL_parameters"]["h"],
                ),
                **supported_entry["hypothesis_metadata"],
            }
            entries.append(supported_entry)

            h_struct = K - 2
            structured_probability = ell4_probability(v, o, ell, r, h_struct, K)
            structured_candidates = (
                (pxl_parameters(K, h_struct, k), structured_probability, target_codimension)
                for k in range(1, h_struct)
            )
            entries.append(
                best_entry(
                    structured_candidates,
                    family="ell4_H_struct",
                    shape=shape,
                    omega=omega,
                    premises=structured_premises,
                    hypothesis_metadata={
                        **ell4_public_premise_metadata(shape, K, h_struct),
                        "solver_premise": {
                            "name": "H_struct",
                            "status": "unproved instancewise structured-family estimate",
                            "degree_model": True,
                            "Hilbert_dimension_model": True,
                            "Macaulay_rank_model": True,
                            "solver_completeness_model": True,
                        },
                    },
                    target_charge_basis=(
                        "stipulated 2^32 sensitivity convention with no target-consistency "
                        "multiplier; this convention is not derived from H_struct"
                    ),
                )
            )
    return entries


def assert_regressions(entries: list[dict[str, Any]]) -> None:
    observed = {}
    for entry in entries:
        shape = entry["shape"]
        key = (
            entry["route_family"],
            shape["level"],
            shape["v"],
            shape["o"],
            entry["omega"],
        )
        parameters = entry["PXL_parameters"]
        observed[key] = (
            parameters["k"],
            parameters["D"],
            round(entry["work"]["reported_work_log2_binary_gates"], 6),
        )

        storage = entry["storage"]
        coefficient_slots = int(storage["exact_F19_coefficient_slots"])
        exact_bits = int(storage["exact_bits"])
        assert storage["bits_per_F19_element"] == (Q - 1).bit_length() == 5
        assert exact_bits == 5 * coefficient_slots
        assert abs(storage["log2_bits"] - log2_int(exact_bits)) < 1e-12
        assert entry["work"]["target_charge_status"].startswith(
            "stipulated unit-leading model convention"
        )

        if entry["route_family"] == "ell4_H_struct":
            premise_text = json.dumps(
                {
                    "route_premises": entry["route_premises"],
                    "hypothesis_metadata": entry["hypothesis_metadata"],
                }
            ).lower()
            assert set(entry["hypothesis_metadata"]) == {
                "exact_public_preflights",
                "fixed_key_spectrum_premise",
                "solver_premise",
            }
            assert entry["hypothesis_metadata"]["solver_premise"]["name"] == "H_struct"
            assert any(
                premise.startswith("H_struct-") for premise in entry["route_premises"]
            )
            assert "graph" not in premise_text
            assert "support" not in premise_text
            assert "injectiv" not in premise_text

        if entry["route_family"] == "ell4_supported":
            metadata = entry["hypothesis_metadata"]
            assert metadata["full_affine_coefficient_support"] is False
            assert metadata["H_off"]["lambda"] == (
                int(shape["M"]) - parameters["K"]
            )
            assert metadata["H_off"]["rank_C_restricted_to_public_vinegar"] == (
                metadata["H_off"]["lambda"]
            )
            assert metadata["H_off"]["A_span_of_kernel_is_public_vinegar"] is True

        if entry["work"]["target_consistency_codimension"] == 0:
            metadata = entry["hypothesis_metadata"]
            exact = metadata["exact_public_preflights"]
            spectrum = metadata["fixed_key_spectrum_premise"]
            assert exact["quotient_rank"]["value"] == parameters["K"]
            assert exact["affine_source_rank"]["value"] == (
                int(shape["M"]) - parameters["K"]
            )
            assert exact["full_affine_rank"]["value"] == int(shape["M"])
            assert spectrum["h"] == parameters["h"]
            assert spectrum["K"] == parameters["K"]
            assert "not an exact public preflight" in spectrum["status"]
    if observed != EXPECTED_ROUNDED:
        missing = EXPECTED_ROUNDED.keys() - observed.keys()
        extra = observed.keys() - EXPECTED_ROUNDED.keys()
        mismatches = {
            key: {"expected": EXPECTED_ROUNDED[key], "observed": observed.get(key)}
            for key in EXPECTED_ROUNDED.keys() & observed.keys()
            if EXPECTED_ROUNDED[key] != observed[key]
        }
        raise AssertionError(
            f"PXL regression mismatch; missing={missing}, extra={extra}, mismatches={mismatches}"
        )


def build_ledger() -> dict[str, Any]:
    entries = ell2_entries() + ell4_entries()
    entries.sort(
        key=lambda entry: (
            entry["route_family"],
            {"I": 1, "III": 3, "V": 5}[entry["shape"]["level"]],
            entry["shape"]["v"],
            entry["shape"]["o"],
            float(entry["omega"]),
        )
    )
    assert_regressions(entries)
    return {
        "schema": "snova-pxl-route-ledger-v2",
        "verification_metadata": {
            "generator": "pxl_final_ledger.py",
            "verification_mode": (
                "independent formula recomputation followed by byte-for-byte JSON comparison"
            ),
            "row_counts": {
                "ell2_supported": 6,
                "ell4_supported": 12,
                "ell4_H_struct": 12,
                "total": 30,
            },
            "manuscript_table_dependency": False,
            "numerical_regression_precision_decimal_places": 6,
        },
        "field": "F_19",
        "labels": {
            "exact_arithmetic": [
                "dimension arithmetic and chosen integer parameters",
                "Hilbert-series coefficients A and B once the PXL model is selected",
                "finite-field probability fractions and retry fractions",
                "integer-valued PXL components",
                "integer evaluation of the stipulated target-charge convention",
                "the stated five-bit fixed-width storage proxy",
            ],
            "estimates_or_premises": [
                "slice support, affine rank, graph injectivity, and spectrum premises",
                "degree-of-regularity and Hilbert/Macaulay transfer to the sliced systems",
                "B^omega and A^omega when omega=2.81",
                "150 binary gates per F_19 operation",
                "the stipulated 2^32 target-generation charge",
                "all total binary-gate work exponents",
                "C_aux, which is left symbolic and excluded from the reported total",
            ],
        },
        "global_formulas": {
            "D": "least D>=2 with [t^D](1-t)^(K-h+k-1)*(1+t)^K <= 1",
            "A": "sum_{d=0}^D max([t^d](1-t)^(K-h+k)*(1+t)^K,0)",
            "B": "binom(h-k+D,D)",
            "C_PXL": "B^omega + k^2*(h-k)*B^2 + k^2*A*B*binom(h+D,D) + 19^k*A^2*binom(k+D,D) + 19^k*A^omega",
            "C_target": "stipulated model convention 2^32*19^(target-consistency codimension)",
            "reported_work": "(150*C_PXL+C_target)/p_trial",
            "augmented_work": "(150*C_PXL+C_target+C_aux)/p_trial",
            "storage": "ceil(log2(19))*B^2*binom(k+D,D) = 5*B^2*binom(k+D,D) bits",
        },
        "manuscript_formula_labels": {
            "ell2_acceptance": "eq:pxl-alpha",
            "spectrum_threshold": "eq:sq-eta128",
            "ell2_trial_probability": "eq:pxl-pdir",
            "ell4_acceptance": "eq:l4-acceptance",
            "supported_slice": "eq:pxl-xz-slice",
            "offset_rank_checks": "eq:pxl-offset-checks",
            "A": "eq:pxl-A",
            "D": "eq:pxl-degree",
            "C1": "eq:pxl-C1",
            "C2": "eq:pxl-C2",
            "C3": "eq:pxl-C3",
            "Cfix": "eq:pxl-Cfix",
            "Clin": "eq:pxl-Clin",
            "target_and_total": "eq:pxl-kernel-cost",
            "reported_work": "eq:pxl-estimated-work",
            "augmented_work": "eq:pxl-augmented-work",
            "storage": "eq:pxl-storage",
        },
        "entries": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write", action="store_true", help="regenerate the JSON ledger")
    action.add_argument("--check", action="store_true", help="check the JSON ledger (default)")
    args = parser.parse_args()

    fresh = build_ledger()
    rendered = json.dumps(fresh, indent=2, sort_keys=True) + "\n"
    if args.write:
        LEDGER.write_text(rendered, encoding="utf-8")
        print(f"wrote {LEDGER.name} ({len(fresh['entries'])} rows)")
        return

    if not LEDGER.exists():
        raise SystemExit(f"missing {LEDGER}; run {Path(__file__).name} --write")
    checked_in = LEDGER.read_text(encoding="utf-8")
    if checked_in != rendered:
        raise SystemExit(
            f"{LEDGER.name} is stale; run {Path(__file__).name} --write and inspect the diff"
        )
    print(f"verified {LEDGER.name} ({len(fresh['entries'])} rows)")


if __name__ == "__main__":
    main()
