#!/usr/bin/env python3
"""Exact conditional-atom ledger for the official ell=2 Frobenius channels.

Raw F_19 symbols are independent with integer weights 14 on residues 0..8
and 13 on residues 9..18, corresponding to one uniform byte reduced modulo
19.  The script diagonalizes the official S matrix in F_19[u]/(u^2+1), uses
the resulting explicit block transforms, and exhaustively computes:

* max_h Pr[H=h | D=d] for every diagonal native 2x2 symmetric block;
* max_z Pr[l(H)=z | D=d] for every off-diagonal native 2x2 block and every
  nonzero projective F_19-linear functional l on H in F_19^2.
"""
from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

Q = 19
WEIGHT = tuple(14 if x <= 8 else 13 for x in range(Q))
HERE = Path(__file__).resolve().parent
ART = HERE.parent


def diagonal_ledger():
    by_d: dict[tuple[int, int], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for a in range(Q):
        for b in range(Q):
            for c in range(Q):
                d_pair = ((a + 7*b + 16*c) % Q, (2*b + 7*c) % Q)
                h = (a + 7*b + 18*c) % Q
                by_d[d_pair][h] += WEIGHT[a] * WEIGHT[b] * WEIGHT[c]
    maximum = Fraction(0)
    witness = None
    for d_pair, dist in by_d.items():
        denominator = sum(dist.values())
        for h, numerator in dist.items():
            value = Fraction(numerator, denominator)
            if value > maximum:
                maximum = value
                witness = {
                    "D": list(d_pair), "H": h,
                    "numerator_weight": numerator,
                    "conditional_denominator_weight": denominator,
                }
    return maximum, witness


def off_diagonal_ledger():
    by_d: dict[tuple[int, int], dict[tuple[int, int], int]] = defaultdict(lambda: defaultdict(int))
    for a in range(Q):
        for b in range(Q):
            for c in range(Q):
                for d in range(Q):
                    d_pair = ((a + 13*b + 13*c + 16*d) % Q,
                              (b + c + 7*d) % Q)
                    h_pair = ((a + 13*b + 13*c + 18*d) % Q,
                              (18*b + c) % Q)
                    by_d[d_pair][h_pair] += WEIGHT[a] * WEIGHT[b] * WEIGHT[c] * WEIGHT[d]

    # Projective representatives for nonzero dual vectors in F_19^2.
    functionals = [(1, beta) for beta in range(Q)] + [(0, 1)]
    maximum = Fraction(0)
    witness = None
    for d_pair, dist in by_d.items():
        denominator = sum(dist.values())
        for alpha, beta in functionals:
            values: dict[int, int] = defaultdict(int)
            for (h0, h1), weight in dist.items():
                values[(alpha*h0 + beta*h1) % Q] += weight
            for z, numerator in values.items():
                value = Fraction(numerator, denominator)
                if value > maximum:
                    maximum = value
                    witness = {
                        "D": list(d_pair),
                        "functional": [alpha, beta],
                        "value": z,
                        "numerator_weight": numerator,
                        "conditional_denominator_weight": denominator,
                    }
    return maximum, witness


def main():
    diag, diag_witness = diagonal_ledger()
    off, off_witness = off_diagonal_ledger()
    assert diag == Fraction(49, 829)
    assert off == Fraction(169246, 2971565)
    assert off < diag
    output = {
        "field": "F_19[u]/(u^2+1)",
        "official_S": [[1, 2], [2, 15]],
        "eigenvalues": ["8+2u", "8-2u"],
        "eigenvectors": ["(1,13+u)", "(1,13-u)"],
        "native_symbol_integer_weights": {"0..8": 14, "9..18": 13, "normalizer": 256},
        "diagonal_block_transform": {
            "D_real": "a+7b+16c", "D_u": "2b+7c", "H": "a+7b+18c"
        },
        "off_diagonal_block_transform": {
            "D_real": "a+13b+13c+16d", "D_u": "b+c+7d",
            "H_real": "a+13b+13c+18d", "H_u": "18b+c"
        },
        "diagonal_max_conditional_atom": {
            "fraction": [diag.numerator, diag.denominator], "float": float(diag),
            "witness": diag_witness,
        },
        "off_diagonal_max_conditional_functional_atom": {
            "fraction": [off.numerator, off.denominator], "float": float(off),
            "witness": off_witness,
        },
        "common_upper_bound": {
            "fraction": [diag.numerator, diag.denominator], "float": float(diag)
        },
    }
    path = ART / "l2_channel_conditional_atom_ledger.json"
    path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
