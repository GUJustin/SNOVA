#!/usr/bin/env python3
"""
Semi-regular degree of regularity for a DETERMINED/OVERDETERMINED system of
quadratic equations over a field of characteristic != 2.

For m quadratics in n variables (m >= n), the semi-regular Hilbert series is
    H(t) = (1 - t^2)^m / (1 - t)^n
and the degree of regularity D_reg is the index of the first coefficient that
is <= 0 (the "first fall").  The solving degree of a Groebner/XL computation on
a genuinely semi-regular system is D_reg.

SCOPE / HONESTY NOTE:
  * This is the standard semi-regular predictor.  It applies to determined or
    overdetermined systems (e.g. the 50-in-52 cross-column normal form, or the
    determined subsystems Hashimoto's reduction produces AFTER specialization).
  * It does NOT run Hashimoto's (a,k) optimization and therefore does NOT
    reproduce the paper's end-to-end gate costs.  It reports the *predicted*
    operating degree, which is the quantity a real solving-degree experiment
    (msolve/Magma) should be compared against.
  * Field equations x^q - x are ignored (irrelevant at these low degrees for q=19).
"""
import argparse
from fractions import Fraction

def hilbert_coeffs(m, n, max_deg):
    # numerator (1 - t^2)^m  -> coefficients
    num = [0]*(2*m+1)
    from math import comb
    for j in range(m+1):
        num[2*j] += (-1)**j * comb(m, j)
    # divide by (1-t)^n : multiply by (1-t)^{-n} = sum comb(n-1+d, d) t^d
    series = [Fraction(0)]*(max_deg+1)
    for d in range(max_deg+1):
        s = Fraction(0)
        for k in range(min(d, 2*m)+1):
            if num[k] == 0:
                continue
            s += num[k]*comb(n-1+(d-k), d-k)
        series[d] = s
    return series

def d_reg(m, n, cap=200):
    if m < n:
        return None  # underdetermined: no first-fall; needs Hashimoto reduction first
    c = hilbert_coeffs(m, n, cap)
    for d in range(1, cap+1):
        if c[d] <= 0:
            return d
    return None

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("m", type=int, help="number of quadratic equations")
    ap.add_argument("n", type=int, help="number of variables")
    a = ap.parse_args()
    dr = d_reg(a.m, a.n)
    if dr is None and a.m < a.n:
        print(f"m={a.m} n={a.n}: UNDERDETERMINED (m<n); apply Hashimoto reduction to a determined subsystem first.")
    else:
        print(f"m={a.m} n={a.n}: predicted semi-regular D_reg = {dr}")
