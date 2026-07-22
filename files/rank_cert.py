#!/usr/bin/env python3
"""
Exact rank of a matrix over F_p (default p=19) by fraction-free / mod-p
Gaussian elimination.  This is the tool for concern #3: certifying that the
symmetry-reduced feature map E_R-bar and the affine-constraint matrix W_R L_R
attain their generic ranks (rho = K, lambda_R = M - K) for the Level-III,
Level-V, and v2.4 parameter sets, not just Level-I.

Input: a 2-D integer array as .npy, or a whitespace/comma-separated text file.
Output: exact rank over F_p.

This computes rank EXACTLY over the prime field (no floating point), so a
reported rank is a certificate, not an estimate.
"""
import argparse, numpy as np

def rank_mod_p(A, p=19):
    A = (np.asarray(A, dtype=np.int64) % p).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c] % p != 0:
                piv = i; break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p-2, p)          # Fermat inverse
        A[r] = (A[r] * inv) % p
        for i in range(rows):
            if i != r and A[i, c] % p != 0:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows:
            break
    return r

def load(path):
    if path.endswith(".npy"):
        return np.load(path)
    return np.loadtxt(path, delimiter="," if path.endswith(".csv") else None)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix", help=".npy / .csv / whitespace-delimited integer matrix")
    ap.add_argument("-p", type=int, default=19)
    a = ap.parse_args()
    M = load(a.matrix)
    print(f"shape={M.shape}  rank_over_F{a.p} = {rank_mod_p(M, a.p)}")
