#!/usr/bin/env python3
"""Degree-one Macaulay full-column-rank check for the Level-I 50-in-102 residual.

The paper claims the degree-one Macaulay matrix has full column rank 5150
(= m*(n+1) = 50*103): the shifted generators {f_i} U {x_j f_i} are F_19-linearly
independent, i.e. NO unexpected degree-one syzygy.

A literal monomial Macaulay matrix is 187460 x 5150 (C(105,3) monomials). Instead
we compute the SAME rank memory-lightly by EVALUATION: since every generator has
total degree <= 3 < q=19, formal independence over F_19 equals functional
independence (a nonzero deg<=3 poly cannot vanish on all of F_19^n). So the rank
of the 5150 x P evaluation matrix at P > 5150 random points equals the formal
column rank. Exact over F_19 (no floating point).
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np

Q = 19
ROOT = Path(__file__).resolve().parents[1]
REPRO = ROOT/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'/'repro'
sys.path.insert(0, str(REPRO)); sys.path.insert(0, str(ROOT/'experiments'))
import symmetry_attack_validation as V   # noqa: E402
import emit_level1_core as EM            # noqa: E402
from task3_lowdegree import rank_mod     # noqa: E402  (vectorized mod-p RREF rank)


def eval_forms(polys, X):
    """Evaluate each quadratic f=(A,lin,const) at every row of X (P x n). -> (m x P)."""
    P = X.shape[0]; m = len(polys)
    F = np.zeros((m, P), dtype=np.int64)
    for i, (A, lin, const) in enumerate(polys):
        XA = (X @ A) % Q                      # P x n
        quad = np.einsum('pj,pj->p', XA, X) % Q
        F[i] = (quad + (X @ lin) + int(const)) % Q
    return F % Q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kat', default=str(ROOT/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'
                                         /'source_snapshots'/'PQCsignKAT_SNOVA_28_5_19_4.txt'))
    ap.add_argument('--seed', type=int, default=20260721)
    ap.add_argument('--points', type=int, default=6000)
    args = ap.parse_args()

    print("[*] building official 50-in-102 residual from KAT key ...", flush=True)
    data = EM.build_official_residual(V, Path(args.kat), args.seed)
    assert EM.verify_polys(data), "official residual FAILED reference verifier"
    polys = data['polys']; n = data['residual_vars']; m = len(polys)
    ncols = m*(n+1)
    print(f"    verifier net PASSED; residual {m}/{n}; expected full column rank = m*(n+1) = {ncols}", flush=True)

    rng = np.random.default_rng(args.seed)
    P = max(args.points, ncols + 200)
    X = rng.integers(0, Q, size=(P, n), dtype=np.int64)
    F = eval_forms(polys, X)                       # m x P  (values of f_i)

    # shifted generators: rows {f_i} then {x_j * f_i}
    G = np.empty((ncols, P), dtype=np.int64)
    G[:m] = F
    r = m
    for j in range(n):
        G[r:r+m] = (F * X[:, j][None, :]) % Q
        r += m
    assert r == ncols

    print(f"[*] evaluation matrix {G.shape} over F_{Q}; computing exact rank ...", flush=True)
    rk = rank_mod(G % Q)
    print(f"[*] degree-one Macaulay column rank = {rk}  (expected {ncols}; "
          f"{'FULL RANK -> no unexpected degree-one syzygy' if rk == ncols else 'DEFICIENT -> syzygy present'})",
          flush=True)


if __name__ == '__main__':
    main()
