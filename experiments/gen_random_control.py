#!/usr/bin/env python3
"""Random control systems for the solving-degree experiment.

Generates N random systems of `m` quadratics in `n` variables over F_19 with:
  * matched quadratic support density (each z_i z_j / z_i^2 term present with prob
    = target density, matching the official core), and
  * a planted solution (constant term set so a random point is a root -> same
    consistency status as the official core, which has a known root).

Written in the same msolve (.ms) and Magma (.mag) format as the official core so
the identical solver command applies.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

Q = 19

def gen_one(m, n, density, rng):
    root = rng.integers(0, Q, size=n, dtype=np.int64)
    polys = []
    for _ in range(m):
        A = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            for j in range(i, n):
                if rng.random() < density:
                    A[i, j] = rng.integers(1, Q)   # nonzero coeff on present monomial
        lin = np.array([rng.integers(0, Q) if rng.random() < density else 0
                        for _ in range(n)], dtype=np.int64)
        # constant so that planted root is a solution
        val = (int(root @ np.triu(A) @ root) + int(lin @ root)) % Q
        const = (-val) % Q
        polys.append((A % Q, lin % Q, const))
    return polys, root

def poly_terms(A, lin, const, n):
    terms = []
    for i in range(n):
        for j in range(i, n):
            c = int(A[i, j]) % Q
            if c: terms.append(f"{c}*x{i}*x{j}")
    for i in range(n):
        c = int(lin[i]) % Q
        if c: terms.append(f"{c}*x{i}")
    if const % Q: terms.append(f"{int(const)%Q}")
    return "+".join(terms) if terms else "0"

def write_ms(polys, n, path):
    vs = ",".join(f"x{i}" for i in range(n))
    body = ",\n".join(poly_terms(A, lin, c, n) for (A, lin, c) in polys)
    path.write_text(f"{vs}\n{Q}\n{body}\n")

def write_mag(polys, n, path):
    head = [f"q := {Q};", f"n := {n};",
            f"R<{','.join('x'+str(i) for i in range(n))}> := PolynomialRing(GF(q), n, \"grevlex\");",
            "F := ["]
    fs = ",\n".join("  " + poly_terms(A, lin, c, n) for (A, lin, c) in polys)
    path.write_text("\n".join(head) + "\n" + fs + "\n];\n")

def support_density(polys, n):
    maxq = n*(n+1)//2
    return float(np.mean([sum(1 for i in range(n) for j in range(i,n) if int(A[i,j])%Q)/maxq
                          for (A,_,_) in polys]))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--m', type=int, required=True)
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--density', type=float, required=True)
    ap.add_argument('--count', type=int, default=5)
    ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--tag', default='control')
    args = ap.parse_args()
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    meta = []
    for k in range(args.count):
        polys, root = gen_one(args.m, args.n, args.density, rng)
        base = outdir / f"{args.tag}_{args.m}in{args.n}_{k:02d}"
        write_ms(polys, args.n, base.with_suffix('.ms'))
        write_mag(polys, args.n, base.with_suffix('.mag'))
        meta.append(dict(index=k, m=args.m, n=args.n,
                         density=support_density(polys, args.n),
                         planted_root=[int(v) for v in root]))
    (outdir / f"{args.tag}_{args.m}in{args.n}_meta.json").write_text(json.dumps(meta, indent=2)+"\n")
    print(f"wrote {args.count} random {args.m}/{args.n} controls, mean density "
          f"{np.mean([mm['density'] for mm in meta]):.4f}")

if __name__ == '__main__':
    main()
