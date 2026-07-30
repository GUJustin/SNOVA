#!/usr/bin/env python3
"""Emit the Level-I l=4 (28,5,19,4) determined core (50 quadratics / 41 vars)
FROM THE OFFICIAL KAT KEY, as an msolve system with a planted solution.

Pipeline (reuses repro/symmetry_attack_validation.py as a library):
  1. Reconstruct the official public key from the KAT seed (byte-for-byte checked).
  2. Build the symmetry-reduced affine system F(x)=f0+Lmat x+E c(x), pick a random
     target from a planted x*, eliminate the 30 affine constraints -> x=x0+T z with
     102 residual free vars z; planted residual z* = x*[free].
  3. Express the 50 residual quadratics g_r(z) explicitly (each feature is a quadratic
     form x^T M_j x; substitute x=x0+T z).  VERIFY the emitted polynomials against the
     reference direct_output on random z (exact over F_19) -- this is the correctness net.
  4. Specialize (102-CORE_N) residual coords to their planted values, leaving 50
     quadratics in CORE_N=41 vars with the same planted root (consistency preserved).
  5. Write msolve (.ms) and Magma (.mag) inputs.

No number is trusted unless step 3's verification passes.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np

Q = 19

def load_lib(repro_dir: Path):
    sys.path.insert(0, str(repro_dir))
    import symmetry_attack_validation as V
    return V

def build_official_residual(V, kat_path: Path, seed: int):
    """Reproduce validate_level1's construction and return the residual-system data."""
    p = V.PARAMS[0]  # I-square-l4 (28,5,4,4)
    kat = V.parse_kat(kat_path)
    seed48 = bytes.fromhex(kat['sk']); expected = bytes.fromhex(kat['pk'])
    pk, P = V.kat_public_key(seed48, p)
    assert pk == expected, "KAT public key mismatch -- not the official key"
    abq = V.reconstruct_abq(p)
    rng = np.random.default_rng(seed)
    rho = np.array([1, 8, 9, 14], dtype=np.int64)
    E, labels = V.build_E(p, abq, rho); er = V.rank_mod(E)
    H = V.nullspace_mod(E.T)
    Voff, fmt = V.safe_offsets(p, rho, rng)
    # affine decomposition F(x) = f0 + Lmat x + E c(x)
    f0 = V.direct_output(p, abq, P, Voff)
    Lmat = np.zeros((p.outputs, p.variables), dtype=np.int64)
    for j in range(p.variables):
        e = np.zeros(p.variables, dtype=np.int64); e[j] = 1
        Lmat[:, j] = (V.direct_output(p, abq, P, V.relation_U(p, e, rho, Voff))
                      - f0 - E @ V.coords(p, P, e, labels)) % Q
    # planted solution and target
    xstar = rng.integers(0, Q, size=p.variables, dtype=np.int64)
    target = V.direct_output(p, abq, P, V.relation_U(p, xstar, rho, Voff))
    # eliminate affine constraints
    C = H @ Lmat % Q; b = H @ (target - f0) % Q
    x0, T, free, piv = V.affine_solve(C, b)
    zstar = xstar[free]
    assert np.array_equal((x0 + T @ zstar) % Q, xstar), "planted residual mismatch"
    rows = V.rref_mod(E.T)[1][:er]  # 50 rows making E square-invertible
    # feature quadratic-form matrices M_j : c_j(x) = x^T M_j x
    Sps = V.powers(V.official_S(p.l), p.l)
    KS = [V.kronS(A, p.n) for A in Sps]           # kron(I_n, S^a)
    Mlist = []
    for (pp, a, bl) in labels:
        Mlist.append((KS[a].T @ P[pp] @ KS[bl]) % Q)   # x^T M x = W[a]^T P W[b]
    # residual quadratics g_r(z), r in rows: substitute x = x0 + T z
    # g_r(x) = (f0[r]-target[r]) + Lmat[r].x + sum_j E[r,j] x^T M_j x
    polys = []  # each: (Qmat 102x102, lin 102, const) with g = z^T Qmat z + lin.z + const
    for r in rows:
        Gr = np.zeros((p.variables, p.variables), dtype=np.int64)
        Erow = E[r]
        for j in range(len(Mlist)):
            c = int(Erow[j]) % Q
            if c: Gr = (Gr + c * Mlist[j]) % Q
        # g_r(x) = x^T Gr x + Lmat[r].x + (f0[r]-target[r])
        A = (T.T @ Gr @ T) % Q                                  # z^T A z
        lin_x = Lmat[r].astype(np.int64)
        # linear in z: from 2*x0^T Gr T (quadratic cross) + lin_x @ T
        lin_z = ((x0 @ (Gr + Gr.T) % Q) @ T + lin_x @ T) % Q
        const = (int(x0 @ Gr @ x0 % Q) + int(lin_x @ x0 % Q)
                 + int(f0[r]) - int(target[r])) % Q
        polys.append((A % Q, lin_z % Q, const % Q))
    return dict(p=p, free=free, x0=x0, T=T, zstar=zstar, rows=rows,
                polys=polys, E=E, Lmat=Lmat, f0=f0, target=target,
                abq=abq, P=P, rho=rho, Voff=Voff, labels=labels, V=V,
                kat_match=True, er=er, left_kernel=int(H.shape[0]),
                lin_rank=len(piv), residual_vars=T.shape[1])

def eval_poly(A, lin, const, z):
    return (int(z @ A @ z) + int(lin @ z) + int(const)) % Q

def verify_polys(data, trials=200):
    """Check emitted polynomials == reference verifier residual on random z, over F_19."""
    V = data['V']; p = data['p']; rng = np.random.default_rng(12345)
    x0, T, rows = data['x0'], data['T'], data['rows']
    E, f0, target = data['E'], data['f0'], data['target']
    ok = True
    for _ in range(trials):
        z = rng.integers(0, Q, size=T.shape[1], dtype=np.int64)
        x = (x0 + T @ z) % Q
        # reference residual g_r via direct verifier: F(x)-target at rows
        g = (V.direct_output(p, data['abq'], data['P'],
                             V.relation_U(p, x, data['rho'], data['Voff'])) - target) % Q
        for idx, r in enumerate(rows):
            A, lin, const = data['polys'][idx]
            if eval_poly(A, lin, const, z) != int(g[r]):
                ok = False; break
        if not ok: break
    return ok

def planted_check(data):
    z = data['zstar']
    return all(eval_poly(A, lin, const, z) == 0 for (A, lin, const) in data['polys'])

def specialize(data, core_n, seed=777):
    """Fix (residual_vars - core_n) coords to planted values -> core_n free vars."""
    rng = np.random.default_rng(seed)
    nvar = data['residual_vars']; zstar = data['zstar']
    keep = sorted(rng.choice(nvar, size=core_n, replace=False).tolist())
    keepset = set(keep)
    fixed = {i: int(zstar[i]) for i in range(nvar) if i not in keepset}
    newpolys = []
    for (A, lin, const) in data['polys']:
        A = A % Q; lin = lin.copy() % Q; const = int(const)
        # substitute fixed z_i = val
        c2 = const; lin2 = np.zeros(core_n, dtype=np.int64)
        Anew = np.zeros((core_n, core_n), dtype=np.int64)
        idx = {v: k for k, v in enumerate(keep)}
        for i in range(nvar):
            for j in range(nvar):
                a = int(A[i, j]) % Q
                if not a: continue
                fi = i in keepset; fj = j in keepset
                if fi and fj:
                    Anew[idx[i], idx[j]] = (Anew[idx[i], idx[j]] + a) % Q
                elif fi and not fj:
                    lin2[idx[i]] = (lin2[idx[i]] + a * fixed[j]) % Q
                elif (not fi) and fj:
                    lin2[idx[j]] = (lin2[idx[j]] + a * fixed[i]) % Q
                else:
                    c2 = (c2 + a * fixed[i] * fixed[j]) % Q
        for i in range(nvar):
            li = int(lin[i]) % Q
            if not li: continue
            if i in keepset: lin2[idx[i]] = (lin2[idx[i]] + li) % Q
            else: c2 = (c2 + li * fixed[i]) % Q
        newpolys.append((Anew % Q, lin2 % Q, c2 % Q))
    zc = np.array([int(zstar[i]) for i in keep], dtype=np.int64)
    # planted root must still satisfy
    assert all(eval_poly(A, lin, const, zc) == 0 for (A, lin, const) in newpolys), \
        "planted root lost under specialization"
    return newpolys, keep, zc

def poly_terms(A, lin, const, nvar):
    """Return list of (coeff, monomial-string) for z^T A z + lin.z + const over F_19."""
    terms = []
    for i in range(nvar):
        for j in range(i, nvar):
            if i == j:
                c = int(A[i, i]) % Q
            else:
                c = (int(A[i, j]) + int(A[j, i])) % Q
            if c:
                terms.append((c, f"x{i}*x{j}"))
    for i in range(nvar):
        c = int(lin[i]) % Q
        if c: terms.append((c, f"x{i}"))
    if const % Q: terms.append((int(const) % Q, None))
    return terms

def to_msolve(polys, nvar, path: Path):
    varnames = ",".join(f"x{i}" for i in range(nvar))
    lines = [varnames, str(Q)]
    body = []
    for (A, lin, const) in polys:
        terms = poly_terms(A, lin, const, nvar)
        if not terms: terms = [(0, None)]
        s = "+".join((f"{c}*{m}" if m else f"{c}") for (c, m) in terms)
        body.append(s)
    path.write_text(",\n".join(body) + "\n" if False else
                    lines[0] + "\n" + lines[1] + "\n" + ",\n".join(body) + "\n")

def to_magma(polys, nvar, path: Path):
    out = [f"q := {Q};", f"n := {nvar};",
           f"R<{','.join('x'+str(i) for i in range(nvar))}> := PolynomialRing(GF(q), n, \"grevlex\");",
           "F := ["]
    fs = []
    for (A, lin, const) in polys:
        terms = poly_terms(A, lin, const, nvar)
        if not terms: terms = [(0, None)]
        s = "+".join((f"{c}*{m}" if m else f"{c}") for (c, m) in terms)
        fs.append("  " + s)
    out.append(",\n".join(fs))
    out.append("];")
    path.write_text("\n".join(out) + "\n")

def main():
    ap = argparse.ArgumentParser()
    root = Path(__file__).resolve().parents[1]
    ap.add_argument('--repro', default=str(root/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'/'repro'))
    ap.add_argument('--kat', default=str(root/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'/'source_snapshots'/'PQCsignKAT_SNOVA_28_5_19_4.txt'))
    ap.add_argument('--core-n', type=int, default=41)
    ap.add_argument('--seed', type=int, default=20260721)
    ap.add_argument('--outdir', default=str(root/'experiments'/'systems'))
    ap.add_argument('--tag', default='levelI_l4')
    args = ap.parse_args()
    V = load_lib(Path(args.repro))
    data = build_official_residual(V, Path(args.kat), args.seed)
    print(f"[*] official key match: {data['kat_match']}  residual = {data['er']} quadratics / {data['residual_vars']} vars")
    print(f"[*] verifying emitted polynomials against reference verifier (200 random points, exact F_19)...")
    assert verify_polys(data), "EMITTED POLYNOMIALS DO NOT MATCH REFERENCE -- abort"
    print("    verification PASSED")
    assert planted_check(data), "planted root does not satisfy residual"
    print(f"    planted root satisfies all {data['er']} residual quadratics")
    # full 50/102 residual
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    to_msolve(data['polys'], data['residual_vars'], outdir/f"{args.tag}_residual_50in102.ms")
    # core 50/CORE_N
    newpolys, keep, zc = specialize(data, args.core_n, seed=args.seed)
    to_msolve(newpolys, args.core_n, outdir/f"{args.tag}_core_50in{args.core_n}.ms")
    to_magma(newpolys, args.core_n, outdir/f"{args.tag}_core_50in{args.core_n}.mag")
    # record planted solution + support density
    def density(polys, nvar):
        maxq = nvar*(nvar+1)//2
        dens = [sum(1 for i in range(nvar) for j in range(i,nvar)
                    if ((int(A[i,i])%Q) if i==j else ((int(A[i,j])+int(A[j,i]))%Q)))/maxq
                for (A,_,_) in polys]
        return float(np.mean(dens))
    meta = dict(tag=args.tag, official_key=True, residual="50/102",
                core=f"50/{args.core_n}", planted_core_solution=[int(v) for v in zc],
                kept_coords=keep, quad_support_density_core=density(newpolys, args.core_n),
                quad_support_density_residual=density(data['polys'], data['residual_vars']))
    (outdir/f"{args.tag}_core_meta.json").write_text(json.dumps(meta, indent=2)+"\n")
    print(f"[*] wrote residual + core systems and meta to {outdir}")
    print(f"[*] core quadratic support density = {meta['quad_support_density_core']:.4f}")

if __name__ == '__main__':
    main()
