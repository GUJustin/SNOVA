#!/usr/bin/env python3
"""Task 3: reproduce the paper's "Low-degree sanity check" statistics for the
official Level-I 50-in-102 residual, from the OFFICIAL key (no shipped script
exists for these numbers, so they are otherwise unverified).

Reuses emit_level1_core.build_official_residual, whose emitted quadratics are
verified against the reference verifier (exact over F_19) before any statistic
is computed. The 20-key arm builds residuals from 20 independent 48-byte seeds
through the SAME code path and re-runs the identical verifier net per key.

Statistics reproduced (paper lines 766-784):
  * planted-root Jacobian rank (claim: 50)
  * 10,000 random coordinate 50-strips invertible (claim: 9,466; rmt 0.9445987429)
  * 20-key arm: every planted root rank 50; strips full rank per key range 9,412-9,491 (mean 9,447)
  * polar matrices full output rank 50 for all 102 coord dirs + root + 10,000 random
  * Hessian rank over all Hamming-weight<=2 coeff vectors (22,100) + 20,000 random
    (claim: minimum observed rank 100)
Every reported number comes from this run; none is copied from the paper.
"""
from __future__ import annotations
import sys, json, argparse
from pathlib import Path
import numpy as np

Q = 19
ROOT = Path(__file__).resolve().parents[1]
REPRO = ROOT/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'/'repro'
sys.path.insert(0, str(REPRO))
sys.path.insert(0, str(ROOT/'experiments'))
import symmetry_attack_validation as V   # noqa: E402
import emit_level1_core as EM            # noqa: E402


def rank_mod(M, p=Q):
    """Exact rank over F_p by vectorized Gaussian elimination (numpy int64)."""
    A = (np.asarray(M, dtype=np.int64) % p).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        nz = np.nonzero(A[r:, c] % p)[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p-2, p)
        A[r] = (A[r]*inv) % p
        col = A[:, c].copy(); col[r] = 0
        A = (A - np.outer(col, A[r])) % p
        r += 1
        if r == rows:
            break
    return r


def build_from_seed(seed48: bytes, rng_seed: int):
    """Mirror emit_level1_core.build_official_residual but from an arbitrary key
    seed (no KAT assert). Returns the same dict shape so EM.verify_polys applies."""
    p = V.PARAMS[0]
    pk, P = V.kat_public_key(seed48, p)
    abq = V.reconstruct_abq(p)
    rng = np.random.default_rng(rng_seed)
    rho = np.array([1, 8, 9, 14], dtype=np.int64)
    E, labels = V.build_E(p, abq, rho); er = V.rank_mod(E)
    H = V.nullspace_mod(E.T)
    Voff, fmt = V.safe_offsets(p, rho, rng)
    f0 = V.direct_output(p, abq, P, Voff)
    Lmat = np.zeros((p.outputs, p.variables), dtype=np.int64)
    for j in range(p.variables):
        e = np.zeros(p.variables, dtype=np.int64); e[j] = 1
        Lmat[:, j] = (V.direct_output(p, abq, P, V.relation_U(p, e, rho, Voff))
                      - f0 - E @ V.coords(p, P, e, labels)) % Q
    xstar = rng.integers(0, Q, size=p.variables, dtype=np.int64)
    target = V.direct_output(p, abq, P, V.relation_U(p, xstar, rho, Voff))
    C = H @ Lmat % Q; b = H @ (target - f0) % Q
    x0, T, free, piv = V.affine_solve(C, b)
    zstar = xstar[free]
    assert np.array_equal((x0 + T @ zstar) % Q, xstar)
    rows = V.rref_mod(E.T)[1][:er]
    Sps = V.powers(V.official_S(p.l), p.l)
    KS = [V.kronS(A, p.n) for A in Sps]
    Mlist = [(KS[a].T @ P[pp] @ KS[bl]) % Q for (pp, a, bl) in labels]
    polys = []
    for r in rows:
        Gr = np.zeros((p.variables, p.variables), dtype=np.int64)
        Erow = E[r]
        for j in range(len(Mlist)):
            c = int(Erow[j]) % Q
            if c: Gr = (Gr + c*Mlist[j]) % Q
        A = (T.T @ Gr @ T) % Q
        lin_x = Lmat[r].astype(np.int64)
        lin_z = ((x0 @ (Gr + Gr.T) % Q) @ T + lin_x @ T) % Q
        const = (int(x0 @ Gr @ x0 % Q) + int(lin_x @ x0 % Q)
                 + int(f0[r]) - int(target[r])) % Q
        polys.append((A % Q, lin_z % Q, const % Q))
    return dict(p=p, free=free, x0=x0, T=T, zstar=zstar, rows=rows, polys=polys,
                E=E, Lmat=Lmat, f0=f0, target=target, abq=abq, P=P, rho=rho,
                Voff=Voff, labels=labels, V=V, pk=pk,
                residual_vars=T.shape[1], er=er)


def symmetric_forms(data):
    """S_r = A_r + A_r^T (the Hessian of z^T A_r z), lin_r, over F_19."""
    S = [ (A + A.T) % Q for (A, lin, c) in data['polys'] ]
    L = [ lin % Q for (A, lin, c) in data['polys'] ]
    return S, L


def jacobian_at(data, S, L, z):
    """Jacobian rows grad g_r(z) = S_r z + lin_r -> (50 x nvar)."""
    return np.array([ (S[r] @ z + L[r]) % Q for r in range(len(S)) ], dtype=np.int64) % Q


def strip_test(J, count, rng):
    """count random coordinate n_eq-subsets of columns; fraction giving invertible
    n_eq x n_eq submatrix. Returns (num_invertible, count)."""
    neq, nvar = J.shape
    good = 0
    for _ in range(count):
        cols = rng.choice(nvar, size=neq, replace=False)
        if rank_mod(J[:, cols]) == neq:
            good += 1
    return good, count


def polar_test(S, dirs):
    """For each direction d, polar matrix rows S_r d (50 x nvar); return list of ranks."""
    ranks = []
    for d in dirs:
        M = np.array([ (S[r] @ d) % Q for r in range(len(S)) ], dtype=np.int64) % Q
        ranks.append(rank_mod(M))
    return ranks


def hessian_min_rank(S, extra_random, rng):
    """Min rank of sum_r c_r S_r over all Hamming-weight<=2 coeff vectors plus
    `extra_random` uniform nonzero vectors. Returns (min_rank, n_tested)."""
    neq = len(S); nvar = S[0].shape[0]
    minr = nvar + 1; n = 0
    # weight 1
    for r in range(neq):
        minr = min(minr, rank_mod(S[r])); n += 1
    # weight 2: positions i<j, coeff (1, t) for t in 1..Q-1  (projective normal form)
    for i in range(neq):
        for j in range(i+1, neq):
            SiSj = S[i]
            for t in range(1, Q):
                H = (SiSj + t*S[j]) % Q
                minr = min(minr, rank_mod(H)); n += 1
    # random
    for _ in range(extra_random):
        c = rng.integers(0, Q, size=neq, dtype=np.int64)
        if not c.any():
            c[0] = 1
        H = np.zeros_like(S[0])
        for r in range(neq):
            if c[r]:
                H = (H + int(c[r])*S[r]) % Q
        minr = min(minr, rank_mod(H)); n += 1
    return minr, n


def rmt_invertible_prob(n, q=Q):
    p = 1.0
    for k in range(1, n+1):
        p *= (1 - q**(-k))
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--kat', default=str(ROOT/'files'/'SNOVA_Symmetry_Quotient_Reproducibility_2026-07-19'/'source_snapshots'/'PQCsignKAT_SNOVA_28_5_19_4.txt'))
    ap.add_argument('--seed', type=int, default=20260721)
    ap.add_argument('--strips', type=int, default=10000)
    ap.add_argument('--rand-dirs', type=int, default=10000)
    ap.add_argument('--rand-hess', type=int, default=20000)
    ap.add_argument('--keys', type=int, default=20)
    ap.add_argument('--out', default=str(ROOT/'results'/'task3_lowdegree.json'))
    ap.add_argument('--skip-hessian', action='store_true')
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    results = {}

    # ---- official key ----
    print("[*] building official 50-in-102 residual from KAT key ...", flush=True)
    off = EM.build_official_residual(V, Path(args.kat), args.seed)
    assert EM.verify_polys(off), "official residual FAILED reference verifier"
    assert EM.planted_check(off), "official planted root not a root"
    print(f"    verifier net PASSED; residual {off['er']}/{off['residual_vars']}", flush=True)
    S, L = symmetric_forms(off)
    nvar = off['residual_vars']; neq = len(S)

    J = jacobian_at(off, S, L, off['zstar'])
    jrank = rank_mod(J)
    print(f"[*] planted-root Jacobian rank = {jrank} (expect {neq})", flush=True)

    good, tot = strip_test(J, args.strips, rng)
    rmt = rmt_invertible_prob(neq)
    print(f"[*] strips: {good}/{tot} invertible; empirical {good/tot:.6f} vs rmt {rmt:.10f}", flush=True)

    # polar: 102 coord dirs + planted-root dir + random dirs
    dirs = [np.eye(nvar, dtype=np.int64)[i] for i in range(nvar)]
    dirs.append(off['zstar'] % Q)
    for _ in range(args.rand_dirs):
        d = rng.integers(0, Q, size=nvar, dtype=np.int64)
        if not d.any(): d[0] = 1
        dirs.append(d)
    pranks = polar_test(S, dirs)
    print(f"[*] polar matrices: {len(pranks)} dirs, min rank {min(pranks)}, max {max(pranks)} (expect all {neq})", flush=True)

    results['official'] = dict(
        key='KAT', residual=f"{neq}/{nvar}", jacobian_rank=int(jrank),
        strips_invertible=int(good), strips_total=int(tot),
        strips_fraction=good/tot, rmt_prob=rmt,
        polar_dirs=len(pranks), polar_min_rank=int(min(pranks)), polar_max_rank=int(max(pranks)),
    )

    if not args.skip_hessian:
        print("[*] Hessian scan (weight<=2 + random) ... this is the long part", flush=True)
        hmin, hn = hessian_min_rank(S, args.rand_hess, rng)
        print(f"[*] Hessian: {hn} combinations, minimum rank = {hmin} in {nvar} vars (expect 100)", flush=True)
        results['official']['hessian_tested'] = int(hn)
        results['official']['hessian_min_rank'] = int(hmin)

    # ---- 20 independently generated official-shaped keys ----
    print(f"[*] {args.keys}-key arm: Jacobian rank + strips per key ...", flush=True)
    keyrows = []
    strips_full_total = 0; strips_total = 0
    for k in range(args.keys):
        seed48 = rng.integers(0, 256, size=48, dtype=np.uint8).tobytes()
        d = build_from_seed(seed48, rng_seed=args.seed + 1000 + k)
        assert EM.verify_polys(d), f"key {k} residual FAILED reference verifier"
        Sk, Lk = symmetric_forms(d)
        Jk = jacobian_at(d, Sk, Lk, d['zstar'])
        jr = rank_mod(Jk)
        g, t = strip_test(Jk, args.strips, rng)
        strips_full_total += g; strips_total += t
        keyrows.append(dict(key=k, jacobian_rank=int(jr), strips_invertible=int(g), strips_total=int(t)))
        print(f"    key {k:2d}: Jrank={jr}  strips {g}/{t}", flush=True)
    results['multikey'] = dict(
        n_keys=args.keys,
        all_jacobian_full=all(r['jacobian_rank'] == neq for r in keyrows),
        min_jacobian_rank=min(r['jacobian_rank'] for r in keyrows),
        strips_full_total=int(strips_full_total), strips_total=int(strips_total),
        strips_fraction=strips_full_total/strips_total if strips_total else None,
        per_key=keyrows,
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2)+"\n")
    print(f"[*] wrote {args.out}", flush=True)


if __name__ == '__main__':
    main()
