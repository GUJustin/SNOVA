#!/usr/bin/env python3
"""Reproduce the SNOVA symmetry-quotient cost tables.

The script implements the finite-field semi-regular series and Hashimoto
underdetermined-MQ reduction used by the public SNOVA analysis.  It emits both
the p1=0 convention matching Version 2.3 Table 11 and the stricter p1=1
convention used for the paper's headline values.
"""
from __future__ import annotations
import argparse, functools, json, math
from pathlib import Path

PRECISION = 300


def reg_coeff(d: int, n: int, m: int, q: int) -> int:
    total = 0
    for c in range(d // (2*q) + 1):
        cc = math.comb(m+c-1, c) if c else 1
        rem1 = d - 2*q*c
        for a in range(min(m, rem1//2) + 1):
            ca = (-1)**a * math.comb(m, a)
            e = rem1 - 2*a
            inner = 0
            for b in range(min(n, e//q) + 1):
                ee = e - b*q
                inner += (-1)**b * math.comb(n, b) * math.comb(n+ee-1, ee)
            total += cc * ca * inner
    return total


@functools.lru_cache(None)
def regdim(n: int, m: int, q: int, p1: int = 0) -> int:
    nn = n + p1
    for d in range(max(PRECISION, nn+1)):
        if reg_coeff(d, nn, m, q) <= 0:
            return d
    return math.inf


@functools.lru_cache(None)
def mq(n: int, m: int, q: int, p1: int = 0) -> int | float:
    d = regdim(n, m, q, p1)
    if not math.isfinite(d) or d >= 999:
        return math.inf
    return 3 * math.comb(n+p1+1, 2) * math.comb(n+p1-1+d, d)**2


@functools.lru_cache(None)
def hybrid_mq(n: int, m: int, q: int, p1: int = 0):
    if m < n:
        raise ValueError(f"HybridMQ requires m>=n, got {m}<{n}")
    best = (math.inf, None)
    for k in range(n):
        val = mq(n-k, m, q, p1) * q**k
        if val < best[0]:
            best = (val, k)
    return best


@functools.lru_cache(None)
def hashimoto_mq(n: int, m: int, q: int, p1: int = 0):
    if n < m:
        raise ValueError(f"HashimotoMQ requires n>=m, got {n}<{m}")
    best = (math.inf, None)
    for k in range(m):
        for a in range(2, m-k):
            if n < (a+1)*(m-k-a+1):
                continue
            if n < a*(m-k) - (a-1)**2 + k:
                continue
            val = q**k*(mq(m-a-k, m-a, q, p1)+mq(a-1, a-1, q, p1))
            val += (m-a-k+1)*mq(a, a, q, p1)
            if val < best[0]:
                best = (val, (a, k))
    return best


def solve_mq(n: int, m: int, q: int, p1: int = 0) -> dict:
    gate_per_field = 2*math.log2(q)**2 + math.log2(q)
    if n > m:
        h, ha = hashimoto_mq(n, m, q, p1)
        s, sa = hybrid_mq(m, m, q, p1)
        core, method, arg = (h, "Hashimoto", ha) if h <= s else (s, "square-specialization", sa)
    else:
        core, arg = hybrid_mq(n, m, q, p1)
        method = "HybridMQ"
    return {
        "n": n, "m": m, "q": q, "p1": p1,
        "degree_model": regdim(min(n,m), m, q, p1),
        "field_operations": core,
        "gate_cost": core * gate_per_field,
        "bits": math.log2(core) + math.log2(gate_per_field),
        "method": method, "arg": arg,
    }


PARAMS_23 = [
    ("I-square-l4",28,5,4,4,143,194), ("I-square-l2",48,16,2,2,143,164),
    ("I-rect-l4xr5",28,4,4,5,143,194), ("III-square-l4",40,7,4,4,207,264),
    ("III-square-l2",72,24,2,2,207,234), ("III-rect-l4xr5",38,5,4,5,207,239),
    ("V-square-l4",50,9,4,4,272,333), ("V-square-l2",96,32,2,2,272,303),
    ("V-rect-l4xr6",52,6,4,6,272,333),
]
PARAMS_24 = [
    ("I-l4",26,5,4,4,143), ("I-l2",43,17,2,2,143),
    ("III-l4",37,7,4,4,207), ("III-l2",64,25,2,2,207),
    ("V-l4xr6",44,6,4,6,272), ("V-l2",90,33,2,2,272),
]


def profile(name: str, v: int, o: int, l: int, r: int, target: int, spec: int | None = None) -> dict:
    n=v+o; m1=math.ceil(o*r/l); M=o*r*l
    old_m=min(M,m1*l*l)-(l-1); old_n=l*n-(l-1)
    new_m=min(M,m1*l*(l+1)//2); new_n=l*n-(M-new_m)
    return {
        "name":name,"v":v,"o":o,"q":19,"l":l,"r":r,"m1":m1,"M":M,
        "target_bits":target,"spec_table_forgery_bits":spec,
        "ordered":{"n":old_n,"m":old_m,"p1_0":solve_mq(old_n,old_m,19,0),"p1_1":solve_mq(old_n,old_m,19,1)},
        "symmetry_quotient":{"n":new_n,"m":new_m,"p1_0":solve_mq(new_n,new_m,19,0),"p1_1":solve_mq(new_n,new_m,19,1)},
    }


def main() -> None:
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default=str(root/'results'/'cost_profiles.json'))
    args=ap.parse_args()
    out={
        "version_2_3":[profile(*x) for x in PARAMS_23],
        "version_2_4_preview":[profile(*x) for x in PARAMS_24],
    }
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2,default=str)+'\n')
    print('Version 2.3 (headline p1=1)')
    for x in out['version_2_3']:
        q=x['symmetry_quotient']; b=q['p1_1']['bits']
        print(f"{x['name']:20s} {q['m']:3d}Q/{q['n']:3d}V  {b:8.2f}  margin {b-x['target_bits']:+7.2f}")
    print('\nVersion 2.4 preview (conditional, p1=1)')
    for x in out['version_2_4_preview']:
        q=x['symmetry_quotient']; b=q['p1_1']['bits']
        print(f"{x['name']:20s} {q['m']:3d}Q/{q['n']:3d}V  {b:8.2f}  margin {b-x['target_bits']:+7.2f}")


if __name__=='__main__':
    main()
