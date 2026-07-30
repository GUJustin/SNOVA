#!/usr/bin/env python3
"""Exact Frobenius-orbit-complete homotopy ledgers for q=19, ell=4 SNOVA.

For the near-square A-linear slices h=K-2, every full-rank Jacobian basis has
one of 55 family-count patterns.  At an F_19-rational (descent) root, Frobenius
cyclically permutes the four eigenblocks and carries a nonzero minor to a
nonzero conjugate minor.  It therefore suffices to solve one representative
from each C4 orbit of patterns.  This script enumerates the patterns and
orbits, computes exact multihomogeneous Bezout/companion quantities, applies
the finite-cardinality separator threshold |k| >= 8 B^2, and prices all
filters and retries.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
from typing import Iterable

getcontext().prec = 120
ROOT = Path(__file__).resolve().parent
Q = 19
D = 4
A_SIZE = Q**D
RHO = Fraction(14, 256)
G4_AXN = 6081
PAIRS = ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))
FAMILIES = ((0,0),(1,1),(2,2),(3,3)) + PAIRS
FAMILY_INDEX = {tuple(sorted(p)): i for i,p in enumerate(FAMILIES)}
ROTATION = tuple(FAMILY_INDEX[tuple(sorted(((a+1)%4,(b+1)%4)))] for a,b in FAMILIES)


def dlog2_fraction(x: Fraction) -> Decimal:
    if x <= 0:
        raise ValueError("logarithm of nonpositive rational")
    return (Decimal(x.numerator).ln() - Decimal(x.denominator).ln()) / Decimal(2).ln()


def dlog2_int(x: int) -> Decimal:
    return Decimal(x).ln() / Decimal(2).ln()


def human_bytes(n: int) -> str:
    units = [(2**80,"YiB"),(2**70,"ZiB"),(2**60,"EiB"),(2**50,"PiB"),
             (2**40,"TiB"),(2**30,"GiB"),(2**20,"MiB"),(2**10,"KiB")]
    for scale,name in units:
        if n >= scale:
            return f"{Decimal(n)/Decimal(scale):.3f} {name}"
    return f"{n} B"


def bounded_compositions(total: int, length: int, cap: int) -> Iterable[tuple[int,...]]:
    vals=[0]*length
    def rec(i: int, rem: int):
        if i==length-1:
            if 0 <= rem <= cap:
                vals[i]=rem
                yield tuple(vals)
            return
        lo=max(0, rem-cap*(length-i-1)); hi=min(cap,rem)
        for x in range(lo,hi+1):
            vals[i]=x
            yield from rec(i+1,rem-x)
    yield from rec(0,total)


def rotate(v: tuple[int,...]) -> tuple[int,...]:
    out=[0]*10
    for i,j in enumerate(ROTATION):
        out[j]=v[i]
    return tuple(out)


def orbit(v: tuple[int,...]) -> tuple[tuple[int,...],...]:
    out=[]; x=v
    for _ in range(4):
        if x not in out:
            out.append(x)
        x=rotate(x)
    return tuple(out)


def orient(edges: tuple[int,...], target: tuple[int,...]) -> int:
    """Coefficient of prod_{ab}(t_a+t_b)^e_ab at target."""
    if sum(edges) != sum(target) or min(target) < 0:
        return 0
    e01,e02,e03,e12,e13,e23=edges
    r0,r1,r2,r3=target
    ans=0
    for k01 in range(e01+1):
        for k02 in range(e02+1):
            k03=r0-k01-k02
            if not (0 <= k03 <= e03):
                continue
            for k12 in range(e12+1):
                k13=r1-(e01-k01)-k12
                if not (0 <= k13 <= e13):
                    continue
                k23=r2-(e02-k02)-(e12-k12)
                if not (0 <= k23 <= e23):
                    continue
                if (e03-k03)+(e13-k13)+(e23-k23) != r3:
                    continue
                ans += (math.comb(e01,k01)*math.comb(e02,k02)*math.comb(e03,k03)
                        *math.comb(e12,k12)*math.comb(e13,k13)*math.comb(e23,k23))
    return ans


def bezout_companion(counts: tuple[int,...], s: int) -> tuple[int,int]:
    loops=counts[:4]; edges=counts[4:]
    target=tuple(s-loops[i] for i in range(4))
    B=(2**sum(loops))*orient(edges,target)
    if B==0:
        return 0,0
    Bp=B
    # Omit one diagonal row.  Its theta_0 contribution has coefficient sum
    # over the four cap-minus-one monomials.
    for a,mult in enumerate(loops):
        if not mult:
            continue
        subtotal=0
        for j in range(4):
            rr=list(target); rr[a]+=1; rr[j]-=1
            subtotal += orient(edges,tuple(rr))
        Bp += mult*(2**(sum(loops)-1))*subtotal
    # Omit one cross row.
    for z,mult in enumerate(edges):
        if not mult:
            continue
        ee=list(edges); ee[z]-=1
        subtotal=0
        for j in range(4):
            rr=list(target); rr[j]-=1
            subtotal += orient(tuple(ee),tuple(rr))
        Bp += mult*(2**sum(loops))*subtotal
    return B,Bp


def extension_degree(threshold: int) -> int:
    r=1; size=A_SIZE
    while size < threshold:
        size *= A_SIZE; r += 1
    return r


def packed_bytes_per_element(r: int) -> int:
    return math.ceil(D*r*math.log2(Q)/8)


@dataclass(frozen=True)
class Level:
    name: str
    reference: int
    m1: int
    K: int
    n_rejection: int
    t_by_shape: tuple[int,int]

LEVELS=(
    Level("I",143,5,50,33,(74,74)),
    Level("III",207,7,70,47,(110,114)),
    Level("V",272,9,90,59,(138,146)),
)


def pattern_record(level: Level, s: int, deficit: tuple[int,...]) -> dict[str,object]:
    m=level.m1; counts=tuple(m-x for x in deficit)
    loops=counts[:4]; edges=counts[4:]
    B,Bp=bezout_companion(counts,s)
    if not B:
        raise ValueError("zero Bezout number in near-square pattern")
    h=4*s
    qdiag=s*(s+1)//2+s+1
    qcross=(s+1)**2
    active_diag=sum(x>0 for x in loops); active_cross=sum(x>0 for x in edges)
    pre=active_diag*s*(s+1)//2+active_cross*s*s
    L=pre+2*sum(loops)*qdiag+2*sum(edges)*qcross
    H=L+2*h+h*h
    group_sums=[2*loops[i] for i in range(4)]
    for c,(a,b) in zip(edges,PAIRS):
        group_sums[a]+=c; group_sums[b]+=c
    # New finite-cardinality form: arbitrary distinct start nodes and a
    # uniformly random coefficient vector as separator.  8 B^2 gives <=1/8
    # separator failure; 8h is included for the projection determinant.
    threshold=max(max(group_sums),8*B*B,8*h)
    r=extension_degree(threshold); c_r=2*r-1
    pre_full=4*s*(s+1)//2+6*s*s
    L_full=pre_full+2*m*(4*qdiag+6*qcross)
    hom=c_r*B*Bp*H
    filt=c_r*B*L_full
    elem_bytes=packed_bytes_per_element(r)
    output_elements=(h+2)*B
    return {
        "deficit":list(deficit),"counts":list(counts),"loops":list(loops),
        "edges":list(edges),"B":B,"B_plus":Bp,"L":L,"H":H,
        "group_degree_sums":group_sums,"threshold_E":threshold,
        "extension_degree_r":r,"pointwise_factor":c_r,
        "homotopy_weight":hom,"filter_weight":filt,
        "output_elements":output_elements,"packed_bytes_per_element":elem_bytes,
        "packed_output_bytes":output_elements*elem_bytes,
    }


def evaluate_level(level: Level) -> dict[str,object]:
    # K=10m is 2 mod 4 for odd m; nearest A-linear square-under slice is K-2.
    s=(level.K-2)//4; h=4*s; tau=level.K-h
    assert tau==2
    patterns=[]
    for deficit in bounded_compositions(tau,10,level.m1):
        patterns.append(pattern_record(level,s,deficit))
    by_def={tuple(p["deficit"]):p for p in patterns}
    seen=set(); orbit_rows=[]
    for deficit in sorted(by_def):
        if deficit in seen:
            continue
        o=orbit(deficit); seen.update(o)
        rep=min(o)
        row=dict(by_def[rep]); row["orbit_size"]=len(o); row["orbit_members"]=[list(x) for x in o]
        # Exact invariance checks.
        for x in o:
            other=by_def[x]
            for key in ("B","B_plus","L","H","extension_degree_r","homotopy_weight","filter_weight"):
                if other[key] != row[key]:
                    raise AssertionError((level.name,deficit,key,row[key],other[key]))
        orbit_rows.append(row)
    if len(patterns)!=55 or len(orbit_rows)!=16:
        raise AssertionError((len(patterns),len(orbit_rows)))

    sum_hom_full=sum(int(p["homotopy_weight"]) for p in patterns)
    sum_filter_full=sum(int(p["filter_weight"]) for p in patterns)
    sum_hom=sum(int(p["homotopy_weight"]) for p in orbit_rows)
    sum_filter=sum(int(p["filter_weight"]) for p in orbit_rows)

    # One orbit-complete sweep, including repetitions for projection and
    # homotopy randomness.  Use the conservative base-field lower bound
    # 1-h/|A| although projections are sampled in the larger extension.
    projection=Fraction(A_SIZE-h,A_SIZE)
    per_sweep=Fraction(8*G4_AXN,7)*Fraction(1,1)/projection*(h*sum_hom+sum_filter)
    per_full=Fraction(8*G4_AXN,7)*Fraction(1,1)/projection*(h*sum_hom_full+sum_filter_full)
    alpha=Fraction((A_SIZE-1)**level.n_rejection,A_SIZE**level.n_rejection)
    costs={}
    for label,eta in (("compact",Fraction(1,Q)),("robust",Fraction(1,2))):
        root_term=alpha-eta/Fraction(Q-1,1)
        if root_term<=0: raise AssertionError("bad root term")
        mu = Fraction(Q**(h-level.K),1) if h >= level.K else Fraction(1,Q**(level.K-h))
        success=mu*root_term*root_term/(1+eta)
        total=per_sweep/success
        full_total=per_full/success
        bits=dlog2_fraction(total)
        costs[label]={
            "eta_numerator":eta.numerator,"eta_denominator":eta.denominator,
            "root_success_numerator":success.numerator,"root_success_denominator":success.denominator,
            "bits":float(bits),"margin":float(Decimal(level.reference)-bits),
            "all_55_bits":float(dlog2_fraction(full_total)),
            "frobenius_gain_bits":float(dlog2_fraction(full_total/total)),
        }

    baseline=Fraction(Q**h-1,Q**level.K)
    tails=[]
    for t in level.t_by_shape:
        row={"t":t}
        for label,eta in (("compact",Fraction(1,Q)),("robust",Fraction(1,2))):
            eps=eta-baseline
            if eps<=0: raise AssertionError("threshold below baseline")
            bound=Fraction(Q**h-1,1)*(1-Fraction(1,Q**level.K))*RHO**t/eps
            row[label]={"numerator":bound.numerator,"denominator":bound.denominator,
                        "log2_bound":float(dlog2_fraction(bound))}
        tails.append(row)

    max_output=max(int(p["packed_output_bytes"]) for p in orbit_rows)
    total_output=sum(int(p["packed_output_bytes"]) for p in orbit_rows)
    return {
        "level":level.name,"reference":level.reference,"m1":level.m1,"K":level.K,
        "s":s,"h":h,"tau":tau,"all_pattern_count":len(patterns),
        "frobenius_orbit_count":len(orbit_rows),
        "orbit_sizes":sorted([int(p["orbit_size"]) for p in orbit_rows]),
        "rotation_on_families":list(ROTATION),
        "sum_homotopy_weight_all_55":sum_hom_full,
        "sum_homotopy_weight_orbits":sum_hom,
        "sum_filter_weight_all_55":sum_filter_full,
        "sum_filter_weight_orbits":sum_filter,
        "per_orbit_sweep_AXN_numerator":per_sweep.numerator,
        "per_orbit_sweep_AXN_denominator":per_sweep.denominator,
        "per_orbit_sweep_bits":float(dlog2_fraction(per_sweep)),
        "costs":costs,"spectrum_tails":tails,
        "peak_individual_output_bytes":max_output,
        "peak_individual_output_human":human_bytes(max_output),
        "sum_sequential_output_bytes":total_output,
        "sum_sequential_output_human":human_bytes(total_output),
        "orbits":orbit_rows,
    }


def write_markdown(data: dict[str,object]) -> None:
    lines=[
      "# Frobenius-orbit-complete homotopy certificate","",
      "For h=K-2 there are 55 family-count patterns.  At a descent root, q-Frobenius cyclically rotates the four eigenblocks and carries every nonzero Jacobian minor to a nonzero conjugate minor.  One representative from each C4 orbit therefore suffices.  The exact orbit count is 16.","",
      "The finite-cardinality threshold uses arbitrary distinct start nodes and a uniformly random coefficient-vector separator: |k| >= max(max_j sum_i d_ij, 8 B^2, 8h).  All omitted equations and descent/rejection/verifier checks are charged as filters.","",
      "| Level | (h,K) | patterns -> orbits | compact AXN | compact margin | robust AXN | robust margin | gain over 55 | peak one-core output |",
      "|:--:|:--:|:--:|--:|--:|--:|--:|--:|--:|",
    ]
    for z in data["levels"]:
        c=z["costs"]["compact"]; r=z["costs"]["robust"]
        lines.append(f"| {z['level']} | ({z['h']},{z['K']}) | {z['all_pattern_count']} -> {z['frobenius_orbit_count']} | {c['bits']:.5f} | {c['margin']:.5f} | {r['bits']:.5f} | {r['margin']:.5f} | {c['frobenius_gain_bits']:.5f} | {z['peak_individual_output_human']} |")
    lines += ["","## Random-XOF excess-spectrum tails","",
      "| Level | t by shape | log2 bad compact | log2 bad robust |",
      "|:--:|:--:|:--:|:--:|"]
    for z in data["levels"]:
        ts="/".join(str(x["t"]) for x in z["spectrum_tails"])
        cs="/".join(f"{x['compact']['log2_bound']:.2f}" for x in z["spectrum_tails"])
        rs="/".join(f"{x['robust']['log2_bound']:.2f}" for x in z["spectrum_tails"])
        lines.append(f"| {z['level']} | {ts} | {cs} | {rs} |")
    lines += ["","The packed-output values are straightforward dense parametrization ceilings, not peak-memory theorems.  Sequential processing needs only the largest individual parametrization plus solver workspace, but these ceilings are still very large.",""]
    (ROOT/"FROBENIUS_ORBIT_SWEEP_CERTIFICATE.md").write_text("\n".join(lines))


def main() -> None:
    levels=[evaluate_level(x) for x in LEVELS]
    data={
      "field":{"q":Q,"A_size":A_SIZE,"rho":"14/256","g4_AXN":G4_AXN},
      "theorem_scope":"near-square A-linear slices h=K-2 for ell=4",
      "levels":levels,
    }
    (ROOT/"frobenius_orbit_sweep_certificate.json").write_text(json.dumps(data,indent=2)+"\n")
    write_markdown(data)
    for z in levels:
        print(z["level"],z["costs"]["compact"]["bits"],z["costs"]["robust"]["bits"],z["costs"]["compact"]["frobenius_gain_bits"])

if __name__=="__main__":
    main()
