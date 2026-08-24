#!/usr/bin/env python3
"""Deterministic formula, circuit, and ledger checks for the SNOVA artifact.

This checker runs the separately coded all-nine recomputation, exhaustively
checks the F_19^2 quadratic-root routine on every discriminant, and recomputes
the Level-I Just Guess finite-gate and cross-channel numbers.  These are
internal consistency checks, not independent cryptanalytic evidence.
"""
from __future__ import annotations

import json
import math
import runpy
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = 19

# All-nine conditional-ledger checks.
runpy.run_path(str(HERE / "verify_ledgers.py"), run_name="__main__")

# Route-specific Level-I homotopy necessary-cost audit.
runpy.run_path(
    str(HERE / "verify_level1_homotopy_operation_counter.py"),
    run_name="__main__",
)

# Final ell=2 and ell=4 five-component PXL route ledgers.
runpy.run_path(str(HERE / "pxl_final_ledger.py"), run_name="__main__")

# ``verify_ledgers.py`` evaluates the committed scalar netlist on all
# 361 valid inputs and recomputes the tower identities and recurrence counts.
# Regeneration is deliberately separate so that verification is nonmutating.

# F_19^2 arithmetic, u^2=-1.
def add(x, y): return ((x[0]+y[0]) % P, (x[1]+y[1]) % P)
def sub(x, y): return ((x[0]-y[0]) % P, (x[1]-y[1]) % P)
def neg(x): return ((-x[0]) % P, (-x[1]) % P)
def mul(x, y): return ((x[0]*y[0]-x[1]*y[1]) % P,
                       (x[0]*y[1]+x[1]*y[0]) % P)
def powe(x, e):
    z=(1,0)
    while e:
        if e & 1: z=mul(z,x)
        x=mul(x,x); e//=2
    return z

ZERO=(0,0); ONE=(1,0); MINUS_ONE=(18,0)
elems=[(a,b) for a in range(P) for b in range(P)]
def legendre(x):
    if x == ZERO: return 0
    y=powe(x,180)
    return 1 if y==ONE else -1 if y==MINUS_ONE else None
NONRES=next(x for x in elems if legendre(x)==-1)
C0=powe(NONRES,45)

class Counter:
    def __init__(self): self.m=0; self.eq=0
    def M(self,x,y): self.m+=1; return mul(x,y)
    def E(self,x,y): self.eq+=1; return x==y

def sqrt_ts(n):
    c=Counter()
    if c.E(n,ZERO): return [ZERO],c
    n2=c.M(n,n); n4=c.M(n2,n2); n8=c.M(n4,n4)
    n16=c.M(n8,n8); n32=c.M(n16,n16)
    t=c.M(c.M(c.M(n32,n8),n4),n)  # n^45
    r=c.M(c.M(c.M(n16,n4),n2),n)  # n^23
    t2=c.M(t,t); leg=c.M(t2,t2)
    if not c.E(leg,ONE): return [],c
    M=3; cc=C0
    while not c.E(t,ONE):
        tt=c.M(t,t); i=1
        while i<M and not c.E(tt,ONE):
            tt=c.M(tt,tt); i+=1
        assert i<M
        b=cc
        for _ in range(M-i-1): b=c.M(b,b)
        r=c.M(r,b); b2=c.M(b,b); t=c.M(t,b2); cc=b2; M=i
    rr=neg(r)
    return ([r] if rr==r else [r,rr]),c

stats=[]
for d in elems:
    roots,c=sqrt_ts(d)
    truth=[x for x in elems if mul(x,x)==d]
    assert set(roots)==set(truth)
    stats.append((c.m,c.eq))
assert max(m for m,_ in stats)==22
assert max(e for _,e in stats)==8

const=json.loads((HERE/'f19_constant_multiplier_ledger.json').read_text())['counts']
G_MUL,G_ADD,G_SUB,G_NEG,G_EQ=692,84,116,32,19
g4=2*int(const['4']); ghalf=2*int(const['10'])
maxsqrt=max(m*G_MUL+e*G_EQ for m,e in stats)
quad=G_MUL+g4+G_SUB+maxsqrt+G_NEG+G_ADD+G_SUB+2*ghalf
assert quad==16886
root_ledger=json.loads((HERE/'f361_quadratic_root_ledger.json').read_text())
assert root_ledger['complete_monic_quadratic_root_upper_bound_AXN']==quad

# Exact re-enumeration of the conditional ell=2 channel atoms.
weights=[14 if x<=8 else 13 for x in range(P)]
diag=defaultdict(lambda: defaultdict(int))
for a in range(P):
    for b in range(P):
        for c in range(P):
            D=((a+7*b+16*c)%P,(2*b+7*c)%P)
            H=(a+7*b+18*c)%P
            diag[D][H]+=weights[a]*weights[b]*weights[c]
diag_max=Fraction(0)
for dist in diag.values():
    den=sum(dist.values())
    diag_max=max(diag_max,*(Fraction(num,den) for num in dist.values()))
assert diag_max==Fraction(49,829)

off=defaultdict(lambda: defaultdict(int))
for a in range(P):
    for b in range(P):
        for c in range(P):
            for d in range(P):
                D=((a+13*b+13*c+16*d)%P,(b+c+7*d)%P)
                H=((a+13*b+13*c+18*d)%P,(18*b+c)%P)
                off[D][H]+=weights[a]*weights[b]*weights[c]*weights[d]
functionals=[(1,beta) for beta in range(P)]+[(0,1)]
off_max=Fraction(0)
for dist in off.values():
    den=sum(dist.values())
    for alpha,beta in functionals:
        vals=defaultdict(int)
        for (h0,h1),num in dist.items():
            vals[(alpha*h0+beta*h1)%P]+=num
        off_max=max(off_max,*(Fraction(num,den) for num in vals.values()))
assert off_max==Fraction(169246,2971565)
assert off_max<diag_max

# Separate recomputation of the Just Guess ledger.
Q=P*P; m,n,k,p=16,64,5,6; ell=m-k-p
def aa(x): return x*(x+1)/2+x, x*(x+1)/2+x-1
def av(x): return x,x-1
def subst(x,nlin,check=False):
    a,b=aa(x); a*=nlin; b*=nlin
    if not check:
        c,d=av(x); a+=c*nlin; b+=d*nlin
    return a,b
fm=fa=0.0
for i in range(p):
    a,b=subst(k+i,1); fm+=a; fa+=b
a,b=subst(k+p,ell); fm+=ell*a; fa+=ell*b
fm+=ell**3/3; fa+=ell**3/3
a,b=aa(m); factor=Q*(1-Q**(-k))/(Q-1)
fm+=a*factor; fa+=b*factor
per_guess=math.ceil(fm)*G_MUL+math.ceil(fa)*G_ADD+p*quad+1_000_000
assert per_guess==3_179_584
base=Q**8*(Q**k*per_guess+2**40)
base_exp=math.log2(base)
assert abs(base_exp-132.04652202347475)<1e-12
candidate_fraction=Fraction(1,4)
beta=(P*diag_max)**16
p_cross=float(candidate_fraction/(1+candidate_fraction*beta))
assert abs(p_cross-0.09613444074026915)<1e-15
cross_bits=math.log2(1/p_cross)
adjusted=base_exp+cross_bits
assert abs(adjusted-135.4253248354897)<1e-12
assert abs((143-adjusted)-7.574675164510299)<1e-12
assert abs((adjusted+4)-139.4253248354897)<1e-12

# Strictly capped binary-tree sensitivity ledger.
capped=json.loads((HERE/'jg_level1_capped_tree_ledger.json').read_text())
assert capped['tree']['quadratic_nodes_per_guess']==63
assert capped['tree']['leaves_per_guess']==64
assert capped['gate_costs']['per_guess']==161548042
assert abs(capped['per_trial']['exponent']-137.71350337514204)<1e-12
assert abs(capped['success_adjusted_before_kappa_JG_cap']['exponent']-141.092306187157)<1e-12
assert abs(capped['break_even']['log2_kappa_JG_cap_upper_bound']-1.9076938128430099)<1e-12

atom_ledger=json.loads((HERE/'l2_channel_conditional_atom_ledger.json').read_text())
assert Fraction(*atom_ledger['diagonal_max_conditional_atom']['fraction'])==diag_max
assert Fraction(*atom_ledger['off_diagonal_max_conditional_functional_atom']['fraction'])==off_max
ledger=json.loads((HERE/'jg_level1_expected_ledger.json').read_text())
assert abs(ledger['per_trial']['exponent']-base_exp)<1e-12
assert abs(ledger['cross_second_moment']['success_probability_lower_bound']-p_cross)<1e-15
assert abs(ledger['success_adjusted_before_kappa_JG']['exponent']-adjusted)<1e-12

print("Artifact consistency checks passed")
print("- all-nine conditional-ledger recomputation")
print("- F19^2 Tonelli-Shanks on all 361 discriminants")
print("- 16,886-gate monic-quadratic all-roots bound")
print("- exact ell=2 conditional atoms 49/829 and 169246/2971565")
print("- 2^132.046522 expected-tree and 2^137.713503 capped-tree Just Guess ledgers")
print("- cross success > 0.096134 and exponent 135.425325 + log2(kappa_JG)")
print("- Level-I break-even log2(kappa_JG) < 7.574675")
print("- final ell=2 and ell=4 five-component PXL route ledgers")
