#!/usr/bin/env python3
"""Regenerate the finite-gate Level-I Just Guess sensitivity ledger."""
from __future__ import annotations
import json, math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ART = HERE.parent
q = 19
Q = q*q
m, n, k, p = 16, 64, 5, 6
ell = m-k-p
G_MUL, G_ADD = 692, 84
G_ROOT = int(json.loads((ART/'f361_quadratic_root_ledger.json').read_text())['complete_monic_quadratic_root_upper_bound_AXN'])

def aa(x): return x*(x+1)/2+x, x*(x+1)/2+x-1
def av(x): return x, x-1
def sub(x, nlin, check=False):
    a,b=aa(x); a*=nlin; b*=nlin
    if not check:
        c,d=av(x); a+=c*nlin; b+=d*nlin
    return a,b

mul=add=0.0
for i in range(p):
    a,b=sub(k+i,1); mul+=a; add+=b
a,b=sub(k+p,ell); mul+=ell*a; add+=ell*b
mul+=ell**3/3; add+=ell**3/3
a,b=aa(m)
fac=Q*(1-Q**(-k))/(Q-1)
mul+=a*fac; add+=b*fac
control=1_000_000
per_guess=math.ceil(mul)*G_MUL+math.ceil(add)*G_ADD+p*G_ROOT+control
per_z_filter=2**40
base=Q**8*(Q**k*per_guess+per_z_filter)
base_exp=math.log2(base)
rho=Fraction(14,256)
c=Fraction(1,4)
beta=float((q*rho)**16)
p_cross=float(c)/(1+float(c)*beta)
cross_mult=1/p_cross
cross_bits=math.log2(cross_mult)
adjusted=base_exp+cross_bits
out={
 'model':'published Just Guess expected-operation model with a fully instantiated per-trial Boolean circuit and an exact cross-channel second moment',
 'parameters':{'q':q,'A_size':Q,'n':n,'m':m,'k':k,'p':p,'linear_unknowns':ell,'z_A_dimension':8},
 'expected_field_counts_per_guess':{'multiplications':mul,'additions':add,'quadratic_root_calls':p},
 'gate_costs':{'F19_2_multiplication':G_MUL,'F19_2_addition':G_ADD,'quadratic_all_roots':G_ROOT,'control_per_guess':control,'per_guess':per_guess,'filter_and_verifier_per_z':per_z_filter},
 'per_trial':{'AXN':base,'exponent':base_exp},
 'candidate_event':{'minimum_sign_canonical_candidates':'19^16/4','regularity_multiplier_symbol':'kappa_JG'},
 'cross_second_moment':{'beta':beta,'success_probability_lower_bound':p_cross,'restart_multiplier_upper_bound':cross_mult,'restart_bits_upper_bound':cross_bits},
 'success_adjusted_before_kappa_JG':{'exponent':adjusted,'headroom_to_143':143-adjusted},
 'break_even':{'log2_kappa_JG_upper_bound':143-adjusted},
 'sensitivity':{}
}
for r in [2,4,8,16,32,64,128,256]:
    e=adjusted+math.log2(r)
    out['sensitivity'][str(r)]={'exponent':e,'headroom_to_143':143-e}
(ART/'jg_level1_revised_ledger.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
