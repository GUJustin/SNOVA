#!/usr/bin/env python3
"""Regenerate the strictly capped Level-I Just Guess tree sensitivity ledger.

The capped schedule explores at most every binary root branch for p=6, hence
63 quadratic nodes and 64 leaves per guess.  A branch is abandoned if the
residual 5x5 linear system is not unique.  This makes the work per trial a
finite upper bound; candidate abundance among surviving branches remains the
separate event measured by kappa_JG^cap.
"""
from __future__ import annotations
import json, math
from pathlib import Path

HERE=Path(__file__).resolve().parent
ART=HERE.parent
Q=361; m=16; k=5; p=6; ell=m-k-p; zdim=8
GM,GA=692,84
GR=int(json.loads((ART/'f361_quadratic_root_ledger.json').read_text())['complete_monic_quadratic_root_upper_bound_AXN'])

def aa(n): return n*(n+1)//2+n, n*(n+1)//2+n-1
def av(n): return n,n-1
def sub(n,nlin,check=False):
    a,b=aa(n); a*=nlin; b*=nlin
    if not check:
        c,d=av(n); a+=c*nlin; b+=d*nlin
    return a,b

mul=add=root_calls=0
for i in range(p):
    nodes=2**i
    a,b=sub(k+i,1)
    mul+=nodes*a; add+=nodes*b; root_calls+=nodes
leaves=2**p
lsm,lsa=sub(k+p,ell); lsm*=ell; lsa*=ell
# Conservative fixed 5x5 Gauss-Jordan bound, including <=16 multiplies/inverse.
gauss_mul=16*ell + ell*ell + ell*(ell-1)*ell
gauss_add=ell*(ell-1)*ell
cm,ca=aa(m); cm*=k; ca*=k
mul+=leaves*(lsm+gauss_mul+cm)
add+=leaves*(lsa+gauss_add+ca)
control=1_000_000
per_guess=mul*GM+add*GA+root_calls*GR+control
filter_per_z=2**40
trial=Q**zdim*(Q**k*per_guess+filter_per_z)
trial_exp=math.log2(trial)
# Same c=1/4 cross-channel second moment as the expected-tree ledger.
rho=14/256
p_cross=.25/(1+.25*(19*rho)**16)
cross_bits=math.log2(1/p_cross)
adjusted=trial_exp+cross_bits
out={
  'model':'strictly capped binary-tree schedule; rank-deficient residual linear branches are abandoned',
  'parameters':{'A_size':Q,'m':m,'k':k,'p':p,'linear_unknowns':ell,'z_A_dimension':zdim},
  'tree':{'quadratic_nodes_per_guess':root_calls,'leaves_per_guess':leaves,'field_multiplications_per_guess':mul,'field_additions_per_guess':add},
  'gate_costs':{'F19_2_multiplication':GM,'F19_2_addition':GA,'quadratic_all_roots':GR,'control_per_guess':control,'filter_and_verifier_per_z':filter_per_z,'per_guess':per_guess},
  'per_trial':{'AXN':trial,'exponent':trial_exp},
  'cross_second_moment':{'success_probability_lower_bound':p_cross,'restart_bits_upper_bound':cross_bits},
  'success_adjusted_before_kappa_JG_cap':{'exponent':adjusted,'headroom_to_143':143-adjusted},
  'break_even':{'log2_kappa_JG_cap_upper_bound':143-adjusted},
}
(ART/'jg_level1_capped_tree_ledger.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
