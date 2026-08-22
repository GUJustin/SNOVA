#!/usr/bin/env python3
from collections import defaultdict
from fractions import Fraction
import math, json
q=19
weights=[14]*9+[13]*10
rho=Fraction(49,829)
# diagonal D in F19^2, H in F19
den=defaultdict(int); num=defaultdict(int)
for a in range(q):
 for b in range(q):
  for c in range(q):
   ww=weights[a]*weights[b]*weights[c]
   D=((a+7*b+16*c)%q,(2*b+7*c)%q)
   H=(a+7*b+18*c)%q
   den[D]+=ww; num[(D,H)]+=ww
diag=max(Fraction(n,den[D]) for (D,H),n in num.items())
assert diag==rho
# off-diagonal D,H in F19^2
Den=defaultdict(int); Joint=defaultdict(int)
rows=[]
for a in range(q):
 for b in range(q):
  for c in range(q):
   for d in range(q):
    ww=weights[a]*weights[b]*weights[c]*weights[d]
    D=((a+13*b+13*c+16*d)%q,(b+c+7*d)%q)
    H=((a+13*b+13*c+18*d)%q,(18*b+c)%q)
    Den[D]+=ww; Joint[(D,H)]+=ww
sigma=max(Fraction(n,Den[D]) for (D,H),n in Joint.items())
assert sigma==Fraction(38416,11886083)
assert sigma < rho*rho
# scalar projective functionals (1,t), plus (0,1)
forms=[(1,t) for t in range(q)]+[(0,1)]
scalar_max=Fraction(0,1); scalar_arg=None
for u,v in forms:
 S=defaultdict(int)
 for (D,H),n in Joint.items():
  z=(u*H[0]+v*H[1])%q
  S[(D,z)]+=n
 m=max((Fraction(n,Den[D]),D,z) for (D,z),n in S.items())
 if m[0]>scalar_max: scalar_max=m[0]; scalar_arg=(u,v,m[1],m[2])
assert scalar_max==Fraction(169246,2971565)
assert scalar_max < rho
# projective union bound
P47=(q**48-1)//(q-1)
log2_bad=2*math.log2(P47)+126*math.log2(float(rho))
assert log2_bad < -114.684
# work ledger
Q=q*q
clean=Q**13*2**21+Q**8*2**40+Q**9*2**40
literal=1947153*Q**13+Q**8*2**40+Q**9*2**40
p_good_frac=Fraction(q**48,2*q**48-1)
p_good=float(p_good_frac)
clean_exp=math.log2(clean)-math.log2(p_good)
literal_exp=math.log2(literal)-math.log2(p_good)
assert clean_exp < 132.447
assert literal_exp < 132.340
out={
 'q':q,
 'diag_conditional_atom':str(diag),
 'offdiag_scalar_functional_atom':str(scalar_max),
 'offdiag_joint_point_atom':str(sigma),
 'rho_squared':str(rho*rho),
 'P47':str(P47),
 'coefficient_rank':126,
 'log2_bad_joint_pencil_union_bound':log2_bad,
 'good_key_trial_success_lower_bound_exact':str(p_good_frac),
 'clean_success_adjusted_log2':clean_exp,
 'literal_success_adjusted_log2':literal_exp,
}
print(json.dumps(out,indent=2))
