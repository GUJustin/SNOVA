#!/usr/bin/env python3
import json, math
from fractions import Fraction
from pathlib import Path
q=19; Q=q*q; GM=692; GA=84; GROOT=16886
m,k,p=15,4,7; ell=m-k-p; zdim=9
assert ell==4 and zdim+k==13

def aa(x): return x*(x+1)//2+x, x*(x+1)//2+x-1
def av(x): return x,x-1
def sub(x,nlin,check=False):
    a,b=aa(x); a*=nlin; b*=nlin
    if not check:
        c,d=av(x); a+=c*nlin; b+=d*nlin
    return a,b
# Exact per-live-node arithmetic components from the published MOR ledger.
qm=qa=0
for i in range(p):
    a,b=sub(k+i,1); qm+=a; qa+=b
quadratic=qm*GM+qa*GA+p*GROOT
child=p*(2**14)
a,b=sub(k+p,ell)
residual=(ell*a)*GM+(ell*b)*GA
# Keep the old 5x5 singular-safe reserve unchanged; 4x4 cannot require more.
singular=397_312
completion=19_400
# Four selected withheld equations plus the one omitted original D-form.
cm,ca=aa(m)
one_check=cm*GM+ca*GA
checks=5*one_check
control=1_000_000
components={
 'seven_quadratic_substitutions_and_root_calls':quadratic,
 'child_generation_reserve':child,
 'four_residual_affine_equations':residual,
 'singular_safe_elimination_kernel_reserve':singular,
 'completion_generation':completion,
 'five_full_diagonal_checks':checks,
 'control_reserve':control,
}
C=sum(components.values())
assert C==3_500_734 and C<2**22
# Full-diagonal rank-17 certificate abundance.
p_abund=1-Fraction(1,Q**16)-4*(Fraction(1,Q)-Fraction(1,Q**17))
rho=Fraction(49,829); c=Fraction(1,4)
p_cross=c/(1+c*(q*rho)**16)
p_trial=p_abund*p_cross
trial_actual=Q**13*C
trial_clean=Q**13*2**22
actual_exp=math.log2(trial_actual)-math.log2(float(p_trial))
clean_exp=13*math.log2(Q)+22-math.log2(float(p_trial))
out={
 'parameters':{'q':q,'A_size':Q,'m_selected':m,'k':k,'p':p,'linear_unknowns':ell,'streamed_A_dimension':zdim,'outer_A_exponent':13},
 'gate_costs':{'A_multiplication_AXN':GM,'A_addition_AXN':GA,'quadratic_root_AXN':GROOT},
 'components_per_expected_A13_aggregate':components,
 'component_total_AXN':C,
 'component_total_log2':math.log2(C),
 'clean_component_bound':'<2^22',
 'clean_bound_margin_AXN':2**22-C,
 'abundance':{
   'mean_final_D_roots':'361^8 = 19^16',
   'rank_certificate':'rank_A K_full(lambda) >= 17 for every nonzero lambda in A^16',
   'failure_upper_bound':float(1-p_abund),
   'success_lower_bound':float(p_abund)},
 'cross_channel':{'rho_cond':[49,829],'c':[1,4],'success_lower_bound':float(p_cross)},
 'combined_trial_success_lower_bound':float(p_trial),
 'restart_bits_upper_bound':-math.log2(float(p_trial)),
 'actual_component_sum':{'trial_exponent':math.log2(trial_actual),'success_adjusted_exponent':actual_exp,'headroom_to_143':143-actual_exp},
 'clean_2^22_component_bound':{'trial_exponent':13*math.log2(Q)+22,'success_adjusted_exponent':clean_exp,'headroom_to_143':143-clean_exp},
 'status':'fixed-key rank certificate + target-averaged aggregate Boolean ledger; key-distribution success is not claimed'
}
Path('/mnt/data/lib/jg_m15_k4_p7_audited_ledger.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps(out,indent=2,sort_keys=True))
