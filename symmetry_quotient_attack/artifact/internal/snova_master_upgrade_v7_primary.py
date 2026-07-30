#!/usr/bin/env python3
"""Defensible v7 ledger for the SNOVA symmetry-quotient upgrade.

Primary ledger inherited from v6, with four strengthening steps:
  * full-domain accepted-root spectrum, with excess 2^-128;
  * the sharp second-moment denominator a+eta rather than 1+eta;
  * the finite-cardinality homotopy theorem's conservative B^2 separator bound;
  * zero-offset chosen-message ell=2 route (no reserved-line key-density charge);
  * explicit 150/692/2628-gate F19/F19^2/F19^4 multiplication circuits;
  * exact public fast-core preflight with either a direct complete-square
    fallback or the lower-work 16-orbit fallback.

The direct complete-square route is intentionally highlighted: it supplies a
single core-complete theorem for all nine public parameter sets and eliminates
the preferred-core and orbit-enumeration premises from the simplest headline.

A separate secondary ledger counts only variable-variable scalar products in
extension multiplication.  All fixed linear transforms remain in kappa_hom;
that secondary ledger is explicitly not the primary headline.
"""
from __future__ import annotations
import importlib.util, json, math, sys
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE=Path(__file__).resolve().parent
BUNDLE=HERE

def load(name:str,path:Path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);sys.modules[name]=m;s.loader.exec_module(m);return m
v5=load('snova_v5_for_v6_defensible',BUNDLE/'snova_master_upgrade_v5.py')
base=v5.base
Q=19;RHO=Fraction(14,256);A2=Q**2;A4=Q**4
G1=150;G2=692;G4=2628
REF={'I':143,'III':207,'V':272}
EXCESS_BITS=128

def lg2(x:int|Fraction|float)->float:
    if isinstance(x,Fraction):return math.log2(x.numerator)-math.log2(x.denominator)
    return math.log2(x)

def eta_threshold(h:int,K:int)->Fraction:
    return Fraction(Q**h-1,Q**K)+Fraction(1,2**EXCESS_BITS)

def accepted_ns(alpha:Fraction,h:int,K:int,eta:Fraction)->Fraction:
    a=alpha-eta/(Q-1)
    if a<=0: return Fraction(0,1)
    # If Y counts accepted nonsingular roots, then E Y >= mu*a and
    # E[Y(Y-1)] <= mu*eta.  Hence E[Y^2] <= E[Y]+mu*eta and the
    # Cauchy--Schwarz quotient is minimized at E Y=mu*a.
    return Fraction(1,Q**(K-h))*a*a/(a+eta)

def full_domain_tail(*,d:int,v:int,h:int,K:int,zero_offset:bool)->Fraction:
    t=d*(v-1 if zero_offset else v-2)
    return Fraction(Q**h-1)*(1-Fraction(1,Q**K))*RHO**t*2**EXCESS_BITS

def exact_homotopy_factor(*,base_size:int,B:int,degree_bound:int,gate_cost:int)->dict[str,Any]:
    bad_pairs=B*B
    r0=1
    while base_size**r0<degree_bound or base_size**r0<=bad_pairs:r0+=1
    best=None
    for r in range(r0,r0+80):
        E=base_size**r
        if 2*r-1>base_size:break
        succ=1-Fraction(bad_pairs,E)
        if succ<=0:continue
        fac=Fraction((2*r-1)*gate_cost,1)/succ
        if best is None or fac<best[0]:best=(fac,r,E,succ)
    if best is None:raise RuntimeError('no extension')
    fac,r,E,succ=best
    return {'extension_degree':r,'field_size':E,'separator_success':[succ.numerator,succ.denominator],
            'separator_success_float':float(succ),'pointwise_products':2*r-1,
            'separator_bad_hyperplanes':bad_pairs,
            'gate_factor':[fac.numerator,fac.denominator],'gate_factor_log2':lg2(fac),
            'base_multiplication_gate_cost':gate_cost}

def ff(info:dict[str,Any])->Fraction:return Fraction(*info['gate_factor'])

def target_filter_log2(c:int)->float:return c*math.log2(Q)+32.0

def combine_logs(a:float,b:float)->float:
    hi=max(a,b);lo=min(a,b);return hi+math.log2(1+2**(lo-hi))

def l2_diag(level:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,m,d,r=params;n=v+m;M=4*m;K=3*m;s=m;h=2*s
    alpha=base.exact_l2_acceptance(n);eta=eta_threshold(h,K);delta=Fraction(s,A2)
    mu=Fraction(1,Q**(K-h));a=alpha-delta
    root=mu*max(a-eta,a*a/(a+eta))
    B=2**s;Bp=B+s*2**(s-1);slp=base.corrected_dense_quadratic_slp(s,s);H=slp+2*s+s*s
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*s,gate_cost=G2)
    W=ff(hom)*B*Bp*H*s/root
    out_bytes=packed_parametrization_bytes(dimension=s,B=B,
                                           absolute_extension_degree=2*hom['extension_degree'])
    pjac=base.p_atom_product(2*c for c in range(1,s+1));eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=True)
    if pjac<=eps:raise RuntimeError('bad key density')
    solve=W/(pjac-eps)
    target=target_filter_log2(M-K);combined=combine_logs(lg2(solve),target)
    return {'level':level,'parameters':list(params),'route':'zero-offset one-eigenblock core',
            's':s,'h':h,'K':K,'eta':[eta.numerator,eta.denominator],
            'root_probability_log2':lg2(root),'jacobian_success_lower':float(pjac),
            'jacobian_failure_upper':float(1-pjac),'spectrum_failure_log2':lg2(eps),
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),'solve_coverage_normalized_log2_AXN':lg2(solve),
            'target_filter_log2_AXN':target,'total_log2_AXN':combined,'headroom':ref-combined,
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp),
            'packed_parametrization_bytes_log2':lg2(out_bytes),
            'packed_parametrization_bytes':[out_bytes.numerator,out_bytes.denominator]}

def l2_complete(level:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,m,d,r=params;n=v+m;M=4*m;K=3*m;s=K//2;h=K
    alpha=base.exact_l2_acceptance(n);eta=eta_threshold(h,K);root=accepted_ns(alpha,h,K,eta)
    B=2**(2*m)*math.comb(m,m//2);Bpf=Fraction(B)*(1+2*m+Fraction(m*m,m+2));assert Bpf.denominator==1;Bp=Bpf.numerator
    D=math.comb(s+1,2);slp=2*D+s*s+4*m*(D+s+1)+2*m*(s+1)**2;H=slp+6*m+K*K
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*s,gate_cost=G2)
    W=ff(hom)*B*Bp*H*K/root
    out_bytes=packed_parametrization_bytes(dimension=K,B=B,
                                           absolute_extension_degree=2*hom['extension_degree'])
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=True);solve=W/(1-eps)
    target=target_filter_log2(M-K);combined=combine_logs(lg2(solve),target)
    return {'level':level,'parameters':list(params),'route':'zero-offset complete two-block core',
            's':s,'h':h,'K':K,'eta':[eta.numerator,eta.denominator],
            'root_probability_log2':lg2(root),'spectrum_failure_log2':lg2(eps),'homotopy':hom,
            'per_good_key_log2_AXN':lg2(W),'solve_coverage_normalized_log2_AXN':lg2(solve),
            'target_filter_log2_AXN':target,'total_log2_AXN':combined,'headroom':ref-combined,
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp),
            'packed_parametrization_bytes_log2':lg2(out_bytes),
            'packed_parametrization_bytes':[out_bytes.numerator,out_bytes.denominator]}

def l4_profile(name:str,params:tuple[int,int,int,int],ref:int,s:int,a:int,b:int)->dict[str,Any]:
    v,o,d,r=params;M=o*r*d;m=(o*r+d-1)//d;K=10*m;h=4*s;n=v+o
    alpha=base.l4_acceptance(n,r);eta=eta_threshold(h,K);delta=Fraction(a+Q*b,A4)
    good_density=alpha-delta
    root=Fraction(1,Q**(K-h))*good_density*good_density/(good_density+eta)
    B=2**a*(Q+1)**b;Bpf=Fraction(B)*(1+Fraction(a,2)+Fraction(b,Q+1));assert Bpf.denominator==1;Bp=Bpf.numerator
    D=math.comb(s+1,2);slp0=0 if a==0 else D+2*a*(D+s+1);slp1=0 if b==0 else 6*s+s*s+2*b*(s*s+2*s+1)
    deg=2*a+(Q+1)*b;H=slp0+slp1+deg+s*s
    hom=exact_homotopy_factor(base_size=A4,B=B,degree_bound=deg,gate_cost=G4)
    W=ff(hom)*B*Bp*H*s/root
    ppre=base.structural_preflight(d,v,o,m,M)
    # Conservative projective right-kernel union bound.  For every fixed
    # nonzero kernel direction, vanishing of all s selected rows imposes 4s
    # independent F19 constraints; union over A-projective directions.
    lines=Fraction(A4**s-1,A4-1)
    jacfail=lines*RHO**(4*s)
    pjac=1-jacfail
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=False);pkey=ppre*(pjac-eps)
    norm=W/pkey
    return {'name':name,'level':name.split('-')[0],'parameters':list(params),'profile':[s,a,b],
            'h':h,'K':K,'eta':[eta.numerator,eta.denominator],'root_probability_log2':lg2(root),
            'jacobian_success_lower':float(pjac),'jacobian_failure_upper':float(jacfail),'jacobian_failure_log2':lg2(jacfail),
            'spectrum_failure_log2':lg2(eps),'structural_probability':float(ppre),
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),'fixed_core_normalized_log2_AXN':lg2(norm),
            'headroom':ref-lg2(norm),'B_log2':lg2(B),'Bplus_log2':lg2(Bp)}

def l4_frontier(name:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,o,d,r=params;m=(o*r+d-1)//d;Hdim=v+o-d*m;rows=[]
    for a in range(m+1):
        for b in range(m+1):
            s=a+b
            if s and s<=Hdim and 4*s<=10*m:rows.append(l4_profile(name,params,ref,s,a,b))
    return {'fast':min(rows,key=lambda z:z['fixed_core_normalized_log2_AXN']),'profiles_tested':len(rows)}

def packed_parametrization_bytes(*,dimension:int,B:int,absolute_extension_degree:int)->Fraction:
    """Packed (dimension+1)B-field-element output ceiling."""
    return Fraction((dimension+1)*B*absolute_extension_degree*Q.bit_length(),8)

def l4_complete_square(name:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    """Direct K-by-K ordinary quadratic fallback on the guaranteed fresh slice.

    The coefficients lie in F19.  We embed them in F19^2 and choose the
    homotopy extension over that base because this gives the best complete
    constructive multiplication ledger among the validated towers.
    """
    v,o,d,r=params;M=o*r*d;m=(o*r+d-1)//d;K=10*m;h=K;n=v+o
    alpha=base.l4_acceptance(n,r);eta=eta_threshold(h,K)
    root=accepted_ns(alpha,h,K,eta)
    B=2**K;Bpf=Fraction(B)*(1+Fraction(K,2));assert Bpf.denominator==1;Bp=Bpf.numerator
    slp=base.corrected_dense_quadratic_slp(K,K);H=slp+2*K+K*K
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*K,gate_cost=G2)
    W=ff(hom)*B*Bp*H*K/root
    ppre=base.structural_preflight(d,v,o,m,M)
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=False)
    norm=W/(ppre*(1-eps))
    absolute_degree=2*hom['extension_degree']
    out_bytes=packed_parametrization_bytes(dimension=K,B=B,absolute_extension_degree=absolute_degree)
    return {'name':name,'level':name.split('-')[0],'parameters':list(params),
            'route':'direct complete ordinary square','h':h,'K':K,
            'eta':[eta.numerator,eta.denominator],'root_probability_log2':lg2(root),
            'spectrum_failure_log2':lg2(eps),'structural_probability':float(ppre),
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),
            'normalized_log2_AXN':lg2(norm),'headroom':ref-lg2(norm),
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp),
            'packed_parametrization_bytes_log2':lg2(out_bytes),
            'packed_parametrization_bytes':[out_bytes.numerator,out_bytes.denominator]}

def p_root_orbit(level:str,eta:Fraction)->Fraction:
    h={'I':48,'III':68,'V':88}[level];K={'I':50,'III':70,'V':90}[level]
    alpha=v5.v4._orbit_acceptance(level)
    return accepted_ns(alpha,h,K,eta)

def p_root_orbit_v5(level:str,eta:Fraction)->Fraction:
    h={'I':48,'III':68,'V':88}[level];K={'I':50,'III':70,'V':90}[level]
    alpha=v5.v4._orbit_acceptance(level);a=alpha-eta/(Q-1)
    return Fraction(1,Q**(K-h))*a*a/(1+eta)

def orbit_fallback(level:str,shape:str)->dict[str,Any]:
    # Re-cost the pinned v5 16-orbit schedule.  The v5 schedule retains the
    # conservative 8/7 separator reserve and is therefore safe to scale by the
    # complete F19^4 multiplication circuit ratio.
    ref=REF[level];h={'I':48,'III':68,'V':88}[level];K={'I':50,'III':70,'V':90}[level]
    params={('I','a'):(28,5,4,4),('I','b'):(28,4,4,5),('III','a'):(40,7,4,4),('III','b'):(38,5,4,5),('V','a'):(50,9,4,4),('V','b'):(52,6,4,6)}[(level,shape)]
    old_eta=Fraction(Q**h-1,Q**K)+Fraction(1,2**32);new_eta=eta_threshold(h,K)
    root_shift=lg2(p_root_orbit_v5(level,old_eta)/p_root_orbit(level,new_eta))
    field_shift=math.log2(G4/3164)
    if level=='III' and shape=='b':
        old=v5.zero_offset_level3b_orbit_v5(); oldnorm=old['solve_normalized_log2_AXN']
        # v5's zero-offset spectrum used the correct full-vinegar theorem but a 2^-32 threshold.
        eps=full_domain_tail(d=4,v=params[0],h=h,K=K,zero_offset=True)
        solve=oldnorm+field_shift+root_shift
        target=old['target_generation_envelope_log2_AXN'];total=combine_logs(solve,target)
        ppre=base.structural_preflight(4,params[0],params[1],(params[1]*params[3]+3)//4,params[1]*params[3]*4)
        raw=2**solve*float(ppre) # recover the per-key solve ledger used for adaptive comparison
        return {'level':level,'shape':shape,'parameters':list(params),'route':'zero-offset 16-orbit fallback',
                'normalized_log2_AXN':total,'solve_normalized_log2_AXN':solve,'per_good_key_log2_AXN':math.log2(raw),
                'target_filter_log2_AXN':target,'spectrum_failure_log2':lg2(eps),'structural_probability':float(ppre),
                'headroom':ref-total,'field_shift_bits':field_shift,'root_shift_bits':root_shift}
    old=v5.orbit_v5(level)['v5_projection_complete_log2_AXN']
    # v5 was normalized by the exact structural preflight; recover raw work.
    v,o,d,r=params;M=o*r*d;m=(o*r+d-1)//d;ppre=base.structural_preflight(d,v,o,m,M)
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=False)
    norm=old+field_shift+root_shift
    raw=2**norm*float(ppre)
    return {'level':level,'shape':shape,'parameters':list(params),'route':'16-orbit core-complete fallback',
            'normalized_log2_AXN':norm,'per_good_key_log2_AXN':math.log2(raw),
            'spectrum_failure_log2':lg2(eps),'structural_probability':float(ppre),'headroom':ref-norm,
            'field_shift_bits':field_shift,'root_shift_bits':root_shift}

def adaptive_l4(fast:dict[str,Any],fallback:dict[str,Any],*,label:str)->dict[str,Any]:
    ppre=Fraction.from_float(fast['structural_probability']);epsf=2**fast['spectrum_failure_log2'];epsb=2**fallback['spectrum_failure_log2'];eps=min(0.5,epsf+epsb)
    Wf=2**fast['per_good_key_log2_AXN'];WF=2**fallback['per_good_key_log2_AXN'];delta=fast['jacobian_failure_upper']
    qfail=min(1.0,delta/(1-eps));W=Wf+qfail*max(0.0,WF-Wf)
    norm=W/(float(ppre)*(1-eps))
    return {'route':label,'normalized_log2_AXN':math.log2(norm),
            'headroom':REF[fast['level']]-math.log2(norm),'fast_preflight_failure_upper':delta,
            'common_spectrum_failure_upper':eps,'fast_per_good_key_log2':fast['per_good_key_log2_AXN'],
            'fallback_per_good_key_log2':fallback['per_good_key_log2_AXN']}

def adaptive_l2(diag:dict[str,Any],comp:dict[str,Any])->dict[str,Any]:
    # Same zero-offset target filter for both routes; solve dominates, so mix solve
    # work first and add the common target-generation ledger once.
    eps=min(0.5,2**diag['spectrum_failure_log2']+2**comp['spectrum_failure_log2'])
    Wf=2**diag['per_good_key_log2_AXN'];WF=2**comp['per_good_key_log2_AXN'];delta=diag['jacobian_failure_upper'];qfail=min(1.0,delta/(1-eps))
    W=Wf+qfail*max(0.0,WF-Wf);solve=math.log2(W/(1-eps));target=diag['target_filter_log2_AXN'];total=combine_logs(solve,target)
    return {'route':'certified one-eigenblock core, then complete two-block fallback','normalized_log2_AXN':total,
            'headroom':REF[diag['level']]-total,'fast_preflight_failure_upper':delta,'common_spectrum_failure_upper':eps,
            'fast_per_good_key_log2':diag['per_good_key_log2_AXN'],'fallback_per_good_key_log2':comp['per_good_key_log2_AXN']}

def output_oriented_l2(diag:dict[str,Any],comp:dict[str,Any])->dict[str,Any]:
    """Use the smaller-output diagonal core whenever its certificate passes."""
    eps=min(0.5,2**diag['spectrum_failure_log2']+2**comp['spectrum_failure_log2'])
    Wf=2**diag['per_good_key_log2_AXN'];WF=2**comp['per_good_key_log2_AXN']
    delta=diag['jacobian_failure_upper'];qfail=min(1.0,delta/(1-eps))
    W=(1-qfail)*Wf+qfail*WF
    solve=math.log2(W/(1-eps));target=diag['target_filter_log2_AXN'];total=combine_logs(solve,target)
    return {'route':'smaller-output one-eigenblock core, complete fallback on preflight failure',
            'normalized_log2_AXN':total,'headroom':REF[diag['level']]-total,
            'fast_preflight_failure_upper':delta,'fast_branch_probability_lower':1-qfail,
            'common_spectrum_failure_upper':eps,'fast_output_bytes_log2':diag['packed_parametrization_bytes_log2'],
            'fallback_output_bytes_log2':comp['packed_parametrization_bytes_log2'],
            'fast_per_good_key_log2':diag['per_good_key_log2_AXN'],
            'fallback_per_good_key_log2':comp['per_good_key_log2_AXN']}

# Secondary multiplicative-core factor: only variable-variable F19 products are
# charged; every fixed linear transform remains in kappa_hom.
def scalar_mult_schedule(N:int)->tuple[int,int,int]|None:
    best=None
    for k in range(1,N+1):
        if N%k:continue
        d=N//k
        if 2*d-1>Q**k or 2*k-1>Q:continue
        c=(2*d-1)*(2*k-1)
        if best is None or c<best[0]:best=(c,k,d)
    return best

def multiplicative_factor(*,absolute_base_degree:int,B:int,degree_bound:int)->dict[str,Any]:
    best=None
    for r in range(1,120):
        N=absolute_base_degree*r;E=Q**N
        if E<degree_bound or E<=B*B:continue
        sc=scalar_mult_schedule(N)
        if sc is None:continue
        c,k,d=sc;succ=1-Fraction(B*B,E);fac=Fraction(c*G1,1)/succ
        if best is None or fac<best[0]:best=(fac,r,N,E,succ,c,k,d)
    if best is None:raise RuntimeError('no bilinear schedule')
    fac,r,N,E,succ,c,k,d=best
    return {'extension_degree':r,'absolute_degree':N,'field_size':E,'scalar_products':c,
            'intermediate_degree':k,'outer_degree':d,'separator_success_float':float(succ),
            'gate_factor_log2':lg2(fac),'gate_factor':[fac.numerator,fac.denominator]}

def secondary_recost(route:dict[str,Any],base_degree:int)->float:
    # Replace only the primary homotopy factor; all other route factors stay.
    B=round(2**route['B_log2']); # exact for powers/combinatorial B may exceed float precision; reconstruct below only used l2 chosen rows.
    raise NotImplementedError

def _level_summary(report:dict[str,Any],selector)->dict[str,Any]:
    out={}
    for L in ['I','III','V']:
        vals=[selector(r) for r in report['rows'].values() if r['level']==L]
        x=max(vals)
        out[L]={'exponent':x,'headroom':REF[L]-x,'plus7_headroom':REF[L]-x-7}
    return out

def main()->None:
    report={
      'constants':{
        'q':Q,'rho':[RHO.numerator,RHO.denominator],
        'spectrum_excess_bits':EXCESS_BITS,
        'separator_bad_set':'at most B^2 forbidden vectors',
        'accepted_root_denominator':'a+eta',
        'multiplication_AXN':{'F19':G1,'F19^2':G2,'F19^4':G4}},
      'rows':{},'orbit_fallbacks':{}
    }
    # ell=2: the complete two-block square gives the simple theorem.  The
    # Level-I optimized route first tries the one-eigenblock core.  A separate
    # output-oriented ledger uses that smaller core at every level.
    for level,params,ref in base.L2_ROWS:
        dg=l2_diag(level,params,ref);cp=l2_complete(level,params,ref)
        opt=adaptive_l2(dg,cp) if level=='I' else {
            'route':cp['route'],'normalized_log2_AXN':cp['total_log2_AXN'],
            'headroom':cp['headroom']}
        compact=output_oriented_l2(dg,cp)
        report['rows'][str(params)]={
            'level':level,'parameters':list(params),'diagonal':dg,'complete':cp,
            'optimized':opt,'compact_output':compact,
            'simple_complete_square_exponent':cp['total_log2_AXN']}

    # ell=4: the direct K-by-K ordinary square is the simple core-complete
    # fallback.  The legacy 16-orbit schedule is retained only as a lower-work
    # optional refinement of the fallback.
    for name,params,ref in base.L4_ROWS:
        fr=l4_frontier(name,params,ref);fast=fr['fast']
        level=name.split('-')[0];shape=name.split('-')[1]
        direct=l4_complete_square(name,params,ref)
        orbit=orbit_fallback(level,shape)
        simple_ad=adaptive_l4(
            fast,direct,label='certified fast channel core, then direct complete ordinary square')
        optimized_ad=adaptive_l4(
            fast,orbit,label='certified fast channel core, then lower-work 16-orbit fallback')
        report['rows'][str(params)]={
            'name':name,'level':level,'parameters':list(params),'fast':fast,
            'complete_square':direct,'orbit_fallback':orbit,
            'simple_adaptive':simple_ad,'optimized':optimized_ad,
            'profiles_tested':fr['profiles_tested'],
            'simple_complete_square_exponent':direct['normalized_log2_AXN']}
        report['orbit_fallbacks'][name]=orbit

    def simple_selector(r): return r['simple_complete_square_exponent']
    def simple_adaptive_selector(r):
        return r['optimized']['normalized_log2_AXN'] if 'complete' in r else r['simple_adaptive']['normalized_log2_AXN']
    def optimized_selector(r): return r['optimized']['normalized_log2_AXN']
    def compact_selector(r):
        return r['compact_output']['normalized_log2_AXN'] if 'complete' in r else r['simple_adaptive']['normalized_log2_AXN']

    report['all_nine_simple_complete_square']=_level_summary(report,simple_selector)
    report['all_nine_adaptive_direct_fallback']=_level_summary(report,simple_adaptive_selector)
    report['all_nine_optimized']=_level_summary(report,optimized_selector)
    report['all_nine_compact_output']=_level_summary(report,compact_selector)

    out=HERE/'snova_master_upgrade_v7_primary.json'
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('v7 primary ledger (audited conservative bounds)')
    for title in ['all_nine_simple_complete_square','all_nine_adaptive_direct_fallback',
                  'all_nine_optimized','all_nine_compact_output']:
        print(title)
        for L,z in report[title].items(): print(L,z)
    print('rows')
    for p,r in report['rows'].items():
        print(p,'simple',r['simple_complete_square_exponent'],
              'optimized',r['optimized']['normalized_log2_AXN'])
    print('wrote',out)
if __name__=='__main__':main()
