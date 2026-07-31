#!/usr/bin/env python3
"""Exact ledger for the SNOVA symmetry-quotient manuscript.

This module is intentionally self-contained.  It implements only the routes
used in the paper: the complete square on every row and the certified
fast-core/direct-fallback optimization.  It assigns no random probability to
SNOVA's fixed outer maps.  The resulting costs are conditional on the exact
public structural preflight stated in the paper.
"""
from __future__ import annotations
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

HERE=Path(__file__).resolve().parent
Q=19;RHO=Fraction(14,256);A2=Q**2;A4=Q**4
G1=150;G2=692;G4=2628
REF={'I':143,'III':207,'V':272}
EXCESS_BITS=128
L2_ROWS=(
    ('I',(48,16,2,2),143),
    ('III',(72,24,2,2),207),
    ('V',(96,32,2,2),272),
)
L4_ROWS=(
    ('I-a',(28,5,4,4),143),
    ('I-b',(28,4,4,5),143),
    ('III-a',(40,7,4,4),207),
    ('III-b',(38,5,4,5),207),
    ('V-a',(50,9,4,4),272),
    ('V-b',(52,6,4,6),272),
)

def product(values:Iterable[Fraction])->Fraction:
    out=Fraction(1,1)
    for value in values:out*=value
    return out

def p_atom_product(exponents:Iterable[int])->Fraction:
    return product(1-RHO**e for e in exponents)

def exact_l2_acceptance(n:int)->Fraction:
    return Fraction(sum(math.comb(n,j)*(Q-1)**(n-j)
                        for j in range(n//4+1)),Q**n)

def l4_acceptance_lower_bound(n:int,r:int)->Fraction:
    return (1-Fraction(1,Q**4))**n if r==4 else Fraction(1,1)

def dense_quadratic_slp(nvars:int,neqs:int)->int:
    quadratic_monomials=math.comb(nvars+1,2)
    terms_per_equation=quadratic_monomials+nvars+1
    return quadratic_monomials+2*neqs*terms_per_equation

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
    e0=1
    while base_size**e0<degree_bound or base_size**e0<=bad_pairs:e0+=1
    best=None
    for e in range(e0,e0+80):
        E=base_size**e
        if 2*e-1>base_size:break
        succ=1-Fraction(bad_pairs,E)
        if succ<=0:continue
        fac=Fraction((2*e-1)*gate_cost,1)/succ
        if best is None or fac<best[0]:best=(fac,e,E,succ)
    if best is None:raise RuntimeError('no extension')
    fac,e,E,succ=best
    return {'extension_degree':e,'field_size':E,'separator_success':[succ.numerator,succ.denominator],
            'separator_success_float':float(succ),'pointwise_products':2*e-1,
            'separator_bad_hyperplanes':bad_pairs,
            'gate_factor':[fac.numerator,fac.denominator],'gate_factor_log2':lg2(fac),
            'base_multiplication_gate_cost':gate_cost}

def ff(info:dict[str,Any])->Fraction:return Fraction(*info['gate_factor'])

def target_filter_log2(c:int)->float:return c*math.log2(Q)+32.0

def combine_logs(a:float,b:float)->float:
    hi=max(a,b);lo=min(a,b);return hi+math.log2(1+2**(lo-hi))

def l2_diag(level:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,m,d,r=params;n=v+m;M=4*m;K=3*m;s=m;h=2*s
    alpha=exact_l2_acceptance(n);eta=eta_threshold(h,K);delta=Fraction(s,A2)
    mu=Fraction(1,Q**(K-h));a=alpha-delta
    root=mu*max(a-eta,a*a/(a+eta))
    B=2**s;Bp=B+s*2**(s-1);slp=dense_quadratic_slp(s,s);H=slp+2*s+s*s
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*s,gate_cost=G2)
    W=ff(hom)*B*Bp*H*s/root
    pjac=p_atom_product(2*c for c in range(1,s+1));eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=True)
    if pjac<=eps:raise RuntimeError('bad key density')
    solve=W/(pjac-eps)
    target=target_filter_log2(M-K);combined=combine_logs(lg2(solve),target)
    return {'level':level,'parameters':list(params),'route':'zero-offset one-eigenblock core',
            's':s,'h':h,'K':K,'eta':[eta.numerator,eta.denominator],
            'root_probability_log2':lg2(root),'jacobian_success_lower':float(pjac),
            'jacobian_failure_upper':float(1-pjac),'spectrum_failure_log2':lg2(eps),
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),'solve_coverage_normalized_log2_AXN':lg2(solve),
            'target_filter_log2_AXN':target,'total_log2_AXN':combined,'headroom':ref-combined,
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp)}

def l2_complete(level:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,m,d,r=params;n=v+m;M=4*m;K=3*m;s=K//2;h=K
    alpha=exact_l2_acceptance(n);eta=eta_threshold(h,K);root=accepted_ns(alpha,h,K,eta)
    B=2**(2*m)*math.comb(m,m//2);Bpf=Fraction(B)*(1+2*m+Fraction(m*m,m+2));assert Bpf.denominator==1;Bp=Bpf.numerator
    D=math.comb(s+1,2);slp=2*D+s*s+4*m*(D+s+1)+2*m*(s+1)**2;H=slp+6*m+K*K
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*s,gate_cost=G2)
    W=ff(hom)*B*Bp*H*K/root
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=True);solve=W/(1-eps)
    target=target_filter_log2(M-K);combined=combine_logs(lg2(solve),target)
    return {'level':level,'parameters':list(params),'route':'zero-offset complete two-block core',
            's':s,'h':h,'K':K,'eta':[eta.numerator,eta.denominator],
            'root_probability_log2':lg2(root),'spectrum_failure_log2':lg2(eps),'homotopy':hom,
            'per_good_key_log2_AXN':lg2(W),'solve_coverage_normalized_log2_AXN':lg2(solve),
            'target_filter_log2_AXN':target,'total_log2_AXN':combined,'headroom':ref-combined,
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp)}

def l4_profile(name:str,params:tuple[int,int,int,int],ref:int,s:int,a:int,b:int)->dict[str,Any]:
    v,o,d,r=params;m=(o*r+d-1)//d;K=10*m;h=4*s;n=v+o
    alpha=l4_acceptance_lower_bound(n,r);eta=eta_threshold(h,K);delta=Fraction(a+Q*b,A4)
    good_density=alpha-delta
    root=Fraction(1,Q**(K-h))*good_density*good_density/(good_density+eta)
    B=2**a*(Q+1)**b;Bpf=Fraction(B)*(1+Fraction(a,2)+Fraction(b,Q+1));assert Bpf.denominator==1;Bp=Bpf.numerator
    D=math.comb(s+1,2);slp0=0 if a==0 else D+2*a*(D+s+1);slp1=0 if b==0 else 6*s+s*s+2*b*(s*s+2*s+1)
    deg=2*a+(Q+1)*b;H=slp0+slp1+deg+s*s
    hom=exact_homotopy_factor(base_size=A4,B=B,degree_bound=deg,gate_cost=G4)
    W=ff(hom)*B*Bp*H*s/root
    # Conservative projective right-kernel union bound.  For every fixed
    # nonzero kernel direction, vanishing of all s selected rows imposes 4s
    # independent F19 constraints; union over A-projective directions.
    lines=Fraction(A4**s-1,A4-1)
    jacfail=lines*RHO**(4*s)
    pjac=1-jacfail
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=False);pkey=pjac-eps
    norm=W/pkey
    return {'name':name,'level':name.split('-')[0],'parameters':list(params),'profile':[s,a,b],
            'h':h,'K':K,'eta':[eta.numerator,eta.denominator],'root_probability_log2':lg2(root),
            'jacobian_success_lower':float(pjac),'jacobian_failure_upper':float(jacfail),'jacobian_failure_log2':lg2(jacfail),
            'spectrum_failure_log2':lg2(eps),'structural_preflight':'required exact public preflight',
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),'fixed_core_normalized_log2_AXN':lg2(norm),
            'headroom':ref-lg2(norm),'B_log2':lg2(B),'Bplus_log2':lg2(Bp)}

def l4_frontier(name:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    v,o,d,r=params;m=(o*r+d-1)//d;Hdim=v+o-d*m;rows=[]
    for a in range(m+1):
        for b in range(m+1):
            s=a+b
            if s and s<=Hdim and 4*s<=10*m:rows.append(l4_profile(name,params,ref,s,a,b))
    return {'fast':min(rows,key=lambda z:z['fixed_core_normalized_log2_AXN']),'profiles_tested':len(rows)}

def l4_complete_square(name:str,params:tuple[int,int,int,int],ref:int)->dict[str,Any]:
    """Direct K-by-K ordinary quadratic fallback on the guaranteed fresh slice.

    The coefficients lie in F19.  We embed them in F19^2 and choose the
    homotopy extension over that base because this gives the best complete
    constructive multiplication ledger among the validated towers.
    """
    v,o,d,r=params;m=(o*r+d-1)//d;K=10*m;h=K;n=v+o
    alpha=l4_acceptance_lower_bound(n,r);eta=eta_threshold(h,K)
    root=accepted_ns(alpha,h,K,eta)
    B=2**K;Bpf=Fraction(B)*(1+Fraction(K,2));assert Bpf.denominator==1;Bp=Bpf.numerator
    slp=dense_quadratic_slp(K,K);H=slp+2*K+K*K
    hom=exact_homotopy_factor(base_size=A2,B=B,degree_bound=2*K,gate_cost=G2)
    W=ff(hom)*B*Bp*H*K/root
    eps=full_domain_tail(d=d,v=v,h=h,K=K,zero_offset=False)
    norm=W/(1-eps)
    return {'name':name,'level':name.split('-')[0],'parameters':list(params),
            'route':'direct complete ordinary square','h':h,'K':K,
            'eta':[eta.numerator,eta.denominator],'root_probability_log2':lg2(root),
            'spectrum_failure_log2':lg2(eps),
            'structural_preflight':'required exact public preflight',
            'homotopy':hom,'per_good_key_log2_AXN':lg2(W),
            'normalized_log2_AXN':lg2(norm),'headroom':ref-lg2(norm),
            'B_log2':lg2(B),'Bplus_log2':lg2(Bp)}

def adaptive_l4(fast:dict[str,Any],fallback:dict[str,Any],*,label:str)->dict[str,Any]:
    epsf=2**fast['spectrum_failure_log2'];epsb=2**fallback['spectrum_failure_log2'];eps=min(0.5,epsf+epsb)
    Wf=2**fast['per_good_key_log2_AXN'];WF=2**fallback['per_good_key_log2_AXN'];delta=fast['jacobian_failure_upper']
    qfail=min(1.0,delta/(1-eps));W=Wf+qfail*max(0.0,WF-Wf)
    norm=W/(1-eps)
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

def _level_summary(report:dict[str,Any],selector)->dict[str,Any]:
    out={}
    for L in ['I','III','V']:
        vals=[selector(r) for r in report['rows'].values() if r['level']==L]
        x=max(vals)
        out[L]={'exponent':x,'headroom':REF[L]-x}
    return out

def main()->None:
    report={
      'constants':{
        'q':Q,'rho_byte':[RHO.numerator,RHO.denominator],
        'spectrum_excess_bits':EXCESS_BITS,
        'separator_bad_set':'at most B^2 forbidden vectors',
        'accepted_root_denominator':'a+eta',
        'multiplication_AXN':{'F19':G1,'F19^2':G2,'F19^4':G4}},
      'scope':(
          'unit-leading AXN multiplication-schedule ledger, conditional on '
          'the exact public structural preflight; kappa_hom excluded'
      ),
      'rows':{}
    }
    # ell=2: the complete two-block square gives the direct theorem.  The
    # Level-I adaptive route first tries the one-eigenblock core.
    for level,params,ref in L2_ROWS:
        dg=l2_diag(level,params,ref);cp=l2_complete(level,params,ref)
        opt=adaptive_l2(dg,cp) if level=='I' else {
            'route':cp['route'],'normalized_log2_AXN':cp['total_log2_AXN'],
            'headroom':cp['headroom']}
        report['rows'][str(params)]={
            'level':level,'parameters':list(params),'diagonal':dg,'complete':cp,
            'adaptive':opt,
            'direct_complete_square_exponent':cp['total_log2_AXN']}

    # ell=4: the direct K-by-K ordinary square is the complete
    # fallback.  The adaptive branch uses only a deterministic fast-core
    # certificate plus that same complete fallback.
    for name,params,ref in L4_ROWS:
        fr=l4_frontier(name,params,ref);fast=fr['fast']
        level=name.split('-')[0]
        direct=l4_complete_square(name,params,ref)
        adaptive=adaptive_l4(
            fast,direct,label='certified fast channel core, then direct complete ordinary square')
        report['rows'][str(params)]={
            'name':name,'level':level,'parameters':list(params),'fast':fast,
            'complete_square':direct,'adaptive':adaptive,
            'profiles_tested':fr['profiles_tested'],
            'direct_complete_square_exponent':direct['normalized_log2_AXN']}

    def direct_selector(r): return r['direct_complete_square_exponent']
    def adaptive_selector(r): return r['adaptive']['normalized_log2_AXN']

    report['all_nine_direct_complete_square']=_level_summary(report,direct_selector)
    report['all_nine_adaptive']=_level_summary(report,adaptive_selector)

    out=HERE/'all_nine_ledger.json'
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print('SNOVA symmetry-quotient ledger (conditional bounds)')
    for title in ['all_nine_direct_complete_square','all_nine_adaptive']:
        print(title)
        for L,z in report[title].items(): print(L,z)
    print('rows')
    for p,r in report['rows'].items():
        print(p,'direct',r['direct_complete_square_exponent'],
              'adaptive',r['adaptive']['normalized_log2_AXN'])
    print('wrote',out)
if __name__=='__main__':main()
