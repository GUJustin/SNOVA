#!/usr/bin/env python3
"""Verify the official Level-I 50-in-52 cross-column certificate.

Requires the Version 2.3 Level-I KAT response file.  The script reconstructs
the public key, computes the official ABQ feature map, verifies the six
80-by-80 cross-column ranks and four 80-by-50 self-column ranks, and checks
the explicit rank-80 x0 certificate against the full verifier.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

from symmetry_attack_validation import (
    PARAMS, Q, affine_solve, direct_output, kat_public_key, kronS,
    nullspace_mod, official_S, parse_kat, powers, rank_mod, reconstruct_abq,
)

EXPLICIT_X0 = np.array([
18,12,6,8,5,7,6,3,1,7,14,16,4,14,7,15,17,8,1,15,11,5,1,11,13,15,16,13,15,11,12,7,
13,18,3,11,12,0,13,6,14,14,0,12,2,0,10,8,18,10,3,8,8,11,14,9,10,18,0,9,14,11,18,14,
3,10,5,15,0,17,5,3,8,9,15,7,1,3,13,11,1,0,8,16,7,8,14,14,7,7,9,3,7,7,18,8,13,9,
13,10,9,3,2,13,0,1,13,18,0,13,2,17,16,9,5,1,14,6,4,8,11,7,11,4,5,5,14,14,15,18,9,18
], dtype=np.int64)


def basis_mul_tensor(S: np.ndarray) -> np.ndarray:
    l=S.shape[0]
    P=powers(S,2*l)
    B=np.stack([P[u].reshape(-1) for u in range(l)],axis=1)%Q
    T=np.zeros((l,l,l),dtype=np.int64)
    for d in range(l):
        for a in range(l):
            x0, K, _, _ = affine_solve(B,P[d+a].reshape(-1),Q)
            if K.shape[1]:
                raise RuntimeError('power basis is not independent')
            T[d,a]=x0
    return T


def mul_matrix(coeff: np.ndarray,T: np.ndarray) -> np.ndarray:
    return np.einsum('d,dau->ua',np.asarray(coeff,dtype=np.int64)%Q,T,optimize=True)%Q


def full_feature_map(p,abq):
    """Public 80 x 680 map for the identity relation on all four columns."""
    r=p.r;l=p.l;c=r
    C=np.zeros((r,c,l),dtype=np.int64)
    for j in range(r): C[j,j,0]=1
    T=basis_mul_tensor(abq.S)
    labels=[(s,a) for s in range(c) for a in range(l)]
    pairs=[(labels[u],labels[v]) for u in range(len(labels)) for v in range(u,len(labels))]
    pair_index={z:i for i,z in enumerate(pairs)}
    E=np.zeros((p.outputs,p.m1*len(pairs)),dtype=np.int64)
    MR=np.empty((r,c,l,l),dtype=np.int64)
    for j in range(r):
        for s in range(c): MR[j,s]=mul_matrix(C[j,s],T)
    for i in range(p.o):
        row0=i*p.r*p.l
        for al in range(p.alpha):
            pp=(i+al)%p.m1; A=abq.A[i,al]; B=abq.B[i,al]
            L=np.einsum('jsua,a->jsu',MR,abq.q1[i,al],optimize=True)%Q
            R=np.einsum('ktvb,b->ktv',MR,abq.q2[i,al],optimize=True)%Q
            for j in range(p.r):
                for k in range(p.r):
                    outer=np.outer(A[:,j],B[k,:]).reshape(-1)%Q
                    if not np.any(outer): continue
                    for s in range(c):
                        for t in range(c):
                            for u in range(l):
                                lu=int(L[j,s,u])
                                if not lu: continue
                                iu=s*l+u
                                for v in range(l):
                                    rv=int(R[k,t,v])
                                    if not rv: continue
                                    iv=t*l+v; uv=(iu,iv) if iu<=iv else (iv,iu)
                                    col=pp*len(pairs)+pair_index[uv]
                                    sl=slice(row0,row0+p.r*p.l)
                                    E[sl,col]=(E[sl,col]+lu*rv*outer)%Q
    return E,pairs


def select_columns(pairs, p, s, t):
    F=len(pairs); out=[]
    for pp in range(p.m1):
        for k,((s0,a),(s1,b)) in enumerate(pairs):
            if s==t:
                take=(s0==s1==s)
            else:
                take=(s0==s and s1==t)
            if take: out.append(pp*F+k)
    return out


def build_raw_forms(P,l=4):
    n=P.shape[1]; Sps=powers(official_S(l),l)
    Sig=[kronS(Sps[a],n//l) for a in range(l)]
    inv2=pow(2,Q-2,Q); H=[]; M=[]
    for pp in range(P.shape[0]):
        for a in range(l):
            for b in range(a,l):
                A=Sig[a]@P[pp]@Sig[b]%Q
                H.append((A+A.T)*inv2%Q)
        for a in range(l):
            for b in range(l):
                M.append(Sig[a]@P[pp]@Sig[b]%Q)
    return np.asarray(H,dtype=np.int64),np.asarray(M,dtype=np.int64)


def self_features(x,H):
    return np.einsum('i,kij,j->k',x,H,x,optimize=True)%Q


def cross_features(x,y,M):
    return np.einsum('i,kij,j->k',x,M,y,optimize=True)%Q


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--kat',default=str(root/'source_snapshots'/'PQCsignKAT_SNOVA_28_5_19_4.txt'))
    ap.add_argument('--out',default=str(root/'results'/'cross_column_certificate.json'))
    ap.add_argument('--trials',type=int,default=10)
    args=ap.parse_args()
    kat_path=Path(args.kat)
    if not kat_path.exists():
        raise FileNotFoundError(f'KAT not found: {kat_path}. See README.md.')

    p=PARAMS[0]; d=parse_kat(kat_path)
    pk,P=kat_public_key(bytes.fromhex(d['sk']),p)
    kat_match=(pk==bytes.fromhex(d['pk']))
    abq=reconstruct_abq(p); Efull,pairs=full_feature_map(p,abq)

    self_maps={}; cross_maps={}
    for s in range(p.r):
        A=Efull[:,select_columns(pairs,p,s,s)]
        self_maps[str(s)]={'shape':list(A.shape),'rank':rank_mod(A)}
    for s in range(p.r):
        for t in range(s+1,p.r):
            A=Efull[:,select_columns(pairs,p,s,t)]
            cross_maps[f'{s}{t}']={'shape':list(A.shape),'rank':rank_mod(A)}

    Es0=Efull[:,select_columns(pairs,p,0,0)]
    Es1=Efull[:,select_columns(pairs,p,1,1)]
    Ec=Efull[:,select_columns(pairs,p,0,1)]
    H,M=build_raw_forms(P,p.l)
    C_raw=np.einsum('i,kij->kj',EXPLICIT_X0,M,optimize=True)%Q
    C_out=Ec@C_raw%Q
    Nrows=nullspace_mod(C_out,Q); N=Nrows.T
    Hr=np.einsum('ia,kij,jb->kab',N,H,N,optimize=True)%Q
    span_rank=rank_mod(Hr.reshape(Hr.shape[0],-1),Q)

    rng=np.random.default_rng(20260719); equality=True
    for _ in range(args.trials):
        x1=rng.integers(0,Q,size=p.variables,dtype=np.int64)
        U=np.zeros((p.variables,p.r),dtype=np.int64);U[:,0]=EXPLICIT_X0;U[:,1]=x1
        y=direct_output(p,abq,P,U)
        rhs=(Es0@self_features(EXPLICIT_X0,H)+Ec@cross_features(EXPLICIT_X0,x1,M)+Es1@self_features(x1,H))%Q
        equality &= np.array_equal(y,rhs)

    out={
        'parameter':{'v':p.v,'o':p.o,'q':Q,'l':p.l,'r':p.r},
        'kat_public_key_match':bool(kat_match),
        'full_feature_map':{'shape':list(Efull.shape),'rank':rank_mod(Efull)},
        'self_column_maps':self_maps,
        'cross_column_maps':cross_maps,
        'explicit_x0':EXPLICIT_X0.tolist(),
        'cross_coefficient_map':{'shape':list(C_out.shape),'rank':rank_mod(C_out),'kernel_dimension':N.shape[1]},
        'restricted_self_quadratic_span_rank':span_rank,
        'direct_verifier_equalities':{'trials':args.trials,'all_passed':bool(equality)},
    }
    path=Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
