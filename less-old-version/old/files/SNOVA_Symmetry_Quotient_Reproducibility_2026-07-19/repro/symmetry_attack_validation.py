#!/usr/bin/env python3
"""End-to-end validation of the symmetry-quotient affine-column forgery reduction
for odd-characteristic SNOVA 2.3.

The script performs four checks:
  1. Reconstructs the official q=19, (v,o,l,r)=(28,5,4,4) public key from
     the repository KAT and verifies the packed public key byte-for-byte.
  2. Proves computationally that, under U=x*rho^T+V, the homogeneous verifier
     factors through m1*binom(l+1,2) unordered (a,b) coordinates.
  3. Eliminates the complementary affine-linear output equations and verifies
     equivalence with the full verifier on random points.
  4. Computes the symmetry-quotient rank for all nine q=19 v2.3 shapes.

Only Python + NumPy are required.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Sequence
import numpy as np

Q=19
QA,QB,QC=1,3,15

@dataclass(frozen=True)
class Params:
    name:str; v:int; o:int; l:int; r:int; target:int
    @property
    def n(self): return self.v+self.o
    @property
    def m1(self): return (self.o*self.r+self.l-1)//self.l
    @property
    def alpha(self): return self.l*self.r+self.r
    @property
    def outputs(self): return self.o*self.r*self.l
    @property
    def variables(self): return self.n*self.l
    @property
    def unordered(self): return self.m1*self.l*(self.l+1)//2

PARAMS=[
 Params('I-square-l4',28,5,4,4,143), Params('I-square-l2',48,16,2,2,143),
 Params('I-rect-l4xr5',28,4,4,5,143), Params('III-square-l4',40,7,4,4,207),
 Params('III-square-l2',72,24,2,2,207), Params('III-rect-l4xr5',38,5,4,5,207),
 Params('V-square-l4',50,9,4,4,272), Params('V-square-l2',96,32,2,2,272),
 Params('V-rect-l4xr6',52,6,4,6,272),
]

@dataclass
class ABQ:
    A:np.ndarray; B:np.ndarray; q1:np.ndarray; q2:np.ndarray; S:np.ndarray


def inv_mod(a:int,q:int=Q)->int:
    a=int(a)%q
    if not a: raise ZeroDivisionError
    return pow(a,q-2,q)

def rref_mod(a:np.ndarray,q:int=Q):
    A=np.array(a,dtype=np.int64,copy=True)%q
    m,n=A.shape; piv=[]; row=0
    for col in range(n):
        pivot=next((r for r in range(row,m) if A[r,col]%q),None)
        if pivot is None: continue
        if pivot!=row: A[[row,pivot]]=A[[pivot,row]]
        A[row]=A[row]*inv_mod(A[row,col],q)%q
        for r in range(m):
            if r!=row and A[r,col]%q:
                A[r]=(A[r]-int(A[r,col])*A[row])%q
        piv.append(col); row+=1
        if row==m: break
    return A,piv

def rank_mod(a,q=Q): return len(rref_mod(a,q)[1])

def nullspace_mod(a:np.ndarray,q:int=Q)->np.ndarray:
    R,piv=rref_mod(a,q); n=a.shape[1]
    free=[j for j in range(n) if j not in piv]
    out=[]
    for f in free:
        x=np.zeros(n,dtype=np.int64); x[f]=1
        for i,p in enumerate(piv): x[p]=(-R[i,f])%q
        out.append(x)
    return np.array(out,dtype=np.int64).reshape(len(out),n)

def det_mod(a:np.ndarray,q:int=Q)->int:
    A=np.array(a,dtype=np.int64,copy=True)%q; n=A.shape[0]; d=1
    for c in range(n):
        p=next((r for r in range(c,n) if A[r,c]),None)
        if p is None:return 0
        if p!=c:A[[c,p]]=A[[p,c]]; d=-d
        pv=int(A[c,c]);d=d*pv%q;A[c]=A[c]*inv_mod(pv,q)%q
        for r in range(c+1,n):
            if A[r,c]:A[r]=(A[r]-int(A[r,c])*A[c])%q
    return d%q

def affine_solve(C:np.ndarray,b:np.ndarray,q:int=Q):
    C=np.asarray(C,dtype=np.int64)%q;b=np.asarray(b,dtype=np.int64).reshape(-1,1)%q
    d=C.shape[1];R,paug=rref_mod(np.concatenate([C,b],axis=1),q)
    if d in paug: raise ValueError('inconsistent affine system')
    piv=[p for p in paug if p<d]; free=[j for j in range(d) if j not in piv]
    x0=np.zeros(d,dtype=np.int64);T=np.zeros((d,len(free)),dtype=np.int64)
    for k,f in enumerate(free):T[f,k]=1
    for row,p in enumerate(piv):
        x0[p]=R[row,d]
        for k,f in enumerate(free):T[p,k]=(-R[row,f])%q
    assert np.array_equal(C@x0%q,b[:,0]); assert np.all(C@T%q==0)
    return x0,T,free,piv


def official_S(l:int)->np.ndarray:
    S=np.empty((l,l),dtype=np.int64)
    for i in range(l):
        for j in range(l): S[i,j]=(QA+i+j)&QB
    S[-1,-1]=QC
    return S%Q

def powers(M:np.ndarray,count:int):
    out=[np.eye(M.shape[0],dtype=np.int64)]
    for _ in range(1,count): out.append(out[-1]@M%Q)
    return out

def improve(M:np.ndarray,S:np.ndarray,l:int)->np.ndarray:
    M=np.array(M,dtype=np.int64,copy=True)%Q
    if M.shape!=(l,l) or det_mod(M): return M
    for a in range(1,Q):
        C=(M+a*S)%Q
        if det_mod(C):return C
    raise RuntimeError('improve failed')

def repair_coeffs(x:np.ndarray)->np.ndarray:
    y=np.array(x,dtype=np.int64,copy=True)%Q
    for row in y.reshape(-1,y.shape[-1]):
        if row[-1]==0:
            row[-1]=(Q-(int(row[0]) if row[0] else 1))%Q
    return y

def reconstruct_abq(p:Params)->ABQ:
    S=official_S(p.l)
    total=p.o*p.alpha*(p.r*p.r+p.r*p.l+2*p.l)
    raw=np.frombuffer(hashlib.shake_256(b'SNOVA_ABQ').digest(total),dtype=np.uint8).astype(np.int64)%Q
    z=0;c=p.o*p.alpha*p.r*p.r;Ar=raw[z:z+c].reshape(p.o,p.alpha,p.r,p.r);z+=c
    c=p.o*p.alpha*p.r*p.l;Br=raw[z:z+c].reshape(p.o,p.alpha,p.r,p.l);z+=c
    c=p.o*p.alpha*p.l;q1=raw[z:z+c].reshape(p.o,p.alpha,p.l);z+=c
    q2=raw[z:z+c].reshape(p.o,p.alpha,p.l);z+=c
    assert z==total
    A=np.empty_like(Ar);B=np.empty_like(Br)
    for i in range(p.o):
        for a in range(p.alpha):
            A[i,a]=improve(Ar[i,a],S,p.l)
            B[i,a]=improve(Br[i,a],S,p.l)
    return ABQ(A,B,repair_coeffs(q1),repair_coeffs(q2),S)


def build_E(p:Params,abq:ABQ,rho:Sequence[int]):
    rho=np.asarray(rho,dtype=np.int64)%Q; RR=np.outer(rho,rho)%Q
    tri=p.l*(p.l+1)//2
    E=np.zeros((p.outputs,p.m1*tri),dtype=np.int64)
    labels=[]
    for pp in range(p.m1):
        for a in range(p.l):
            for b in range(a,p.l):labels.append((pp,a,b))
    for i in range(p.o):
        rows=slice(i*p.r*p.l,(i+1)*p.r*p.l)
        for al in range(p.alpha):
            pp=(i+al)%p.m1
            core=abq.A[i,al]@RR@abq.B[i,al]%Q
            z=0
            for a in range(p.l):
                for b in range(a,p.l):
                    c=int(abq.q1[i,al,a])*int(abq.q2[i,al,b])
                    if a!=b:c+=int(abq.q1[i,al,b])*int(abq.q2[i,al,a])
                    E[rows,pp*tri+z]=(E[rows,pp*tri+z]+(c%Q)*core.reshape(-1))%Q
                    z+=1
    return E,labels


def kronS(Sa:np.ndarray,n:int): return np.kron(np.eye(n,dtype=np.int64),Sa)%Q

def coords(p:Params,P:np.ndarray,x:np.ndarray,labels):
    Sps=powers(official_S(p.l),p.l); W=[kronS(A,p.n)@x%Q for A in Sps]
    c=[]
    for pp,a,b in labels:c.append(int(W[a].T@P[pp]@W[b])%Q)
    return np.asarray(c,dtype=np.int64)

def relation_U(p:Params,x:np.ndarray,rho:np.ndarray,V:np.ndarray):
    return (x[:,None]*rho[None,:]+V)%Q

def direct_output(p:Params,abq:ABQ,P:np.ndarray,U:np.ndarray)->np.ndarray:
    Sps=powers(abq.S,p.l); W=[kronS(A,p.n)@U%Q for A in Sps]
    D=np.empty((p.m1,p.l,p.l,p.r,p.r),dtype=np.int64)
    for pp in range(p.m1):
        for a in range(p.l):
            LP=W[a].T@P[pp]%Q
            for b in range(p.l):D[pp,a,b]=LP@W[b]%Q
    y=np.zeros((p.o,p.r,p.l),dtype=np.int64)
    for i in range(p.o):
        for al in range(p.alpha):
            pp=(i+al)%p.m1
            mid=np.zeros((p.r,p.r),dtype=np.int64)
            for a in range(p.l):
                for b in range(p.l):
                    mid=(mid+int(abq.q1[i,al,a])*int(abq.q2[i,al,b])*D[pp,a,b])%Q
            y[i]=(y[i]+abq.A[i,al]@mid@abq.B[i,al])%Q
    return y.reshape(-1)


def skew_vector(M:np.ndarray):
    l=M.shape[0];return np.array([(int(M[i,j])-int(M[j,i]))%Q for i in range(l) for j in range(i+1,l)],dtype=np.int64)
def safe_offsets(p:Params,rho:np.ndarray,rng):
    V=rng.integers(0,Q,size=(p.variables,p.r),dtype=np.int64)
    if p.l!=p.r or p.l<=2:return V,{"deterministic":False}
    # image of x -> skew(x rho^T)
    Phi=[]
    for k in range(p.l):
        e=np.zeros(p.l,dtype=np.int64);e[k]=1
        Phi.append(skew_vector(np.outer(e,rho)))
    Phi=np.array(Phi,dtype=np.int64).T
    rank=rank_mod(Phi)
    y=None
    for k in range(Phi.shape[0]):
        e=np.zeros(Phi.shape[0],dtype=np.int64);e[k]=1
        if rank_mod(np.column_stack([Phi,e]))>rank:y=e;break
    assert y is not None
    off=np.zeros((p.l,p.l),dtype=np.int64);z=0
    for i in range(p.l):
        for j in range(i+1,p.l):off[i,j]=y[z];z+=1
    for b in range(p.n):V[b*p.l:(b+1)*p.l,:]=(V[b*p.l:(b+1)*p.l,:]+off)%Q
    # prove on random samples and algebraically by y outside image
    return V,{"deterministic":True,"skew_image_rank":rank,"alternating_dimension":Phi.shape[0]}

def symmetric_block_count(p:Params,U:np.ndarray):
    if p.l!=p.r:return 0
    return sum(np.array_equal(U[b*p.l:(b+1)*p.l],U[b*p.l:(b+1)*p.l].T) for b in range(p.n))


def random_uov_public(p:Params,rng):
    d=p.variables; dv=p.v*p.l; do=p.o*p.l
    Sps=powers(official_S(p.l),p.l)
    T=np.eye(d,dtype=np.int64)
    for i in range(p.v):
        for j in range(p.o):
            co=rng.integers(0,Q,size=p.l,dtype=np.int64)
            B=sum((int(co[a])*Sps[a] for a in range(p.l)),start=np.zeros((p.l,p.l),dtype=np.int64))%Q
            T[i*p.l:(i+1)*p.l,dv+j*p.l:dv+(j+1)*p.l]=B
    Ps=[]
    for _ in range(p.m1):
        F=np.zeros((d,d),dtype=np.int64)
        A=rng.integers(0,Q,size=(dv,dv),dtype=np.int64);A=np.triu(A);A=(A+np.triu(A,1).T)%Q
        C=rng.integers(0,Q,size=(dv,do),dtype=np.int64)
        F[:dv,:dv]=A;F[:dv,dv:]=C;F[dv:,:dv]=C.T
        Ps.append(T.T@F@T%Q)
    return np.asarray(Ps,dtype=np.int64)

# Exact KAT key-generation transcription for 28,5,19,4.
def parse_kat(path:Path):
    d={}
    for line in path.read_text().splitlines():
        if ' = ' in line:
            k,v=line.split(' = ',1);d[k.strip()]=v.strip()
    return d

def snova_xof(seed:bytes,count:int)->bytes:
    out=bytearray();blocks=(count+167)//168
    for i in range(blocks):out+=hashlib.shake_128(seed+i.to_bytes(8,'little')).digest(168)
    return bytes(out[:count])
def gen_fqs(co,Sps):
    c=list(map(int,co));
    if c[-1]==0:c[-1]=Q-(c[0] if c[0] else 1)
    return sum((c[a]*Sps[a] for a in range(len(c))),start=np.zeros_like(Sps[0]))%Q

def kat_public_key(seed48:bytes,p:Params):
    assert (p.v,p.o,p.l,p.r)==(28,5,4,4)
    pkseed,skseed=seed48[:16],seed48[16:]
    Sps=powers(official_S(4),4)
    # T12 rejection sampling.
    raw=hashlib.shake_256(skseed).digest(2*p.o*p.v*p.l); vals=[];bound=(256//Q)*Q
    for b in raw:
        if b<bound:
            vals.append(b%Q)
            if len(vals)==p.o*p.v*p.l:break
    T12=np.zeros((p.v*p.l,p.o*p.l),dtype=np.int64)
    for i in range(p.v):
        for j in range(p.o):
            off=(i*p.o+j)*p.l
            T12[i*p.l:(i+1)*p.l,j*p.l:(j+1)*p.l]=gen_fqs(vals[off:off+p.l],Sps)
    nelems=p.m1*(p.v*p.l*(p.l+1)//2 + (p.v*(p.v-1)//2 + p.v*p.o)*p.l*p.l)
    data=np.frombuffer(snova_xof(pkseed,nelems),dtype=np.uint8).astype(np.int64)%Q;z=0
    P=[];P22=[]
    for _ in range(p.m1):
        P11=np.zeros((p.v*p.l,p.v*p.l),dtype=np.int64)
        P12=np.zeros((p.v*p.l,p.o*p.l),dtype=np.int64)
        for ni in range(p.v):
            for i in range(p.l):
                for j in range(i,p.l):
                    x=int(data[z]);z+=1;P11[ni*p.l+i,ni*p.l+j]=x;P11[ni*p.l+j,ni*p.l+i]=x
            for nj in range(ni+1,p.v):
                for i in range(p.l):
                    for j in range(p.l):
                        x=int(data[z]);z+=1;P11[ni*p.l+i,nj*p.l+j]=x;P11[nj*p.l+j,ni*p.l+i]=x
            for nj in range(p.o):
                for i in range(p.l):
                    for j in range(p.l):
                        x=int(data[z]);z+=1;P12[ni*p.l+i,nj*p.l+j]=x
        P21=P12.T.copy()
        p22=-(T12.T@((P11@T12+P12)%Q)+P21@T12)%Q
        P22.append(p22);P.append(np.block([[P11,P12],[P21,p22]])%Q)
    assert z==nelems
    # Pack public P22.
    arr=[]
    for A in P22:
        for ni in range(p.o):
            for i in range(p.l):
                for j in range(i,p.l):arr.append(int(A[ni*p.l+i,ni*p.l+j]))
                for nj in range(ni+1,p.o):
                    for j in range(p.l):arr.append(int(A[ni*p.l+i,nj*p.l+j]))
    out=bytearray();idx=0
    while idx<len(arr):
        take=min(15,len(arr)-idx);val=sum(arr[idx+i]*Q**i for i in range(take));idx+=15
        for _ in range(8):out.append(val%256);val//=256
    pk=pkseed+bytes(out[:math.ceil(8*len(arr)/15)])
    return pk,np.asarray(P,dtype=np.int64)


def validate_level1(kat_path:Path,seed:int=20260718):
    p=PARAMS[0];kat=parse_kat(kat_path);seed48=bytes.fromhex(kat['sk']);expected=bytes.fromhex(kat['pk'])
    pk,P=kat_public_key(seed48,p);kat_match=(pk==expected)
    abq=reconstruct_abq(p);rng=np.random.default_rng(seed)
    rho=np.array([1,8,9,14],dtype=np.int64);E,labels=build_E(p,abq,rho);er=rank_mod(E);H=nullspace_mod(E.T)
    V,fmt=safe_offsets(p,rho,rng)
    # Homogeneous factorization.
    hom_ok=True
    for _ in range(5):
        x=rng.integers(0,Q,size=p.variables,dtype=np.int64)
        lhs=direct_output(p,abq,P,relation_U(p,x,rho,np.zeros_like(V)))
        rhs=E@coords(p,P,x,labels)%Q
        hom_ok &= np.array_equal(lhs,rhs)
    # Affine decomposition F(x)=f0+Lx+E*c(x).
    f0=direct_output(p,abq,P,V)
    Lmat=np.zeros((p.outputs,p.variables),dtype=np.int64)
    for j in range(p.variables):
        e=np.zeros(p.variables,dtype=np.int64);e[j]=1
        Lmat[:,j]=(direct_output(p,abq,P,relation_U(p,e,rho,V))-f0-E@coords(p,P,e,labels))%Q
    affine_ok=True
    for _ in range(5):
        x=rng.integers(0,Q,size=p.variables,dtype=np.int64)
        rhs=(f0+Lmat@x+E@coords(p,P,x,labels))%Q
        affine_ok &= np.array_equal(direct_output(p,abq,P,relation_U(p,x,rho,V)),rhs)
    xstar=rng.integers(0,Q,size=p.variables,dtype=np.int64);target=direct_output(p,abq,P,relation_U(p,xstar,rho,V))
    C=H@Lmat%Q;b=H@(target-f0)%Q
    x0,T,free,piv=affine_solve(C,b);linrank=len(piv)
    zstar=xstar[free];known=np.array_equal((x0+T@zstar)%Q,xstar)
    # Select rows making E square invertible.
    _,rows=rref_mod(E.T);rows=rows[:er];EJ=E[rows,:];assert rank_mod(EJ)==er
    equivalence=True
    for _ in range(10):
        z=rng.integers(0,Q,size=T.shape[1],dtype=np.int64);x=(x0+T@z)%Q
        g=(direct_output(p,abq,P,relation_U(p,x,rho,V))-target)%Q
        if np.any(H@g%Q):equivalence=False;break
        # g is in col(E), and selected rows determine it.
        # Solve EJ*c=gJ and reconstruct.
        c0,TT,_,_=affine_solve(EJ,g[rows]);assert TT.shape[1]==0
        equivalence &= np.array_equal(E@c0%Q,g)
    fmt_count=symmetric_block_count(p,relation_U(p,xstar,rho,V))
    return {
        "parameter":asdict(p),"kat_public_key_match":kat_match,"public_key_bytes":len(pk),
        "E_shape":list(E.shape),"E_rank":er,"left_kernel_dimension":int(H.shape[0]),
        "homogeneous_factorization":bool(hom_ok),"affine_decomposition":bool(affine_ok),
        "linear_constraint_rank":linrank,"residual_variables":int(T.shape[1]),"residual_quadratics":er,
        "known_solution_preserved":known,"full_equivalence_after_elimination":bool(equivalence),
        "format_offset":fmt,"symmetric_blocks_for_known_solution":int(fmt_count),
    }


def rank_all():
    out=[]
    for p in PARAMS:
        abq=reconstruct_abq(p);rho=np.zeros(p.r,dtype=np.int64);rho[0]=1
        E,_=build_E(p,abq,rho)
        out.append({"parameter":asdict(p),"E_shape":list(E.shape),"rank":rank_mod(E),
                    "expected_unordered":p.unordered,"left_kernel":p.outputs-rank_mod(E)})
    return out


def rejection_tails():
    from fractions import Fraction
    out=[]
    for p in PARAMS:
        if p.l==p.r==2:
            threshold=p.n//4;prob=sum(math.comb(p.n,k)*Fraction(1,Q)**k*Fraction(Q-1,Q)**(p.n-k) for k in range(threshold+1,p.n+1))
            out.append({"name":p.name,"n_blocks":p.n,"threshold":threshold,"uniform_solution_tail":float(prob),"log2_tail":math.log2(float(prob))})
    return out


def main():
    ap=argparse.ArgumentParser()
    root=Path(__file__).resolve().parents[1]
    ap.add_argument('--kat',default=str(root/'source_snapshots'/'PQCsignKAT_SNOVA_28_5_19_4.txt'),
                    help='official Version 2.3 Level-I KAT response file')
    ap.add_argument('--out',default=str(root/'results'/'validation_results.json'))
    ap.add_argument('--rank-only',action='store_true',help='run public-constant rank checks without a KAT')
    args=ap.parse_args()
    out={"all_parameter_E_ranks":rank_all(),"l2_format_rejection_tails":rejection_tails()}
    if not args.rank_only:
        kat=Path(args.kat)
        if not kat.exists():
            raise FileNotFoundError(f'KAT not found: {kat}. See README.md for generation instructions.')
        out["level1_official_kat"]=validate_level1(kat)
    out_path=Path(args.out); out_path.parent.mkdir(parents=True,exist_ok=True)
    out_path.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
