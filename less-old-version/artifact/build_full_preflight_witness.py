#!/usr/bin/env python3
"""Exact disjoint-coordinate witnesses for all structured public rank preflights.

For each official ell=4 shape, this script verifies:
  1. the relation-only symmetric quadratic expansion X=R_R D has rank K;
  2. the source matrix [X F], restricted to all K quotient columns and selected
     fresh V' columns, has rank M;
  3. U_w restricted to J=W+O_pub has full column rank.

The source-rank minor depends only on W x V' public-form coefficients.  The
projection-injectivity minor depends only on W x (W+O_pub) coefficients.
These coordinate sets are disjoint, and both are disjoint from the internal
V' x V' coefficients used by the eigenblock-core witness.  Exact blockwise
nonuniform Schwartz-Zippel bounds are computed with rho=14/256.

We use R_j=S^j, Gamma_0=1, Gamma_1=0, Gamma_j=S^(j-2) for j>=2.  The (1,0)
raw column pair then contains every row w^T Sigma_xi P_i Sigma_eta, proving
row(F)=U_w.
"""
from __future__ import annotations
import json, math, random, sys
sys.set_int_max_str_digits(100000)
from fractions import Fraction
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec=100
ROOT=Path(__file__).resolve().parent
Q=19; D=4; RHO=Fraction(14,256)
S=[[1,2,3,0],[2,3,0,1],[3,0,1,2],[0,1,2,15]]
SHAPES=[
 {"key":"I-a","level":"I","v":28,"o":5,"r":4,"m1":5,"M":80,"K":50},
 {"key":"I-b","level":"I","v":28,"o":4,"r":5,"m1":5,"M":80,"K":50},
 {"key":"III-a","level":"III","v":40,"o":7,"r":4,"m1":7,"M":112,"K":70},
 {"key":"III-b","level":"III","v":38,"o":5,"r":5,"m1":7,"M":100,"K":70},
 {"key":"V-a","level":"V","v":50,"o":9,"r":4,"m1":9,"M":144,"K":90},
 {"key":"V-b","level":"V","v":52,"o":6,"r":6,"m1":9,"M":144,"K":90},
]

def eye(n):return [[1 if i==j else 0 for j in range(n)] for i in range(n)]
def zero(n):return [[0]*n for _ in range(n)]
def mm(A,B):return [[sum(A[i][k]*B[k][j] for k in range(len(B)))%Q for j in range(len(B[0]))] for i in range(len(A))]
POW=[eye(D)]
for _ in range(11):POW.append(mm(POW[-1],S))

def solve_coords(A):
 # express 4x4 A in basis I,S,S2,S3
 B=POW[:4]; aug=[]
 for i in range(D):
  for j in range(D):aug.append([B[k][i][j] for k in range(D)]+[A[i][j]])
 r=0
 for c in range(D):
  p=next(i for i in range(r,len(aug)) if aug[i][c]%Q)
  aug[r],aug[p]=aug[p],aug[r]; inv=pow(aug[r][c]%Q,-1,Q);aug[r]=[(x*inv)%Q for x in aug[r]]
  for i in range(len(aug)):
   if i!=r and aug[i][c]%Q:
    f=aug[i][c]%Q;aug[i]=[(x-f*y)%Q for x,y in zip(aug[i],aug[r])]
  r+=1
 return [aug[i][-1]%Q for i in range(D)]
COORD={e:solve_coords(POW[e]) for e in range(len(POW))}
PAIRS=[(a,b) for a in range(D) for b in range(a,D)]

def sym_feature(x,y):
 out=[]
 for a,b in PAIRS:
  if a==b:out.append(x[a]*y[b]%Q)
  else:out.append((x[a]*y[b]+x[b]*y[a])%Q)
 return out

def rand_B(rng,nblocks):
 n=D*nblocks; B=[[rng.randrange(Q) for _ in range(n)] for _ in range(D)]
 for i in range(D):
  for j in range(i,D):
   x=rng.randrange(Q);B[i][j]=B[j][i]=x
 return B

def left_offset_row(B,G,L,nblocks):
 # (G e0)^T B diag(L)
 v=[G[i][0] for i in range(D)]
 mid=[sum(v[i]*B[i][j] for i in range(D))%Q for j in range(D*nblocks)]
 out=[0]*(D*nblocks)
 for blk in range(nblocks):
  z=blk*D
  for j in range(D):out[z+j]=sum(mid[z+k]*L[k][j] for k in range(D))%Q
 return out

def addrow(a,b):return [(x+y)%Q for x,y in zip(a,b)]

def rank_select(M,ncols,limit=None):
 A=[row[:] for row in M]; labels=list(range(len(A))); r=0; piv=[]; sel=[]; det=1
 for c in range(ncols):
  p=next((i for i in range(r,len(A)) if A[i][c]%Q),None)
  if p is None:continue
  A[r],A[p]=A[p],A[r];labels[r],labels[p]=labels[p],labels[r]
  pv=A[r][c]%Q;det=det*pv%Q;inv=pow(pv,-1,Q);A[r]=[(x*inv)%Q for x in A[r]]
  for i in range(r+1,len(A)):
   if A[i][c]%Q:
    f=A[i][c]%Q;A[i]=[(x-f*y)%Q for x,y in zip(A[i],A[r])]
  piv.append(c);sel.append(labels[r]);r+=1
  if limit and r==limit:break
 return r,sel,piv,det

def beta_from_degrees(deg):
 b=Fraction(1,1)
 for d in deg:
  f=1-d*RHO
  if f<=0:return Fraction(0,1)
  b*=f
 return b

def log2f(x):return float((Decimal(x.numerator).ln()-Decimal(x.denominator).ln())/Decimal(2).ln())

def make_rows(sh,Bs):
 v,o,r,m1=sh['v'],sh['o'],sh['r'],sh['m1'];n=v+o;K=sh['K']
 R=[POW[j] for j in range(r)]
 Gamma=[eye(D),zero(D)]+[POW[j-2] for j in range(2,r)]
 rows=[];meta=[]
 # Round-robin forms at each raw label/column tuple.
 for a in range(D):
  for b in range(D):
   for j in range(r):
    for k in range(r):
     LA=mm(POW[a],R[j]); LB=mm(POW[b],R[k])
     GA=mm(POW[a],Gamma[j]); GB=mm(POW[b],Gamma[k])
     x=sym_feature(solve_coords(LA),solve_coords(LB))
     for i,B in enumerate(Bs):
      xfull=[0]*K;xfull[10*i:10*(i+1)]=x
      f=addrow(left_offset_row(B,GB,LA,n),left_offset_row(B,GA,LB,n))
      # Source-rank witness deliberately exposes only V' columns (blocks 1..v-1).
      fv=f[D:D*v]
      rows.append(xfull+fv);meta.append((i,a,b,j,k))
 return rows,meta

def make_U_rows(sh,Bs):
 v,o,m1=sh['v'],sh['o'],sh['m1'];n=v+o
 Jblocks=[0]+list(range(v,n)); cols=[]
 for blk in Jblocks:cols.extend(range(D*blk,D*(blk+1)))
 rows=[];meta=[]
 for a in range(D):
  for b in range(D):
   for i,B in enumerate(Bs):
    full=left_offset_row(B,POW[a],POW[b],n)
    rows.append([full[c] for c in cols]);meta.append((i,a,b))
 return rows,meta

def certify(sh):
 n=sh['v']+sh['o']; m1=sh['m1']; M=sh['M'];K=sh['K']; source_cols=K+D*(sh['v']-1)
 for seed in range(1000):
  rng=random.Random(0x46554c4c ^ seed ^ sum((i+1)*ord(c) for i,c in enumerate(sh['key'])))
  Bs=[rand_B(rng,n) for _ in range(m1)]
  Z,zmeta=make_rows(sh,Bs)
  zr,zsel,zpiv,zdet=rank_select(Z,source_cols,limit=M)
  if zr<M or not all(c in zpiv for c in range(K)):continue
  U,umeta=make_U_rows(sh,Bs); udim=D*(1+sh['o'])
  ur,usel,upiv,udet=rank_select(U,udim,limit=udim)
  if ur<udim:continue
  # X-only rank.
  X=[row[:K] for row in Z]
  xr,_,_,_=rank_select(X,K,limit=K)
  if xr<K:continue
  zchosen=[zmeta[i] for i in zsel]; uchosen=[umeta[i] for i in usel]
  zrows=[0]*m1;urows=[0]*m1
  for i,*_ in zchosen:zrows[i]+=1
  for i,*_ in uchosen:urows[i]+=1
  # All ten X columns of form i can be filled only by rows of form i.  Since
  # every quotient column is in the minor, the determinant has exact P_i degree
  # zrows[i]-10 in the V' coordinate block.
  zdeg=[x-10 for x in zrows]
  if min(zdeg)<0 or sum(zdeg)!=M-K:continue
  udeg=urows
  bz=beta_from_degrees(zdeg);bu=beta_from_degrees(udeg)
  if not bz or not bu:continue
  return {**sh,"seed":seed,"X_rank":xr,"source_rank":zr,"source_minor_det_mod_19":zdet,
    "source_selected_rows":[list(x) for x in zchosen],"source_pivot_columns":zpiv,
    "source_rows_per_form":zrows,"source_degree_vector_WxVprime":zdeg,
    "source_probability_lower_numerator":bz.numerator,"source_probability_lower_denominator":bz.denominator,
    "source_probability_lower_decimal":float(bz),"source_log2_lower":log2f(bz),
    "projection_rank":ur,"projection_minor_det_mod_19":udet,"projection_selected_rows":[list(x) for x in uchosen],
    "projection_rows_per_form":urows,"projection_degree_vector_WxJ":udeg,
    "projection_probability_lower_numerator":bu.numerator,"projection_probability_lower_denominator":bu.denominator,
    "projection_probability_lower_decimal":float(bu),"projection_log2_lower":log2f(bu),
    "joint_disjoint_probability_lower_decimal":float(bz*bu),"joint_disjoint_log2_lower":log2f(bz*bu),
    "W_row_assignments":Bs}
 raise RuntimeError(sh['key'])

def outer_success(dim):
 # Independent rows under any surjective linear image: codimension-c subspace
 # has probability <=rho^c.  Sequential independence gives this product.
 p=Fraction(1,1)
 for c in range(1,dim+1):p*=1-RHO**c
 return p

def main():
 rec=[certify(s) for s in SHAPES]
 pK={k:outer_success(k) for k in (50,70,90)};pM={m:outer_success(m) for m in (80,100,112,144)}
 for x in rec:
  # Since X is a subset of the columns of Z=[X F], invertibility of E Z
  # already implies rank(E X)=K.  No union bound over two rank events is needed.
  joint_outer=pM[x['M']]
  x['outer_joint_probability_lower_numerator']=joint_outer.numerator;x['outer_joint_probability_lower_denominator']=joint_outer.denominator
  x['outer_joint_probability_lower_decimal']=float(joint_outer);x['outer_joint_log2_lower']=log2f(joint_outer)
  cross=Fraction(x['source_probability_lower_numerator'],x['source_probability_lower_denominator'])*Fraction(x['projection_probability_lower_numerator'],x['projection_probability_lower_denominator'])*joint_outer
  x['complete_cross_preflight_probability_lower_numerator']=cross.numerator;x['complete_cross_preflight_probability_lower_denominator']=cross.denominator
  x['complete_cross_preflight_probability_lower_decimal']=float(cross);x['complete_cross_preflight_log2_lower']=log2f(cross)
 data={"field":"F_19","rho":{"numerator":14,"denominator":256},"official_S":S,
  "relation":"R_j=S^j","offsets":"Gamma_0=1, Gamma_1=0, Gamma_j=S^(j-2) for j>=2",
  "disjointness":"source minor uses W x V'; projection minor uses W x (W+O); core witness uses V' x V'.",
  "outer_map_bound":"For a surjective linear image of one independent product-distributed row, membership in codimension c has probability at most rho^c. Sequential rows give full rank of E[X F]. Because X is a subset of the columns of [X F], this single event also gives rank(EX)=K; no union bound is required.",
  "shapes":rec}
 (ROOT/'full_preflight_witness.json').write_text(json.dumps(data,indent=2)+'\n')
 md=["# Exact disjoint-coordinate structured-preflight witnesses","", "Generated by `build_full_preflight_witness.py`; all ranks and determinants are exact modulo 19.","",
 "The source minor uses only `W x V'` coefficients, the projection minor only `W x (W+O_pub)`, and the core-Jacobian witness only `V' x V'`. These byte regions are disjoint. The random outer map is treated separately.","",
 "| shape | rank X | rank [X F] | source degrees | source lower | projection degrees | projection lower | outer lower | complete cross-preflight lower |",
 "|:--:|--:|--:|:--|--:|:--|--:|--:|--:|"]
 for x in rec:
  md.append(f"| {x['key']} | {x['X_rank']} | {x['source_rank']} | `{x['source_degree_vector_WxVprime']}` | {x['source_probability_lower_decimal']:.6f} | `{x['projection_degree_vector_WxJ']}` | {x['projection_probability_lower_decimal']:.6f} | {x['outer_joint_probability_lower_decimal']:.6f} | {x['complete_cross_preflight_probability_lower_decimal']:.8f} ($2^{{{x['complete_cross_preflight_log2_lower']:.3f}}}$) |")
 md += ["", "Full rank of `[X F]` and of the random outer image implies `rank Q_R=K` and `rank C_{R,v}=M-K`. Full rank of `U_w|_{W+O}` implies `H intersect (W+O)=0`, hence projection injectivity on the entire stable kernel. The bounds are information-theoretic in the random-XOF idealization; every rank is checked exactly for a fixed public key."]
 (ROOT/'FULL_PREFLIGHT_WITNESS.md').write_text('\n'.join(md)+'\n')
 for x in rec:print(x['key'],x['X_rank'],x['source_rank'],x['source_degree_vector_WxVprime'],x['projection_degree_vector_WxJ'],f"cross={x['complete_cross_preflight_probability_lower_decimal']:.8f}")
if __name__=='__main__':main()
