import numpy as np, time, sys
P=19;Q=361
xs=np.arange(Q);aa=xs%P;bb=xs//P
ADD=((aa[:,None]+aa[None,:])%P + P*((bb[:,None]+bb[None,:])%P)).astype(np.int16)
MUL=((aa[:,None]*aa[None,:]-bb[:,None]*bb[None,:])%P + P*((aa[:,None]*bb[None,:]+bb[:,None]*aa[None,:])%P)).astype(np.int16)
NEG=np.array([((-x%P)%P)+P*((-(x//P))%P) for x in range(Q)],dtype=np.int16)
INV=np.zeros(Q,dtype=np.int16)
for x in range(1,Q):
 a,b=x%P,x//P;den=(a*a+b*b)%P;di=pow(int(den),-1,P);INV[x]=(a*di%P)+P*((-b*di)%P)
def gfmm(A,B):
 Ac=(A%P).astype(np.float64)+1j*(A//P).astype(np.float64)
 Bc=(B%P).astype(np.float64)+1j*(B//P).astype(np.float64)
 C=Ac@Bc
 re=np.rint(C.real).astype(np.int64)%P; im=np.rint(C.imag).astype(np.int64)%P
 return (re+P*im).astype(np.int16)
def invmat(A):
 A=A.copy();n=A.shape[0];I=np.zeros((n,n),dtype=np.int16);np.fill_diagonal(I,1);r=0
 for c in range(n):
  nz=np.flatnonzero(A[r:,c]);
  if not len(nz): raise ValueError('sing')
  p=r+int(nz[0]);
  if p!=r:A[[r,p]]=A[[p,r]];I[[r,p]]=I[[p,r]]
  iv=INV[A[r,c]];A[r]=MUL[iv,A[r]];I[r]=MUL[iv,I[r]]
  inds=np.flatnonzero(A[:,c]);inds=inds[inds!=r]
  for lo in range(0,len(inds),128):
   jj=inds[lo:lo+128];ff=A[jj,c].copy(); A[jj]=ADD[A[jj],NEG[MUL[ff[:,None],A[r][None,:]]]]; I[jj]=ADD[I[jj],NEG[MUL[ff[:,None],I[r][None,:]]]]
  r+=1
 return I
def pivot_cols(T,b):
 # T first b rows x m, return b independent columns
 A=T.copy();nr,nc=A.shape;r=0;pivs=[]
 for c in range(nc):
  nz=np.flatnonzero(A[r:,c]);
  if not len(nz):continue
  p=r+int(nz[0]);
  if p!=r:A[[r,p]]=A[[p,r]]
  iv=INV[A[r,c]];A[r]=MUL[iv,A[r]]
  inds=np.flatnonzero(A[r+1:,c])+r+1
  if len(inds):
   ff=A[inds,c].copy();A[inds]=ADD[A[inds],NEG[MUL[ff[:,None],A[r][None,:]]]]
  pivs.append(c);r+=1
  if r==b:return pivs
 raise ValueError(f'row rank only {r}/{b}')
def block_rank(M,bs=256):
 # rectangular nr x nc, prove full row rank by Schur elimination
 M=M.copy();rank=0;t0=time.time()
 while M.shape[0]:
  nr,nc=M.shape;b=min(bs,nr)
  piv=pivot_cols(M[:b,:],b)
  rest=[j for j in range(nc) if j not in set(piv)]
  perm=piv+rest;M=M[:,perm]
  Piv=M[:b,:b];Pinv=invmat(Piv)
  if nr==b:
   rank+=b;print('done rank',rank,'time',time.time()-t0,flush=True);return rank
  B=M[:b,b:];C=M[b:,:b];D=M[b:,b:]
  L=gfmm(C,Pinv)
  # chunked Schur update to cap complex-temporary memory
  New=np.empty_like(D); chunk=512
  for jj in range(0,B.shape[1],chunk):
   prod=gfmm(L,B[:,jj:jj+chunk]); New[:,jj:jj+chunk]=ADD[D[:,jj:jj+chunk],NEG[prod]]
  M=New
  rank+=b
  print('rank',rank,'remain',M.shape,'density',np.count_nonzero(M)/M.size,'time',time.time()-t0,flush=True)
 return rank
if __name__=='__main__':
 M=np.load(sys.argv[1]);bs=int(sys.argv[2]) if len(sys.argv)>2 else 256
 print('input',M.shape,'density',np.count_nonzero(M)/M.size,flush=True)
 print('RANK',block_rank(M,bs),'/',M.shape[0])
