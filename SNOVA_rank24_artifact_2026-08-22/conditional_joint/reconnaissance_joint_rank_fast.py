import numpy as np,time
from numba import njit
q=19; BAS=np.load('/mnt/data/kat0_cert/joint_pencil_basis.npz')['BAS'].astype(np.int16)
@njit(cache=True)
def inv(a):
 for x in range(1,19):
  if (a*x)%19==1:return x
 return 0
@njit(cache=True)
def rankm(A):
 M=A.copy();nr,nc=M.shape;r=0
 for c in range(nc):
  p=-1
  for i in range(r,nr):
   if M[i,c]%19!=0:p=i;break
  if p<0:continue
  if p!=r:
   for j in range(c,nc):t=M[r,j];M[r,j]=M[p,j];M[p,j]=t
  iv=inv(M[r,c]%19)
  for j in range(c,nc):M[r,j]=(M[r,j]*iv)%19
  for i in range(r+1,nr):
   f=M[i,c]%19
   if f:
    for j in range(c,nc):M[i,j]=(M[i,j]-f*M[r,j])%19
  r+=1
  if r==nc:return r
 return r
@njit(cache=True)
def combo(c,B):
 nr=B.shape[1];nc=B.shape[2];M=np.zeros((nr,nc),np.int16)
 for k in range(B.shape[0]):
  ck=int(c[k])
  if ck:
   for i in range(nr):
    for j in range(nc):M[i,j]=(M[i,j]+ck*B[k,i,j])%19
 return M
# compile
c=np.zeros(48,np.int16);c[0]=1;print('compile rank',rankm(combo(c,BAS)))
rng=np.random.default_rng(3)
for typ,N in [('D',20000),('Hactive',100000)]:
 hist={};mn=99;ex=None;t=time.time()
 for z in range(N):
  c=rng.integers(0,19,size=48,dtype=np.int16)
  if typ=='D':c[32:]=0
  elif not np.any(c[32:]):c[32]=1
  r=rankm(combo(c,BAS));hist[r]=hist.get(r,0)+1
  if r<mn:mn=r;ex=c.copy();print(typ,'newmin',r,z,flush=True)
 print(typ,'hist',hist,'min',mn,'sec',time.time()-t,flush=True)
 np.savez_compressed('/mnt/data/kat0_cert/joint_probe_'+typ+'.npz',example=ex,hist_keys=np.array(list(hist)),hist_vals=np.array(list(hist.values())))
