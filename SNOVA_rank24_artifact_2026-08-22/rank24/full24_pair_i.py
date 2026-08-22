from pathlib import Path
import numpy as np,sys,contextlib,io,time,json
BASE=Path(__file__).resolve().parent;sys.path.insert(0,str(BASE));import block_rank_gf361 as G
chart=int(sys.argv[1]);i=int(sys.argv[2]);z=np.load(BASE/f'full24_prep_e{chart}.npz');Cs=z['Cs'];s=len(Cs);bad=[];t=time.time();checks=0
C=Cs[i]
try: piv=G.pivot_cols(C,300)
except Exception as e:
 print(json.dumps({'chart':chart,'i':i,'checks':1,'bad':[('single',i,str(e))],'seconds':time.time()-t}));raise SystemExit(1)
checks+=1; pset=set(piv);rest=[j for j in range(C.shape[1]) if j not in pset];Pinv=G.invmat(C[:,piv]);R=G.gfmm(Pinv,C[:,rest])
for j in range(i+1,s):
 L=G.gfmm(Cs[j][:,piv],R);Sch=G.ADD[Cs[j][:,rest],G.NEG[L]]
 with contextlib.redirect_stdout(io.StringIO()):r=G.block_rank(Sch,256)
 checks+=1
 if r!=300:bad.append(('pair',i,j,int(r)))
print(json.dumps({'chart':chart,'i':i,'checks':checks,'bad':bad,'seconds':time.time()-t}))
raise SystemExit(1 if bad else 0)
