from pathlib import Path
import subprocess,sys,json
B=Path(__file__).resolve().parent
bad=[];checks=0
for c in range(15):
    s=15-c
    for i in range(s):
        p=subprocess.run([sys.executable,str(B/'full24_pair_i.py'),str(c),str(i)],cwd=B,capture_output=True,text=True)
        if not p.stdout.strip():
            bad.append((c,i,'no output',p.stderr));continue
        d=json.loads(p.stdout.strip().splitlines()[-1]);checks+=d['checks'];bad.extend(d['bad'])
        print(c,i,d['checks'],len(d['bad']),flush=True)
print('TOTAL',checks,'BAD',len(bad));
if bad: print(bad[:20]); raise SystemExit(1)
assert checks==680
