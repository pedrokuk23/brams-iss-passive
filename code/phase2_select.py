# Phase 2 file selection by expected yield: each pass contributes its Phase-1 detection probability for its class and
# sigma_min bin; recordings are ranked by the sum per receiver and taken in decreasing order. ISS-complex slots are
# always taken at all receivers. Selection uses only predicted geometry, never the recordings.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import json, collections, numpy as np
P=json.load(open('../results/phase2_plan_passes.json')); B=json.load(open('../results/blind_summary.json'))['pd_vs_sigma_min']
def pd(c,s):
    rows=B.get(c,[])
    for r in rows:
        if r['lo']<=s<r['hi']: return max(r['k']/r['n']-1e-3,0)
    return 0.0
ex=collections.defaultdict(float); iss=set(); ncls=collections.defaultdict(lambda: collections.Counter())
for p in P:
    ex[p['file']]+=pd(p['cls'],p['sigma_min']); ncls[p['file']][p['cls']]+=1
    if p['cls']=='ISS complex': iss.add(p['file'][:13])
rank={}
for st in ('BEOUDS','BEHUMA'):
    fs=sorted([f for f in ex if f.endswith(st)],key=lambda f:-ex[f]); rank[st]=fs
issf=sorted(f for f in ex if f[:13] in iss)
print('ISS-complex slots',len(iss),'files at all receivers',len(issf),'expected det %.1f'%sum(ex[f] for f in issf))
for N in (200,400,600,800,1200,2000,4320):
    sel=set(issf)|set(rank['BEOUDS'][:N])|set(rank['BEHUMA'][:N])
    print('N per receiver %4d: files %5d (%.1f GB) expected detections %.0f'%(N,len(sel),len(sel)*3.66e-3,sum(ex[f] for f in sel)))
json.dump(dict(rank=rank,iss_files=issf,expected={f:ex[f] for f in ex}),open('../results/phase2_rank.json','w'))
