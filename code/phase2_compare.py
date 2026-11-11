# Phase 2 versus the prediction made from Phase 1 alone: each searched pass carries the Phase 1 detection probability of
# its class and sigma_min bin (results/blind_summary.json); the sum is the predicted number of detections. For the ISS
# complex the unit is one recording (any member detected), binned by the sigma_min of the station itself (ISS (ZARYA)).
import json, numpy as np
from scipy.stats import beta, poisson
A=json.load(open('../results/blind_summary.json')); B=json.load(open('../results/phase2_summary.json'))
def cp(k,n,a=0.05): return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)), 1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
exec(open('phase2_analyze.py').read().split("th={st:7.0 for st in STS}")[0].split("# echo track of an owner")[0])
det={(d['file'],d['obj']) for d in B['detections']}
for r in P: r['det']=(r['file'],r['obj']) in det
def p1(c,s):
    for q in A['pd_vs_sigma_min'].get(c,[]):
        if q['lo']<=s<q['hi']: return q['k']/q['n'],q['n']
    return None,0
out={}
for c in sorted(set(r['class'] for r in P)):
    if c=='ISS complex': continue
    pp=[r for r in P if r['class']==c]; cov=[r for r in pp if p1(c,r['sigma_min'])[1]>=10]
    pred=sum(p1(c,r['sigma_min'])[0] for r in cov); k=sum(r['det'] for r in cov)
    out[c]=dict(n=len(pp),n_covered=len(cov),k_covered=k,predicted=pred,k_all=sum(r['det'] for r in pp))
# ISS complex: per recording, sigma_min of the station
g={}
for r in P:
    if r['class']=='ISS complex': g.setdefault(r['file'],[]).append(r)
rows=[]
for f,u in g.items():
    st=[r for r in u if r['obj'].startswith('ISS')]
    s=min(r['sigma_min'] for r in (st or u)); rows.append((s,any(r['det'] for r in u),f))
bins=[0,100,300,1e3,3e3,1e4,1e9]; iss=[]
for a,b in zip(bins[:-1],bins[1:]):
    q=[x for x in rows if a<=x[0]<b]
    if q: k=sum(x[1] for x in q); iss.append(dict(lo=a,hi=b,n=len(q),k=k,ci=cp(k,len(q))))
out['ISS complex by station sigma_min']=iss
# Phase 1 ISS units with the same binning are in blind_summary only per member; recompute from Phase 1 records
json.dump(out,open('../results/phase2_compare.json','w'),indent=1)
for c,v in out.items(): print(c,v if not isinstance(v,list) else '')
for x in iss: print('  ISS %g-%g: %d/%d [%.2f,%.2f]'%(x['lo'],x['hi'],x['k'],x['n'],*x['ci']))
