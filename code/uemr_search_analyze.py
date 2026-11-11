# Summary of uemr_search.py: thresholds from the one-way nulls of each receiver (rule of the blind search, quantile 1e-3),
# exceedances along the one-way history, and whether they coincide with bistatic detections of the same pass (a strong
# reflection also raises the one-way statistic near its peak). Split by class and by sunlight/shadow (eclipse flag of the
# Phase 2 units of stats_cluster.py).
import json, glob, numpy as np
from far_checks import rule, cp
R={}
import os
for g in sorted(glob.glob(os.environ.get('UEMRGLOB','../results/uemr_search_[0-9]*.json'))):
    for r in json.load(open(g)):
        if r.get('case') in ('pass','null'): R[(r['case'],r['file'],r['obj'])]=r
R=list(R.values()); inside=lambda r: abs(r['tau'])<6 and abs(r['nu'])<12
S=json.load(open('../results/phase2_summary.json')); DET={(d['file'],d['obj']) for d in S['detections']}
exec(open('stats_cluster.py').read().split("out=dict(n_boot")[0])       # U2m with the sunlight group
GRP={(u['file'],u['r']['obj']):u['grp'] for u in U2m}
files={r['file'] for r in R}; out={'n_files':len(files)}
th={}
for st in sorted({r['file'][-6:] for r in R}):
    z=[r['z'] for r in R if r['case']=='null' and r['file'].endswith(st) and inside(r)]
    n=sum(1 for r in R if r['case']=='null' and r['file'].endswith(st)); th[st]=rule(z,n,1e-3)
out['thresholds']=th
ex=[r for r in R if r['case']=='pass' and inside(r) and r['z']>=th[r['file'][-6:]]]
nul=[r for r in R if r['case']=='null' and inside(r) and r['z']>=th[r['file'][-6:]]]
P=[r for r in R if r['case']=='pass']; N=[r for r in R if r['case']=='null']
out['all']=dict(n_pass=len(P),exceed=len(ex),exceed_also_bistatic=sum((r['file'],r['obj']) in DET for r in ex),
                n_null=len(N),null_exceed=len(nul),bistatic_detections_in_these_passes=sum((r['file'],r['obj']) in DET for r in P))
for c in ('Starlink V2-mini','Starlink DTC'):
    for g in ('plate_shadow','plate_sunlit'):
        p=[r for r in P if r['cls']==c and GRP.get((r['file'],r['obj']))==g]; e=[r for r in ex if r in p]
        out[f'{c}|{g}']=dict(n=len(p),exceed=len(e),exceed_not_bistatic=sum((r['file'],r['obj']) not in DET for r in e),
                             bistatic_detections=sum((r['file'],r['obj']) in DET for r in p),
                             rate_not_bistatic=sum((r['file'],r['obj']) not in DET for r in e)/len(p) if p else None,
                             ci=cp(sum((r['file'],r['obj']) not in DET for r in e),len(p)) if p else None)
# chance exceedances expected with these thresholds: with the rule of far_checks (k = floor(1e-3 n) exceedances allowed among
# n nulls), a pass whose statistic is exchangeable with the nulls exceeds the threshold with probability (k+1)/(n+1);
# exact distribution of the total count by simulation of the ranks (exchangeability), and the Poisson approximation.
from scipy.stats import poisson, chi2, beta as _b
rng=np.random.default_rng(5); E=0.0; sims=np.zeros(20000,int)
for st in th:
    n=sum(1 for r in N if r['file'].endswith(st)); m=sum(1 for r in P if r['file'].endswith(st)); k=int(np.floor(1e-3*n))
    E+=(k+1)*m/(n+1)
    # rank of the (k+1)-th largest null among n+1 slots ~ Beta: exceedance probability q ~ Beta(k+1, n-k); count | q ~ Binomial(m, q)
    q=rng.beta(k+1,n-k,size=sims.size); sims+=rng.binomial(m,q)
k_obs=out['all']['exceed']
out['chance']=dict(expected=float(E),p_ge_obs_poisson=float(poisson.sf(k_obs-1,E)),p_ge_obs_exact=float(np.mean(sims>=k_obs)),
                   sim_mean=float(sims.mean()),sim_var=float(sims.var()),
                   ul95_total_poisson=float(chi2.ppf(0.95,2*(k_obs+1))/2),
                   ul95_excess_frac=float((chi2.ppf(0.95,2*(k_obs+1))/2-E)/out['all']['n_pass']),
                   ul95_frac_no_subtraction=float(chi2.ppf(0.95,2*(k_obs+1))/2/out['all']['n_pass']),
                   one_sided_ul95_shadow_1_of_380=float(_b.ppf(0.95,2,379)))
json.dump(out,open('../results/uemr_search.json','w'),indent=1); print(json.dumps(out,indent=1))
