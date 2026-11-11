# Cross-receiver coincidence of the Phase 2 blind-search detections (Sec. 10.2). For every satellite pass recorded at two
# or more receivers in the same five-minute slot, count the receivers at which the pass was detected; do the same for the
# sign-reversed null searches of the same passes (null exceedances: z above the receiver threshold inside the search
# window). Under independence, the expected number of passes detected at two or more receivers follows from the
# per-receiver detection fractions; the time offset tau of a real echo must also agree between receivers, since it is a
# property of the satellite, not of the receiver. Reads the Phase 2 shards and results/phase2_summary.json; writes
# results/coincidence.json.
import os, glob, json, itertools, numpy as np
S=json.load(open('../results/phase2_summary.json'))
TH={st:v['z_th'] for st,v in S['thresholds'].items()}
seen={}
for f in sorted(glob.glob(os.environ.get('P2GLOB','../results/p2_[cd][le]*.json'))):
    for r in json.load(open(f)):
        if r['case'] in ('pass','null'): seen[(r['case'],r['file'],r['obj'])]=r
inside=lambda r: abs(r['tau'])<6 and abs(r['nu'])<12
DET={(d['file'],d['obj']):d for d in S['detections']}
CLS={}
for d in S['detections']: CLS[d['obj']]=d['class']
for (c,f,o),r in seen.items():
    if c=='pass': CLS.setdefault(o,r['cls'])
def group(case):
    g={}
    for (c,f,o),r in seen.items():
        if c!=case: continue
        g.setdefault((f[:13],o),[]).append(r)
    return g
out={}
for case in ('pass','null'):
    g=group(case); res={}
    for cls in ('Starlink V2-mini','Starlink DTC','Starlink V1.x','AST SpaceMobile','all non-ISS'):
        rows=[(k,v) for k,v in g.items() if len(v)>=2 and (CLS.get(k[1])==cls or (cls=='all non-ISS' and CLS.get(k[1])!='ISS complex'))]
        if not rows: continue
        def isdet(r):
            if case=='pass': return (r['file'],r['obj']) in DET
            return r['z']>=TH[r['file'][-6:]] and inside(r)
        n2=0; k1=0; kge2=0; exp2=0.0; dtau=[]
        for k,v in rows:
            d=[isdet(r) for r in v]; n2+=1; k1+=sum(d)>=1; kge2+=sum(d)>=2
            taus=[(DET[(r['file'],r['obj'])]['tau'] if case=='pass' else r['tau']) for r,x in zip(v,d) if x]
            dtau+=[abs(a-b) for a,b in itertools.combinations(taus,2)]
        # expected >=2 under independence, from the per-receiver detection fraction in these groups
        pr={}
        for k,v in rows:
            for r in v: pr.setdefault(r['file'][-6:],[]).append(isdet(r))
        p={st:np.mean(x) for st,x in pr.items()}
        for k,v in rows:
            ps=[p[r['file'][-6:]] for r in v]; p0=np.prod([1-q for q in ps]); p1=sum(q*np.prod([1-w for j,w in enumerate(ps) if j!=i]) for i,q in enumerate(ps))
            exp2+=1-p0-p1
        res[cls]=dict(n_groups=n2,detected_at_least_1=int(k1),detected_at_least_2=int(kge2),expected_at_least_2_if_independent=float(exp2),
                      per_receiver_fraction={st:float(v) for st,v in p.items()},
                      dtau_s_median=float(np.median(dtau)) if dtau else None,dtau_s_frac_le1=float(np.mean(np.array(dtau)<=1.0)) if dtau else None,n_pairs=len(dtau))
    out[case]=res
json.dump(out,open('../results/coincidence.json','w'),indent=1)
for case in out:
    for c,v in out[case].items(): print(case,c,{k:(round(x,3) if isinstance(x,float) else x) for k,x in v.items() if k!='per_receiver_fraction'})
