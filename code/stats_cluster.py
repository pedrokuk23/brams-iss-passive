# Uncertainties that respect the grouping of the data (Table 7 uses the two-way interval; the ISS complex the interval
# over physical passes). Passes searched in the same recording share its noise and
# interference (and Starlink satellites of one train pass together), so detection probabilities are given with 95%
# intervals from a bootstrap over recordings (files), next to the binomial Clopper-Pearson interval. The feasibility
# prediction of Phase 2 from Phase 1 gets a predictive interval that includes the uncertainty of the Phase 1 curve
# (bootstrap over Phase 1 recordings) and the Bernoulli outcome of each Phase 2 pass, and its error is decomposed by
# attitude group, class and receiver.
import json, collections, numpy as np
exec(open('feasibility_validate.py').read().split("# 1. receiver gain corrections")[0])   # U1, U2 (Phase 1 thresholds), fit_curve, pd_of
C0=json.load(open('../results/feasibility_curve.json')); DG.update(C0['gain_correction_dB'])
D2main={(d['file'],d['obj']) for d in json.load(open('../results/phase2_summary.json'))['detections']}
U2m=units(P2,D2main,sat2)                     # Phase 2 with its own thresholds (primary result)
for U in (U1,U2,U2m): setM(U)
rng=np.random.default_rng(2026); NB=1000
def files(U):
    g=collections.defaultdict(list)
    for u in U: g[u['file']].append(u)
    return list(g.values())
def boot_rate(U):
    F=files(U); k=np.array([sum(u['det'] for u in f) for f in F]); n=np.array([len(f) for f in F]); r=[]
    for _ in range(NB):
        i=rng.integers(0,len(F),len(F)); r.append(k[i].sum()/n[i].sum())
    return [float(x) for x in np.percentile(r,[2.5,97.5])]
def boot_obj(U):
    """Bootstrap over objects (all passes of one satellite together), for the correlation of attitude and orbit."""
    g=collections.defaultdict(list)
    for u in U: g[u['r']['obj'] if 'r' in u else u['file']].append(u)
    G=list(g.values()); k=np.array([sum(u['det'] for u in x) for x in G]); n=np.array([len(x) for x in G]); r=[]
    for _ in range(NB):
        i=rng.integers(0,len(G),len(G)); r.append(k[i].sum()/n[i].sum())
    return [float(x) for x in np.percentile(r,[2.5,97.5])]
def tclusters(U):
    """Time clusters: five-minute slots merged when consecutive (an ISS pass can span two slots), shared by receivers."""
    sl=sorted({u['file'][:13] for u in U}); import datetime as _d
    ts=[_d.datetime.strptime(s,'%Y%m%d_%H%M') for s in sl]; lab={}; c=0
    for i,s in enumerate(sl):
        if i and (ts[i]-ts[i-1]).total_seconds()>300: c+=1
        lab[s]=c
    return {u['file'][:13]:lab[u['file'][:13]] for u in U}
def boot_twoway(U):
    """Two-way (pigeonhole) bootstrap: satellites and time clusters resampled independently, each pass weighted by the
    product of the multiplicities of its satellite and its time cluster (Owen 2007)."""
    tc=tclusters(U); so=sorted({(u['r']['obj'] if 'r' in u else 'x') if u['cls']!='ISS complex' else 'ISS' for u in U}); si={s:i for i,s in enumerate(so)}
    a=np.array([si[(u['r']['obj'] if 'r' in u else 'x') if u['cls']!='ISS complex' else 'ISS'] for u in U]); b=np.array([tc[u['file'][:13]] for u in U])
    d=np.array([u['det'] for u in U],float); nA=len(so); nB=b.max()+1; r=[]
    for _ in range(NB):
        wa=np.bincount(rng.integers(0,nA,nA),minlength=nA); wb=np.bincount(rng.integers(0,nB,nB),minlength=nB); w=wa[a]*wb[b]
        if w.sum()>0: r.append((w*d).sum()/w.sum())
    return [float(x) for x in np.percentile(r,[2.5,97.5])], int(nA), int(nB)
def boot_pass(U):
    """ISS complex: bootstrap over physical passes (time clusters), all receivers and recordings of a pass together."""
    tc=tclusters(U); b=np.array([tc[u['file'][:13]] for u in U]); d=np.array([u['det'] for u in U],float); nB=b.max()+1; r=[]
    k=np.bincount(b,weights=d,minlength=nB); n=np.bincount(b,minlength=nB)
    for _ in range(NB):
        i=rng.integers(0,nB,nB); r.append(k[i].sum()/n[i].sum())
    return [float(x) for x in np.percentile(r,[2.5,97.5])], int(nB)
def row(U):
    k=sum(u['det'] for u in U); n=len(U)
    out=dict(n=n,k=int(k),pd=k/n if n else None,ci_binomial=cp(k,n) if n else None,ci_cluster=boot_rate(U) if n else None,
                ci_object=boot_obj(U) if n else None,n_files=len(files(U)))
    if n:
        if all(u['cls']=='ISS complex' for u in U): out['ci_pass'],out['n_passes']=boot_pass(U)
        out['ci_twoway'],out['n_sat'],out['n_timeclusters']=boot_twoway(U)
    return out
out=dict(n_boot=NB,classes={},eclipse={},feasibility={},classes_phase2_without_maubeuge={})
for ph,U in (('phase1',U1),('phase2',U2m)):
    out['classes'][ph]={c:row([u for u in U if u['cls']==c]) for c in CL}
    out['eclipse'][ph]={f'{c}|{g}':row([u for u in U if u['cls']==c and u['grp']==g]) for c in CL if c!='ISS complex' for g in ('plate_shadow','plate_sunlit')}
out['classes_phase2_without_maubeuge']={c:row([u for u in U2m if u['cls']==c and not u['file'].endswith('FRMAUB')]) for c in CL}
# feasibility: curve from Phase 1, prediction of Phase 2 (Phase 1 thresholds in both phases)
C=fit_curve(U1); F1=files(U1); F2=files(U2); P=[]; O=[]
pred=sum(pd_of(C,u) for u in U2); obs=sum(u['det'] for u in U2)
for _ in range(NB):
    Cb=fit_curve([u for i in rng.integers(0,len(F1),len(F1)) for u in F1[i]])
    p=np.array([pd_of(Cb,u) for u in U2]); P.append(int((rng.random(len(p))<p).sum()))
    O.append(sum(u['det'] for i in rng.integers(0,len(F2),len(F2)) for u in F2[i]))
Dd=np.array(O)-np.array(P)       # observed minus predicted, both resampled: Phase 2 recordings and Phase 1 curve + outcomes
fe=dict(observed=int(obs),predicted=round(float(pred),1),pred_95=[float(x) for x in np.percentile(P,[2.5,97.5])],obs_95=[float(x) for x in np.percentile(O,[2.5,97.5])],
        ratio=round(obs/pred,3),p_obs_ge=float(np.mean(np.array(P)>=obs)),diff_95=[float(x) for x in np.percentile(Dd,[2.5,97.5])],
        p_diff_two_sided=float(2*min(np.mean(Dd<=0),np.mean(Dd>=0))))
def dec(sel): return dict(n=len(sel),observed=int(sum(u['det'] for u in sel)),predicted=round(float(sum(pd_of(C,u) for u in sel)),1))
fe['by_group']={g:dec([u for u in U2 if u['grp']==g]) for g in ('structure','plate_shadow','plate_sunlit')}
fe['by_class']={c:dec([u for u in U2 if u['cls']==c]) for c in CL}
fe['by_receiver']={st:dec([u for u in U2 if u['file'].endswith(st)]) for st in STS}
fe['by_margin']={f'{g}|{a}-{b}':dec([u for u in U2 if u['grp']==g and a<=u['M']<b]) for g in ('structure','plate_shadow','plate_sunlit') for a,b in zip(EDGES[:-1],EDGES[1:])}
fe['M_ge_10']={g:dec([u for u in U2 if u['grp']==g and u['M']>=10]) for g in ('structure','plate_shadow','plate_sunlit')}
out['feasibility']=fe
json.dump(out,open('../results/cluster_stats.json','w'),indent=1)
for ph in ('phase1','phase2'):
    for c,v in out['classes'][ph].items(): print(ph,'%-17s %4d/%5d = %.3f  binomial [%.3f,%.3f] cluster [%.3f,%.3f]  files %d'%(c,v['k'],v['n'],v['pd'],*v['ci_binomial'],*v['ci_cluster'],v['n_files']))
    for c,v in out['eclipse'][ph].items(): print(ph,'%-30s %4d/%5d = %.3f cluster [%.3f,%.3f]'%(c,v['k'],v['n'],v['pd'],*v['ci_cluster']))
print('feasibility',{k:fe[k] for k in ('observed','predicted','pred_95','obs_95','ratio','p_obs_ge')})
for k in ('by_group','by_class','by_receiver','M_ge_10'): print(k,fe[k])
