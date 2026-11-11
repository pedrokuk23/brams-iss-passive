# Phase 2 analysis (PROTOCOL.md item 6), identical to blind_analyze.py (items 3-5): thresholds from the sign-reversed nulls, association of candidates to
# echoes (also applied to the nulls), iteration to a fixed point; detection probability per object class and receiver with
# Clopper-Pearson intervals; expected false alarms; detection probability versus the cross section needed by the model.
import os, json, glob, re, numpy as np
from scipy.stats import beta
# Phase 2 (PROTOCOL item 6): the same rule, applied unchanged to the validation campaign, 19 Oct - 2 Nov 2026.
_seen={}
for f in sorted(glob.glob(os.environ.get('P2GLOB','../results/p2_[cd][le]*.json'))):   # shards of one run only (run_all.sh sets P2GLOB)
    for r in json.load(open(f)):
        if r['case'] in ('pass','null'): _seen[(r['case'],r['file'],r['obj'])]=r
R=list(_seen.values())
inside=lambda r: abs(r['tau'])<6 and abs(r['nu'])<12
STS=['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']; FAR=1e-3
def cp(k,n,a=0.05): return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)), 1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
P=[r for r in R if r['case']=='pass']; N=[r for r in R if r['case']=='null']
for r in P: r['class']=r['cls']          # sigma_min from phase2_plan.py (same model)
# PROTOCOL item 6b (correction of the class assignment of phase2_plan.py, 9 Nov 2026): objects named SPACEMOBILE-SSRnn DEB
# are launch debris, not BlueBirds, and are reported as 'Other'; BlueBird Block 1 (SPACEMOBILE-001 to -005) carries a
# 64-m2 array like BlueWalker 3, not 223 m2, so its passes with sigma_min above 4 pi A^2/lambda^2 for 64 m2 are outside the
# physical limit and are reported separately (they were searched under the 223-m2 limit).
LAM=299792458/49.97e6
for r in P:
    if r['class']=='AST SpaceMobile' and 'DEB' in r['obj']: r['class']='Other'
    elif r['class']=='AST SpaceMobile' and (r['obj'].startswith('BLUEWALKER') or r['obj'][-3:] in ('001','002','003','004','005')):
        if r['sigma_min']>4*np.pi*64**2/LAM**2: r['class']='AST beyond physical limit'
# echo track of an owner: raw frequency f(t) = f_D(t + tau) + nu predicted for that object in that recording
import os, glob as _g, bwav, fastsearch as F
from sgp4.api import Satrec
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '): sets.setdefault(L[k][2:].strip(),[]).append((L[k+1],L[k+2]))
D=os.environ.get('P2DIR','../data/raw_wav_phase2/')
PATH={os.path.basename(p)[11:31]:p for d in D.split(':') if d for p in _g.glob(os.path.join(d,'*.wav'))}
_md={}
META=json.load(open('../results/phase2_meta.json')) if os.path.exists('../results/phase2_meta.json') else {}   # WAV header fields of the 788 recordings; read from the WAV files when absent
def meta(f):
    if f not in _md:
        md=META[f] if f in META else bwav.read(PATH[f])[0]; _md[f]=(md['start_us']/1e6,F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000),F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000))
    return _md[f]
EP={}
def satrec(nm,jd0):
    if nm not in EP: EP[nm]=np.array([Satrec.twoline2rv(*q).jdsatepoch+Satrec.twoline2rv(*q).jdsatepochF for q in sets[nm]])
    return Satrec.twoline2rv(*sets[nm][int(np.argmin(abs(EP[nm]-jd0)))])
_trk={}
def track(o):
    key=(o['file'],o['obj'])
    if key not in _trk:
        t0,TX,RX=meta(o['file']); tg,fdg,rr=F.doppler(satrec(o['obj'],t0/86400+2440587.5),t0,TX,RX); _trk[key]=(tg,fdg)
    return _trk[key]
def near(o,r):          # does the peak of r fall on the echo track of owner o? (3 Hz, allowing 1.5 s for the block timing,
    tg,fdg=track(o)     # since near closest approach the track can sweep tens of hertz per second)
    tt=r['t_pk']+np.linspace(-1.5,1.5,61); return np.min(np.abs(r['f_raw']-(np.interp(tt+o['tau'],tg,fdg)+o['nu'])))<3.0
def thresholds(explained):
    th={}
    for st in STS:
        nn=[r for r in N if r['file'].endswith(st)]; nz=np.sort([r['z'] for r in nn if inside(r) and not explained.get(id(r))])[::-1]
        kmax=max(int(np.floor(FAR*len(nn))),0)          # allowed exceedances
        th[st]=float(nz[kmax]+1e-9) if len(nz)>kmax else 0.0
    return th
def associate(th):
    for r in P: r['det']=False; r.pop('ambiguous_with',None)
    owners={}
    for r in sorted([r for r in P if inside(r) and r['z']>=th[r['file'][-6:]]],key=lambda r:-r['z']):
        o=next((q for q in owners.setdefault(r['file'],[]) if near(q,r)),None)
        if o is None: owners[r['file']].append(r); r['det']=True
        else: r['ambiguous_with']=o['obj']
    explained={id(n):True for n in N if any(near(o,n) for o in owners.get(n['file'],[]))}
    return owners,explained
th={st:7.0 for st in STS}; hist=[th]          # start from below (z = 7) so that strong echoes can explain the nulls that cross them
for it in range(10):
    owners,expl=associate(th); th2=thresholds(expl); hist.append(th2)
    if all(abs(th2[s]-th[s])<1e-6 for s in STS): break
    th=th2
owners,expl=associate(th)
out=dict(iterations=len(hist),thresholds={},classes={},detections=[],ambiguous=[])
for st in STS:
    nn=[r for r in N if r['file'].endswith(st)]; ne=sum(1 for r in nn if inside(r) and not expl.get(id(r)) and r['z']>=th[st])
    nz7=sum(1 for r in nn if inside(r) and not expl.get(id(r)) and r['z']>=7)
    out['thresholds'][st]=dict(n_null=len(nn),z_th=th[st],exceed=ne,far=ne/len(nn),far_ci=cp(ne,len(nn)),far_at_7=nz7/len(nn),
        n_null_explained=sum(1 for r in nn if expl.get(id(r))))
# per class: an object pass in a recording counts once; for the ISS complex, one pass per recording (any member detected)
def units(c):
    pp=[r for r in P if r['class']==c]
    if c!='ISS complex': return [[r] for r in pp]
    g={}
    for r in pp: g.setdefault(r['file'],[]).append(r)
    return list(g.values())
for c in sorted(set(r['class'] for r in P)):
    U=units(c)
    for st in STS+['all']:
        uu=[u for u in U if st=='all' or u[0]['file'].endswith(st)]
        if not uu: continue
        k=sum(any(r['det'] for r in u) for u in uu); n=len(uu)
        out['classes'].setdefault(c,{})[st]=dict(n=n,k=k,pd=k/n,ci=cp(k,n),exp_false=n*FAR)
for r in P:
    if r['det']: out['detections'].append({k:r[k] for k in ('file','obj','class','z','tau','nu','rrmin','sigma_min','t_pk','f_raw')})
    if r.get('ambiguous_with'): out['ambiguous'].append(dict(file=r['file'],obj=r['obj'],z=r['z'],attributed_to=r['ambiguous_with']))
bins=[0,30,100,300,1e3,3e3,1e4,1e9]; out['pd_vs_sigma_min']={}
for c in sorted(set(r['class'] for r in P)):
    rows=[]
    for a,b in zip(bins[:-1],bins[1:]):
        pp=[r for r in P if r['class']==c and a<=r['sigma_min']<b]
        if pp: k=sum(r['det'] for r in pp); rows.append(dict(lo=a,hi=b,n=len(pp),k=k,ci=cp(k,len(pp))))
    out['pd_vs_sigma_min'][c]=rows
out['n_pass']=len(P); out['n_null']=len(N)
# checks of the false-alarm rule (far_checks.py): out-of-sample rate by split halves, effective rate, and the null set
from far_checks import far_checks
out['far_checks']=far_checks(N,expl,th,STS,inside)
out['nulls']=[[r['file'][-6:],r['file'],round(r['z'],4),bool(inside(r)),bool(expl.get(id(r)))] for r in N]
# sensitivity: thresholds from all nulls, without removing the nulls that fall on real echoes
dets0=[r['det'] for r in P]; th_all=thresholds({}); associate(th_all); out['no_null_removal']={'thresholds':th_all,'classes':{}}
for c in sorted(set(r['class'] for r in P)):
    U=units(c); k=sum(any(r['det'] for r in u) for u in U); out['no_null_removal']['classes'][c]=dict(n=len(U),k=k,ci=cp(k,len(U)))
for r,d in zip(P,dets0): r['det']=d
# sensitivity: thresholds frozen at their Phase 1 values (results/blind_summary.json) instead of re-estimated
th1={k:v['z_th'] for k,v in json.load(open('../results/blind_summary.json'))['thresholds'].items()}
dets=[r['det'] for r in P]; owners1,expl1=associate(th1); out['frozen_phase1_thresholds']={'thresholds':th1,'classes':{},'far':{}}
for st in STS:
    nn=[r for r in N if r['file'].endswith(st)]; ne=sum(1 for r in nn if inside(r) and not expl1.get(id(r)) and r['z']>=th1[st])
    out['frozen_phase1_thresholds']['far'][st]=dict(n_null=len(nn),exceed=ne,far_ci=cp(ne,len(nn)))
for c in sorted(set(r['class'] for r in P)):
    U=units(c); k=sum(any(r['det'] for r in u) for u in U); out['frozen_phase1_thresholds']['classes'][c]=dict(n=len(U),k=k,ci=cp(k,len(U)))
out['frozen_phase1_thresholds']['detections']=[dict(file=r['file'],obj=r['obj']) for r in P if r['det']]
for r,d in zip(P,dets): r['det']=d
json.dump(out,open('../results/phase2_summary.json','w'),indent=1)
print('passes',len(P),'nulls',len(N),'iterations',len(hist))
for st,v in out['thresholds'].items(): print(st,'n_null %d z_th %.2f far %.4f [%.4f,%.4f] far@7 %.4f explained %d'%(v['n_null'],v['z_th'],v['far'],*v['far_ci'],v['far_at_7'],v['n_null_explained']))
for c,d in out['classes'].items():
    a=d['all']; print('%-17s n=%5d k=%4d Pd=%.3f [%.3f,%.3f] exp.false=%.1f'%(c,a['n'],a['k'],a['pd'],a['ci'][0],a['ci'][1],a['exp_false']),{st:f"{d[st]['k']}/{d[st]['n']}" for st in STS if st in d})
print('ambiguous',len(out['ambiguous']))
