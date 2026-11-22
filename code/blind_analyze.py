# Blind search, analysis (PROTOCOL.md items 3-5): thresholds from the sign-reversed nulls, association of candidates to
# echoes (also applied to the nulls), iteration to a fixed point; detection probability per object class and receiver with
# Clopper-Pearson intervals; expected false alarms; detection probability versus the cross section needed by the model.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os, json, glob, re, numpy as np
from scipy.stats import beta
# results of all runs (cloud shards, runs on the author's computer, Envisat recordings); one record per (case, file, object)
_seen={}; 
_G=os.environ.get('BLINDGLOB')                     # one run only (run_all.sh sets it); default: the published runs
for f in sorted(glob.glob(_G) if _G else glob.glob('../results/blind_shard*.json')+glob.glob('../results/blind_dev*.json')+glob.glob('../results/blind_ev*.json')):
    for r in json.load(open(f)):
        if r['case'] in ('pass','null'): _seen[(r['case'],r['file'],r['obj'])]=r
R=list(_seen.values())
inside=lambda r: abs(r['tau'])<6 and abs(r['nu'])<12
STS=['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']; FAR=1e-3
def cp(k,n,a=0.05): return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)), 1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
ISSC=('ISS','ISS (','POISK','PROGRESS-MS','SOYUZ-MS','CREW DRAGON','DRAGON CRS','CYGNUS','HTV-X')
def cls(n):
    if n=='ISS' or any(n.startswith(p) for p in ISSC[1:]): return 'ISS complex'
    if n.startswith('SPACEMOBILE') or n.startswith('BLUEWALKER'): return 'AST SpaceMobile'
    if '[DTC]' in n: return 'Starlink DTC'
    m=re.match(r'STARLINK-(\d+)',n)
    if m: return 'Starlink V2-mini' if int(m.group(1))>=30000 else 'Starlink V1.x'
    return 'Other'
T=json.load(open('../data/planning/planning_model.json')); G={k:10**(v/10) for k,v in T['Geff_dB'].items()}
def sigma_min(r):
    g=G.get(r['file'][-6:],G['BEOTTI'])
    cn0=T['Pt_W']*g*T['lam_m']**2/((4*np.pi)**3*(r['rrmin']*1e3)**4*1.380649e-23*T['Tsky_K'])
    return 70/min(r['n_valid_s'],50)/cn0
P=[r for r in R if r['case']=='pass']; N=[r for r in R if r['case']=='null']
for r in P: r['class']=cls(r['obj']); r['sigma_min']=sigma_min(r)
# echo track of an owner: raw frequency f(t) = f_D(t + tau) + nu predicted for that object in that recording
import os, glob as _g, bwav, fastsearch as F
from sgp4.api import Satrec
exec(open('blind_enumerate.py').read().split("names=[x[0]")[0])
TLEc={nm:(l1,l2) for nm,l1,l2 in cat}
D=os.environ.get('WAVDIR','../data/raw_wav/'); EV=os.environ.get('EVDIR','../data/raw_wav_envisat/')
PATH={os.path.basename(p)[11:31]:p for p in _g.glob(D+'*.wav')+_g.glob(EV+'*.wav')}
_md={}
def meta(f):
    if f not in _md:
        md,_,_=bwav.read(PATH[f]); _md[f]=(md['start_us']/1e6,F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000),F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000))
    return _md[f]
def satrec(obj,jd0):
    if obj in special:
        sets=special[obj]; ep=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets])
        return Satrec.twoline2rv(*sets[int(np.argmin(abs(ep-jd0)))])
    return Satrec.twoline2rv(*TLEc[obj])
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
# checks of the false-alarm rule (far_checks.py): out-of-sample rate by split halves, effective rate, and the null set
from far_checks import far_checks
out['far_checks']=far_checks(N,expl,th,STS,inside)
out['nulls']=[[r['file'][-6:],r['file'],round(r['z'],4),bool(inside(r)),bool(expl.get(id(r)))] for r in N]
# sensitivity: thresholds from all nulls, without removing the nulls that fall on real echoes
dets0=[r['det'] for r in P]; th_all=thresholds({}); associate(th_all); out['no_null_removal']={'thresholds':th_all,'classes':{}}
for c in sorted(set(r['class'] for r in P)):
    U=units(c); k=sum(any(r['det'] for r in u) for u in U); out['no_null_removal']['classes'][c]=dict(n=len(U),k=k,ci=cp(k,len(U)))
for r,d in zip(P,dets0): r['det']=d
out['n_pass']=len(P); out['n_null']=len(N)
json.dump(out,open('../results/blind_summary.json','w'),indent=1)
print('passes',len(P),'nulls',len(N),'iterations',len(hist))
for st,v in out['thresholds'].items(): print(st,'n_null %d z_th %.2f far %.4f [%.4f,%.4f] far@7 %.4f explained %d'%(v['n_null'],v['z_th'],v['far'],*v['far_ci'],v['far_at_7'],v['n_null_explained']))
for c,d in out['classes'].items():
    a=d['all']; print('%-17s n=%5d k=%4d Pd=%.3f [%.3f,%.3f] exp.false=%.1f'%(c,a['n'],a['k'],a['pd'],a['ci'][0],a['ci'][1],a['exp_false']),{st:f"{d[st]['k']}/{d[st]['n']}" for st in STS if st in d})
print('ambiguous',len(out['ambiguous']))
