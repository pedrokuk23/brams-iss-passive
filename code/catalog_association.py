# Could a Phase 2 detection belong to an object that was not searched? (Sec. 10.2)
# Phase 2 searched only the target classes; debris, rocket bodies and other payloads were not searched, so the
# association rule could not attribute an echo to them. For every Phase 2 detection, the echo track of the attributed
# object (its predicted bistatic Doppler shifted by the fitted tau and nu) is compared over +-10 s around the echo peak
# with the predicted Doppler of every other object of the Phase 2 element sets (all objects with mean motion above
# 14.85 rev/day, debris included). Another object is compatible when, for some tau' in [-6, 6] s and nu' in [-12, 12] Hz
# (the search grid), its track matches the echo track to better than 3 Hz rms. Writes results/catalog_association.json.
import os, sys, glob, json, numpy as np
from sgp4.api import Satrec, SatrecArray
from sgp4.propagation import gstime
import fastsearch as F
S=json.load(open('../results/phase2_summary.json')); META=json.load(open('../results/phase2_meta.json'))
LAM=F.LAM
# element sets: per object, all sets (epoch, line1, line2, name)
SETS={}
files=sorted(glob.glob('../data/tle/phase2/*.txt'))
if not files: sys.exit('catalog_association: no element sets in data/tle/phase2')
for f in files:
    L=open(f).read().splitlines()
    for i in range(len(L)-2):
        if L[i].startswith('0 ') and L[i+1].startswith('1 ') and L[i+2].startswith('2 '):
            s=Satrec.twoline2rv(L[i+1],L[i+2]); SETS.setdefault(L[i+1][2:7],[]).append((s.jdsatepoch+s.jdsatepochF,L[i+1],L[i+2],L[i][2:].strip()))
for k in SETS: SETS[k].sort()
BYNAME={}
for k,v in SETS.items():
    for x in v: BYNAME.setdefault(x[3],k)
CACHE={}
def catalog(jd0):
    # latest set of every object with epoch before the start of the recording's day (cached per day)
    key=int(jd0-0.5)
    if key not in CACHE:
        recs=[];names=[]
        for k,v in SETS.items():
            c=[x for x in v if x[0]<=jd0+0.5]
            if not c: continue
            e=c[-1]
            if jd0-e[0]>7: continue
            recs.append(Satrec.twoline2rv(e[1],e[2])); names.append(e[3])
        CACHE[key]=(SatrecArray(recs),names,recs)
    return CACHE[key]
def ecef_doppler(SA,t):
    jd=t/86400.0+2440587.5; jdi=np.floor(jd-0.5)+0.5; fr=jd-jdi
    e,r,v=SA.sgp4(jdi,fr); g=np.array([gstime(a+b) for a,b in zip(jdi,fr)])
    X=np.stack([r[...,0]*np.cos(g)+r[...,1]*np.sin(g),-r[...,0]*np.sin(g)+r[...,1]*np.cos(g),r[...,2]],-1)
    return X,e
def bist(X,TX,RX,tt):
    R1=np.linalg.norm(X-TX,axis=-1); R2=np.linalg.norm(X-RX,axis=-1)
    return -np.gradient(R1+R2,tt,axis=-1)*1e3/LAM, np.sqrt(R1*R2)
def same_object(a,b):
    a=a.upper(); b=b.upper()
    if a==b: return True
    iss=('ISS','DRAGON','SOYUZ','PROGRESS','CYGNUS','HTV','DREAM CHASER','STARLINER','NAUKA','POISK','ZVEZDA','ZARYA','UNITY','DESTINY')
    return any(x in a for x in iss) and any(x in b for x in iss)
out=[]; dets=[d for d in S['detections'] if d['class']!='ISS complex']
if os.environ.get('NMAX'): dets=dets[:int(os.environ['NMAX'])]
for n,d in enumerate(dets):
    md=META[d['file']]; t0=md['start_us']/1e6
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    SA,names,recs=catalog(t0/86400+2440587.5)
    # echo track of the attributed object
    rec=SETS[BYNAME[d['obj']]]
    jd0=t0/86400+2440587.5; best=min(rec,key=lambda x:abs(x[0]-jd0)); sat=Satrec.twoline2rv(best[1],best[2])
    tg,fdg,_=F.doppler(sat,t0,TX,RX)
    tw=d['t_pk']+np.arange(-10,10.01,0.5); trk=np.interp(tw+d['tau'],tg,fdg)+d['nu']
    ok=(np.abs(trk)>30)&(np.abs(trk)<900)
    if ok.sum()<5: out.append(dict(file=d['file'],obj=d['obj'],cls=d['class'],n_compatible=None)); continue
    # prefilter: every object at the peak time, range and Doppler
    tt=t0+d['t_pk']+np.array([-0.5,0.0,0.5]); X,e=ecef_doppler(SA,tt); fd,rr=bist(X,TX,RX,tt)
    cand=np.flatnonzero((e[:,1]==0)&(rr[:,1]<1500)&(np.abs(fd[:,1]-trk[20])<400))
    comp=[]
    if len(cand):
        tt=t0+d['t_pk']+np.arange(-16,16.01,0.5)
        Xa,ea=ecef_doppler(SatrecArray([recs[c] for c in cand]),tt); fda,_=bist(Xa,TX,RX,tt)
        for j,c in enumerate(cand):
            if same_object(names[c],d['obj']): continue
            bestr=None; adm=[]
            for tp in F.TAUS:
                f2=np.interp(tw+tp,tt-t0,fda[j]); dif=(f2-trk)[ok]; m=np.mean(dif)
                if abs(m)>12: continue
                r=np.sqrt(np.mean((dif-m)**2))
                if r<3.0: adm.append(float(tp))        # every tau of the grid at which this object fits (nu profiled)
                if bestr is None or r<bestr[0]: bestr=(r,tp,m)
            if bestr is not None and bestr[0]<3.0: comp.append(dict(name=names[c],rms_hz=float(bestr[0]),tau=float(bestr[1]),nu=float(bestr[2]),taus_ok=adm))
    out.append(dict(file=d['file'],obj=d['obj'],cls=d['class'],z=d['z'],tau=d['tau'],n_compatible=len(comp),compatible=comp))
    if n%100==0: print(n,len(dets),flush=True)
res=[o for o in out if o['n_compatible'] is not None]
# multi-receiver test: a real echo of the attributed object has the same tau at every receiver of the slot (tau is a
# property of the orbit); an alternative object must also fit every receiver with one common tau' (within 1 s), tested on
# the full sets of admissible tau (rms < 3 Hz with nu profiled), not on the best-fitting tau alone
G={}
for o in res: G.setdefault((o['file'][:13],o['obj']),[]).append(o)
for k,g in G.items():
    if len(g)<2: continue
    alts=set.intersection(*[{a['name'] for a in o['compatible']} for o in g])
    surv=[]
    for nm in alts:
        # admissible sets of tau at each receiver; the alternative survives if one 1-s window holds an admissible tau of
        # every receiver (a common timing offset), not only if the best-fitting taus agree
        S_=[[a['taus_ok'] for a in o['compatible'] if a['name']==nm][0] for o in g]
        if any(all(any(t0<=t<=t0+1.0 for t in s) for s in S_) for s0 in S_ for t0 in s0): surv.append(nm)
    for o in g: o['multi_rx']=len(g); o['alternatives_surviving_multi_rx']=surv
multi=[o for o in res if o.get('multi_rx',0)>=2]
summ=dict(n_detections=len(out),n_checked=len(res),n_with_alternative=sum(o['n_compatible']>0 for o in res),
          by_class={c:dict(n=sum(o['cls']==c for o in res),with_alternative=sum(o['cls']==c and o['n_compatible']>0 for o in res)) for c in sorted({o['cls'] for o in res})},
          n_multi_rx=len(multi),n_multi_rx_with_alternative=sum(o['n_compatible']>0 for o in multi),
          n_multi_rx_alternative_surviving=sum(len(o['alternatives_surviving_multi_rx'])>0 for o in multi),
          ambiguous_by_class={c:sum(o['cls']==c and o['n_compatible']>0 and (o.get('multi_rx',0)<2 or len(o['alternatives_surviving_multi_rx'])>0) for o in res) for c in sorted({o['cls'] for o in res})},
          n_objects_in_catalog=len(SETS),
          alt_debris=sum(any(' DEB' in a['name'] or 'R/B' in a['name'] for a in o['compatible']) for o in res if o['n_compatible']))
json.dump(dict(summary=summ,rows=out),open('../results/catalog_association.json','w'),indent=1); print(json.dumps(summ,indent=1))
