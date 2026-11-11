import json,numpy as np,glob,os,datetime as dt
from scipy import stats
def ecef(lat,lon,h):
    a=6378.137; f=1/298.257223563; e2=f*(2-f); la,lo=np.radians(lat),np.radians(lon); N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.stack([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)],-1)
rows=[]
for day,fn in [('20261103','../data/adsb/adsb_be_nov03.jsonl'),('20261104','../data/adsb/adsb_be_nov04.jsonl')]:
    AC={}
    for l in open(fn):
        d=json.loads(l); P=[p for p in d['pts'] if p[3] not in (None,'ground')]
        if len(P)<3: continue
        ts=np.array([p[0] for p in P],float); alt=np.array([(p[6] if isinstance(p[6],(int,float)) else p[3]) for p in P],float)*0.3048/1000
        AC[d.get('r') or d.get('icao')]=(ts,ecef(np.array([p[1] for p in P]),np.array([p[2] for p in P]),alt))
    for f in glob.glob(f'../results/uni_*_{day}.json'):
        for r in json.load(open(f)): 
            if r['z']>=6 and abs(r['tau_ref'])<=1 and r['ac'] in AC: rows.append((r,AC[r['ac']]))
# station medians
med={}
for r,_ in rows: med.setdefault((r['st'],r['day']),[]).append(r['tau_ref'])
med={k:np.median(v) for k,v in med.items()}
out=[]
for r,(ts,X) in rows:
    t0=dt.datetime.strptime(r['day']+r['file'],'%Y%m%d%H%M').replace(tzinfo=dt.timezone.utc).timestamp()
    w=(ts>=t0-30)&(ts<=t0+330); idx=np.nonzero(w)[0]
    if len(idx)<4: continue
    gap=np.median(np.diff(ts[idx])); 
    errs=[]
    for i in idx:
        if i==0 or i==len(ts)-1: continue
        a=(ts[i]-ts[i-1])/(ts[i+1]-ts[i-1]); Xi=X[i-1]+a*(X[i+1]-X[i-1])
        v=(X[i+1]-X[i-1])/(ts[i+1]-ts[i-1]); sp=np.linalg.norm(v)
        if sp<0.05: continue
        errs.append(np.dot(Xi-X[i],v/sp)/sp)
    errs=np.array(errs)
    out.append(dict(st=r['st'],day=r['day'],ac=r['ac'],dev=r['tau_ref']-med[(r['st'],r['day'])],gap=gap,loo=np.sqrt(np.mean(errs**2)) if len(errs) else np.nan,z=r['z']))
dev=np.array([o['dev'] for o in out]); gap=np.array([o['gap'] for o in out]); loo=np.array([o['loo'] for o in out])
print('n',len(out),'gap median',np.median(gap),'IQR',np.percentile(gap,[25,75]),'max',gap.max())
print('LOO along-track timing rms (s): median',np.nanmedian(loo),'IQR',np.nanpercentile(loo,[25,75]))
print('|dev| median',np.median(abs(dev)),'robust std',1.4826*np.median(abs(dev-np.median(dev))))
print('spearman |dev| vs gap',stats.spearmanr(abs(dev),gap)); 
m=np.isfinite(loo); print('spearman |dev| vs loo',stats.spearmanr(abs(dev[m]),loo[m]))
for th in [8,12,20]:
    s=gap<=th; print(f'gap<={th}s n={s.sum()} robust std={1.4826*np.median(abs(dev[s]-np.median(dev[s]))):.3f}  | gap>{th} n={(~s).sum()} robust std={1.4826*np.median(abs(dev[~s]-np.median(dev[~s]))) if (~s).sum()>3 else np.nan:.3f}')
json.dump(out,open('../results/adsb_latency.json','w'),indent=0)
