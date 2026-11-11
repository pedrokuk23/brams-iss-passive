import numpy as np, json, os, glob, sys, datetime as dt
from sgp4.propagation import gstime
sys.argv=['x','X']
exec(open('station_run.py').read().split('print(f"== {ST} ==")')[0])
C=json.load(open('../data/stations/coords.json')); TXg=C['tx']
def enu(lat,lon,X0,X):
    la,lo=np.radians(lat),np.radians(lon); d=X-X0
    E=np.array([-np.sin(lo),np.cos(lo),0]); N=np.array([-np.sin(la)*np.cos(lo),-np.sin(la)*np.sin(lo),np.cos(la)]); U=np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
    e,n,u=d@E,d@N,d@U; return np.degrees(np.arctan2(e,n))%360, np.degrees(np.arctan2(u,np.hypot(e,n)))
TX=ecef(TXg[0],TXg[1],0.236)
# ---- aircraft: elevation at closest approach vs z
out={}
for day,fn in [('20261103','../data/adsb/adsb_be_nov03.jsonl'),('20261104','../data/adsb/adsb_be_nov04.jsonl')]:
    os.environ['ADSB']=fn; exec(open('air_overlay.py').read().split("files=")[0])
    tr={a[0]:a for a in AC}
    for f in glob.glob(f'../results/uni_*_{day}.json'):
        for r in json.load(open(f)):
            st=r['st']; la0,lo0=C['stations'][st]; RX=ecef(la0,lo0,0.1)
            if r['ac'] not in tr: continue
            name,typ,ts,la,lo,al=tr[r['ac']]
            t0=dt.datetime.strptime(day+r['file'],'%Y%m%d%H%M').replace(tzinfo=dt.timezone.utc).timestamp()
            tt=np.arange(t0,t0+300,2.0); s=(tt>=ts.min())&(tt<=ts.max())
            if s.sum()<5: continue
            X=ecef(np.interp(tt[s],ts,la),np.interp(tt[s],ts,lo),np.interp(tt[s],ts,al))
            g=np.linalg.norm(X-TX,axis=1)*np.linalg.norm(X-RX,axis=1); i=np.argmin(g)
            az,el=enu(la0,lo0,RX,X[i]); azt,elt=enu(TXg[0],TXg[1],TX,X[i])
            out.setdefault(st,[]).append((el,az,r["z"],elt,r["rr_km"]))
for st,v in out.items():
    v=np.array(v); 
    line=f"{st}: n={len(v)}"
    for lo_,hi_ in [(0,10),(10,20),(20,90)]:
        m=(v[:,0]>=lo_)&(v[:,0]<hi_)
        if m.sum(): line+=f" | el {lo_}-{hi_}: n={m.sum()} det={np.mean(v[m,2]>=6)*100:.0f}% medz={np.median(v[m,2]):.1f}"
    print(line)
json.dump({k:np.array(v).tolist() for k,v in out.items()},open('../results/aircraft_elev.json','w'))
# ---- ISS geometry in band
passes={'20261103_0040':None,'20261103_0220':None,'20261103_2355':None,'20261104_0130':None,'20261104_2305':None}
for p in passes:
    t0=dt.datetime.strptime(p,'%Y%m%d_%H%M').replace(tzinfo=dt.timezone.utc).timestamp()
    t=np.arange(t0-300,t0+600,1.0); X=iss_ecef(t)
    for st,(la0,lo0) in C['stations'].items():
        RX=ecef(la0,lo0,0.1); Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam
        m=np.abs(fd)<2700
        if not m.any(): print(p,st,'no in-band'); continue
        az,el=enu(la0,lo0,RX,X[m].T.T) if False else zip(*[enu(la0,lo0,RX,x) for x in X[m]])
        az=np.array(az); el=np.array(el); azt,elt=zip(*[enu(TXg[0],TXg[1],TX,x) for x in X[m]])
        print(f"{p} {st}: in-band {m.sum()} s, ISS el from RX {el.min():.0f}-{el.max():.0f}°, az {az.min():.0f}-{az.max():.0f}°, el from TX {min(elt):.0f}-{max(elt):.0f}°")
