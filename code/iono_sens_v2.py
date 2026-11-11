import numpy as np, json, sys
from sgp4.propagation import gstime
sys.argv=['x','X']
exec(open('station_run.py').read().split('print(f"== {ST} ==")')[0])
import bwav
import os
D=os.environ.get('WAVDIR','../data/raw_wav/'); tr=json.load(open('../results/iono_tracks_v2.json'))
PPS={r['file'][11:31]:r['off0_ms'] for r in json.load(open('../results/pps_audit.json'))}
c=299792458.0; f0=49.97e6; K=40.3/(c*f0); RE=6371.0
O=json.load(open('../results/oem_cmp_v2.json'))
def chap(h,hm,H): z=(h-hm)/H; return np.exp(0.5*(1-z-np.exp(-z)))
def stec(A,B,hm,H):
    s=np.linspace(0,1,400)[:,None]; P=A[None,:]+s*(B-A)[None,:]; h=np.linalg.norm(P,axis=1)-RE
    return np.trapezoid(chap(h,hm,H),dx=1/399)*np.linalg.norm(B-A)*1e3
def cleaned(a): a=np.array(a); med=np.median(a[:,1]); return a[abs(a[:,1]-med)<4]
PK={'0040 03/11':'20261103_0040','0220 03/11':'20261103_0220','2355 03/11':'20261103_2355','0130 04/11':'20261104_0130','2305 04/11':'20261104_2305'}
passes={p:[k for k in tr if k.startswith(v) and len(cleaned(tr[k]))>=19] for p,v in PK.items()}
iono={'2355 03/11':(5.9,623),'0130 04/11':(5.9,649),'2305 04/11':(5.35,418)}
print({p:[k[14:] for k in v] for p,v in passes.items()})
geo={}
for p,keys in passes.items():
    for key in keys:
        if key not in tr or len(tr[key])<5: continue
        md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
        TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
        a=cleaned(tr[key]).copy(); tc=a[:,0]; t=np.arange(tc.min()-5,tc.max()+5,0.5); X=iss_ecef(t0+t)
        Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
        a[:,1]=a[:,1]+np.interp(tc,t,dfd)*PPS[key]*1e-3
        geo[key]=(a,tc,t,X,TX,RX,dfd)
def fit(p,hm,H=55,Nfix=None):
    obs=[];cols=[];iv=[];keys=[k for k in passes[p] if k in geo]
    for j,key in enumerate(keys):
        a,tc,t,X,TX,RX,dfd=geo[key]
        S=np.array([stec(TX,q,hm,H)+stec(q,RX,hm,H) for q in X]); dS=np.interp(tc,t,K*np.gradient(S,t))
        obs+=list(a[:,1]); iv+=list(dS); c=np.zeros((len(tc),len(keys))); c[:,j]=-np.interp(tc,t,dfd); cols.append(c)
    obs=np.array(obs); A=np.vstack(cols); iv=np.array(iv)
    if Nfix is None: A=np.column_stack([A,iv]); y=obs
    else: y=obs-Nfix*iv
    q,*_=np.linalg.lstsq(A,y,rcond=None); r=y-A@q
    dts={k:q[i]*1e3 for i,k in enumerate(keys)}; Nm=q[-1] if Nfix is None else Nfix
    dev=[dts[k]-O[k+('_1104' if '2305' in p else '_1102')]['dt_ms'] for k in keys]
    return dts,Nm,np.sqrt(np.mean(r**2)),dev
out={}; maxchg=0; rmsfree=[]; rmsfix=[]
for p in passes:
    base=fit(p,300)[3]
    for hm in [300,350,420,500,670]:
        d,Nm,rms,dev=fit(p,hm); ch=np.max(np.abs(np.array(dev)-base)); maxchg=max(maxchg,ch)
        if hm==300: rmsfree.append(rms)
        print(f"{p} hm={hm} free: Nm={Nm:.2e} foF2={np.sign(Nm)*np.sqrt(abs(Nm)/1.24e10):+.1f} rms={rms:.3f} dt-OEM={np.round(dev)} change={ch:.1f}")
    if p in iono:
        fo,hm=iono[p]; N=1.24e10*fo**2
        for h in [300,350,420,hm]:
            d,Nm,rms,dev=fit(p,h,Nfix=N); ch=np.max(np.abs(np.array(dev)-base)); maxchg=max(maxchg,ch); rmsfix.append(rms)
            print(f"{p}   FIXED ionosonde foF2={fo} hm={h}: rms={rms:.3f} dt-OEM={np.round(dev)} change={ch:.1f}")
print('max change',round(maxchg,1),'rms free',np.round([min(rmsfree),max(rmsfree)],2),'rms fixed',np.round([min(rmsfix),max(rmsfix)],2))
json.dump(dict(max_change_ms=float(maxchg),rms_free=[float(min(rmsfree)),float(max(rmsfree))],rms_fixed=[float(min(rmsfix)),float(max(rmsfix))]),open('../results/iono_sens_v2.json','w'))
