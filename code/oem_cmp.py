import numpy as np, json, sys, datetime as dt
from sgp4.propagation import gstime
from astropy.utils import iers; iers.conf.auto_download=False; iers.conf.iers_degraded_accuracy='warn'
from astropy.time import Time
from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, CartesianDifferential
import astropy.units as u
sys.argv=['x','X']
exec(open('station_run.py').read().split('print(f"== {ST} ==")')[0])
import bwav
def readoem(fn):
    T=[];S=[]
    for l in open(fn):
        p=l.split()
        if not (len(p)==7 and p[0][:2]=='20' and 'T' in p[0]): continue
        T.append(dt.datetime.fromisoformat(p[0]).replace(tzinfo=dt.timezone.utc).timestamp()); S.append([float(v) for v in p[1:]])
    return np.array(T),np.array(S)
def lagrange(T,S,t,npt=10):
    # npt-point Lagrange interpolation of the OEM positions (a cubic Hermite on the 4-min nodes moved
    # delta t_OEM by up to ~12 ms)
    out=np.empty((len(t),3))
    for j,x in enumerate(t):
        i=int(np.clip(np.searchsorted(T,x)-npt//2,0,len(T)-npt)); idx=np.arange(i,i+npt); u_=(T[idx]-x)/240.
        w=np.array([np.prod([(0-u_[m])/(u_[k]-u_[m]) for m in range(npt) if m!=k]) for k in range(npt)])
        out[j]=w@S[idx,:3]
    return out
def oem_ecef(T,S,t):
    r=lagrange(T,S,t); tt=Time(t,format='unix',scale='utc')
    g=GCRS(CartesianRepresentation(r.T*u.km),obstime=tt); it=g.transform_to(ITRS(obstime=tt))
    return np.array([it.x.to(u.km).value,it.y.to(u.km).value,it.z.to(u.km).value]).T
# full files from the NASA archive (data/oem_archive); data/oem/ holds only the excerpts around the passes, too short
# for a 10-point interpolation at their ends
O2=readoem('../data/oem_archive/ISS_OEM_20261102T1517.txt'); O4=readoem('../data/oem_archive/ISS_OEM_20261104T1610.txt')
import os
tr=json.load(open(os.environ.get('TRACKS','../results/iono_tracks.json'))); D=os.environ.get('WAVDIR','../data/raw_wav/')
keys=sorted(tr) if os.environ.get('TRACKS') else ['20261103_0040_BEOTTI','20261103_0220_BEOTTI','20261103_0220_BEOUDS','20261103_2355_BEOTTI','20261103_2355_BEOUDS','20261104_0130_BEOUDS','20261104_0130_BEOTTI','20261104_2305_BEOUDS','20261104_2305_BEOTTI']
out={}
for key in keys:
    if key not in tr or len(tr[key])<5: continue
    md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    a=np.array(tr[key]); tc=a[:,0]; t=np.arange(tc.min()-5,tc.max()+5,0.5)
    Xt=iss_ecef(t0+t)
    for lab,(T,S) in [('1102',O2),('1104',O4)]:
        if t0+t[0]<T[0] or t0+t[-1]>T[-1]: continue
        Xo=oem_ecef(T,S,t0+t)
        f=lambda X: -np.gradient(np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1),t)/lam
        ft,fo=f(Xt),f(Xo); dfd=np.gradient(ft,t)
        m=(t>=tc.min())&(t<=tc.max())
        y=(fo-ft)[m]; A=-dfd[m]; dtt=np.dot(A,y)/np.dot(A,A); rms=np.sqrt(np.mean((y-A*dtt)**2))
        v=np.gradient(Xo,t,axis=0); vh=v/np.linalg.norm(v,axis=1)[:,None]
        al=np.sum((Xt-Xo)*vh,1)[m].mean(); dpos=np.linalg.norm(Xt-Xo,axis=1)[m].mean()
        out[key+'_'+lab]=dict(dt_ms=dtt*1e3,rms=rms,along_km=al,along_ms=al/np.linalg.norm(v,axis=1)[m].mean()*1e3,dpos_km=dpos)
        print(f"{key} OEM{lab}: δt_OEM={dtt*1e3:+7.1f} ms (fit rms {rms:.3f} Hz) | TLE-OEM along {al*1e3:+6.0f} m = {out[key+'_'+lab]['along_ms']:+6.1f} ms | |Δr| {dpos*1e3:.0f} m")
json.dump(out,open(os.environ.get('OUT','../results/oem_cmp.json'),'w'),indent=1)
