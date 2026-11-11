import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import json,numpy as np,datetime as dt,bwav,glob
from scipy import signal
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def ecef(lat,lon,h):
    a=6378.137; f=1/298.257223563; e2=f*(2-f); la,lo=np.radians(lat),np.radians(lon); N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.stack([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)],-1)
lam=299792.458/49.97e6
AC=[]
import os
for l in open(os.environ.get('ADSB','../data/adsb/adsb_be_nov03.jsonl')):
    d=json.loads(l); P=[p for p in d['pts'] if p[3] not in (None,'ground')]
    if len(P)<3: continue
    ts=np.array([p[0] for p in P],float); 
    alt=np.array([(p[6] if isinstance(p[6],(int,float)) else p[3]) for p in P],float)*0.3048/1000
    AC.append((d.get('r') or d.get('icao'),d.get('t'),ts,np.array([p[1] for p in P]),np.array([p[2] for p in P]),alt))
files=sorted(glob.glob(WAVDIR_+'*_10[0-2]?_BEOTTI*.wav'))
fig,ax=plt.subplots(len(files),1,figsize=(15,3.4*len(files)))
for a,p in zip(ax,files):
    md,x,pps=bwav.read(p); fs=md['samplerate']; t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>900)&(ff<1400); fb=ff[m][np.argmax(F[m])]
    f,t,S=signal.spectrogram(x,fs,nperseg=8192,noverlap=8192-1024); k=(f>fb-110)&(f<fb+110)
    Z=10*np.log10(S[k]); Z-=np.median(Z,axis=1,keepdims=True)
    a.pcolormesh(t,f[k]-fb,Z,vmin=0,vmax=30,shading='auto',cmap='Greys')
    tt=t0+np.arange(0,300,1.0); cands=[]
    for name,typ,ts,la,lo,al in AC:
        if ts.min()>tt[-1] or ts.max()<tt[0]: continue
        sel=(tt>=ts.min())&(tt<=ts.max())
        if sel.sum()<10: continue
        X=ecef(np.interp(tt[sel],ts,la),np.interp(tt[sel],ts,lo),np.interp(tt[sel],ts,al))
        R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1)
        fd=-np.gradient(R1+R2,tt[sel])/lam
        strength=-10*np.log10((R1*R2)**2).max()   # relativo
        cands.append((strength,name,typ,tt[sel]-t0,fd,(R1*R2).min()))
    cands.sort(key=lambda c:-c[0])
    cols=plt.cm.tab10(np.arange(10))
    for i,(s,name,typ,tl,fd,_) in enumerate(cands[:8]):
        a.plot(tl,fd,'--',color=cols[i],lw=1.4,label=f'{name} {typ}')
    a.set_ylim(-110,110); a.set_xlim(0,300); a.legend(fontsize=7,ncol=4,loc='lower left')
    a.set_title(f"{dt.datetime.fromtimestamp(t0,dt.timezone.utc):%Y-%m-%d %H:%M} UTC  Ottignies — tracejado: Doppler bistático previsto pelo ADS-B (8 aviões mais próximos)",fontsize=9)
    a.set_ylabel('Doppler (Hz)')
plt.tight_layout(); plt.savefig('aircraft_overlay.png',dpi=65)
