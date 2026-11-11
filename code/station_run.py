import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import json,numpy as np,datetime as dt,bwav,glob,sys
from scipy import signal
from sgp4.api import Satrec, jday
from sgp4.propagation import gstime
ST=sys.argv[1]; DAY=sys.argv[2] if len(sys.argv)>2 else '20261103'
D=WAVDIR_
def ecef(lat,lon,h):
    a=6378.137; f=1/298.257223563; e2=f*(2-f); la,lo=np.radians(lat),np.radians(lon); N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.stack([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)],-1)
lam=299792.458/49.97e6
L=[l.strip() for l in open('../data/tle/iss.txt')]
for k in range(0,len(L)-2,3):
    if L[k].startswith('ISS (ZARYA)'): sat=Satrec.twoline2rv(L[k+1],L[k+2]); break
def iss_ecef(ts):
    out=[]
    for t in ts:
        d=dt.datetime.fromtimestamp(t,dt.timezone.utc); jd,fr=jday(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond/1e6)
        e,r,v=sat.sgp4(jd,fr); g=gstime(jd+fr)   # GMST (IAU-82), como exige o TEME do SGP4
        out.append([r[0]*np.cos(g)+r[1]*np.sin(g),-r[0]*np.sin(g)+r[1]*np.cos(g),r[2]])
    return np.array(out)
def load(fn):
    md,x,pps=bwav.read(fn); fs=md['samplerate']; t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
    n=len(x); ts=np.arange(n)/fs
    bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-signal.oaconvolve(bb,np.ones(w)/w,mode='same')
    return md,fs,t0,TX,RX,fb,bb,ts
def scan(bb,ts,fs,fdfun,taus,B=0.5,fmin=30,fmax=1100):
    sc=[]
    for tau in taus:
        fd=fdfun(ts+tau); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs); valid=(np.abs(fd)>fmin)&(np.abs(fd)<fmax)
        Bn=int(B*fs); nb=len(y)//Bn; P=np.abs(y[:nb*Bn].reshape(nb,Bn).mean(1))**2; vb=valid[:nb*Bn].reshape(nb,Bn).all(1)
        sc.append(P[vb].mean() if vb.sum()>5 else np.nan)
    sc=np.array(sc); z=(sc-np.nanmedian(sc))/(1.4826*np.nanmedian(abs(sc-np.nanmedian(sc)))); return z
print(f"== {ST} ==")
# ISS
for fnk in ['0040','0220','2355']:
    fn=D+f'RAD_BEDOUR_{DAY}_{fnk}_{ST}_SYS001.wav'
    md,fs,t0,TX,RX,fb,bb,ts=load(fn)
    tg=np.arange(-60,360,0.5); X=iss_ecef(t0+tg); Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fdg=-np.gradient(Rb,tg)/lam
    taus=np.arange(-30,30.01,0.25); z=scan(bb,ts,fs,lambda t: np.interp(t,tg,fdg),taus,fmax=0.9*min(fb,fs/2-fb))
    i=np.nanargmax(z); print(f"ISS {fnk}: pico z={z[i]:.1f} em tau={taus[i]:+.2f}s  (fs={fs:.0f}, beacon {fb:.1f} Hz, z(-3.5)={z[np.argmin(abs(taus+3.5))]:.1f})")
