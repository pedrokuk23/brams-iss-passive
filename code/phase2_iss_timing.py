# ISS timing in Phase 2: the analysis of Sections 6-7 (time-resolved residual Doppler, joint
# fit of the timing term and a Chapman layer, time-base correction from the file's PPS record) applied to every Phase 2
# recording at hand in which the blind search detected the ISS complex, and compared with the archived NASA ephemerides
# (OEM files of 19 Oct - 6 Nov 2026 from nasa-public-data.s3.amazonaws.com/iss-coords/<date>/ISS_OEM/, data/oem_archive/).
# Element set: the ISS (ZARYA) set of the Phase 2 history with the latest epoch before the recording (a prediction).
# Reference: delta t_OEM from the latest OEM issued before the pass (a prediction, as in Phase 1) and, where a later file
# also covers the pass, from that file.
import os, sys, json, glob, datetime as dtm, numpy as np
from scipy import signal as sg
from sgp4.api import Satrec, jday
from sgp4.propagation import gstime
from astropy.utils import iers; iers.conf.auto_download=False; iers.conf.iers_degraded_accuracy='warn'
from astropy.time import Time
from astropy.coordinates import GCRS, ITRS, CartesianRepresentation
import astropy.units as u
import bwav
from linefilt import maybe_clean
WD=os.environ.get('P2DIR','../data/raw_wav_phase2')
path={os.path.basename(p)[11:31]:p for d in WD.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
if not path: sys.exit('phase2_iss_timing: no WAV files found in '+WD)
lam=299792.458/49.97e6; c=299792458.0; f0=49.97e6; K=40.3/(c*f0); RE=6371.0
def ecef(lat,lon,h):
    a=6378.137; f=1/298.257223563; e2=f*(2-f); la,lo=np.radians(lat),np.radians(lon); N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.stack([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)],-1)
# element sets
ISS=[]
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ISS (ZARYA)') and L[k+1].startswith('1 '): s=Satrec.twoline2rv(L[k+1],L[k+2]); ISS.append((s.jdsatepoch+s.jdsatepochF,L[k+1],L[k+2]))
ISS=sorted(set(ISS))
if not ISS: sys.exit('phase2_iss_timing: no ISS element sets in data/tle/phase2 (unzip auxiliary_data_phase2.zip)')
def sat_before(t0):
    jd=t0/86400+2440587.5; c_=[e for e in ISS if e[0]<=jd]; e=c_[-1]; return Satrec.twoline2rv(e[1],e[2]),(jd-e[0])*24
def tle_ecef(sat,ts):
    out=[]
    for t in ts:
        d=dtm.datetime.fromtimestamp(t,dtm.timezone.utc); jd,fr=jday(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond/1e6)
        e,r,v=sat.sgp4(jd,fr); g=gstime(jd+fr); out.append([r[0]*np.cos(g)+r[1]*np.sin(g),-r[0]*np.sin(g)+r[1]*np.cos(g),r[2]])
    return np.array(out)
# OEM archive
def readoem(fn):
    T=[];S=[];cr=None
    for l in open(fn):
        if l.startswith('CREATION_DATE'): cr=dtm.datetime.fromisoformat(l.split('=')[1].strip()).replace(tzinfo=dtm.timezone.utc).timestamp()
        p=l.split()
        if len(p)==7 and p[0][:2]=='20' and 'T' in p[0]:
            T.append(dtm.datetime.fromisoformat(p[0]).replace(tzinfo=dtm.timezone.utc).timestamp()); S.append([float(v) for v in p[1:]])
    return cr,np.array(T),np.array(S)
OEMS=sorted([(os.path.basename(f)[8:21],)+readoem(f) for f in glob.glob('../data/oem_archive/ISS_OEM_*.txt')],key=lambda x:x[1])
def lagrange(T,S,t,npt=10):
    # npt-point Lagrange interpolation of the OEM positions (standard practice for OEM files; a cubic
    # Hermite on the 4-min nodes moves delta t_OEM by up to ~12 ms)
    out=np.empty((len(t),3))
    for j,x in enumerate(t):
        i=int(np.clip(np.searchsorted(T,x)-npt//2,0,len(T)-npt)); idx=np.arange(i,i+npt); u_=(T[idx]-x)/240.
        w=np.array([np.prod([(0-u_[m])/(u_[k]-u_[m]) for m in range(npt) if m!=k]) for k in range(npt)])
        out[j]=w@S[idx,:3]
    return out
def oem_ecef(T,S,t):
    r=lagrange(T,S,t); tt=Time(t,format='unix',scale='utc'); it=GCRS(CartesianRepresentation(r.T*u.km),obstime=tt).transform_to(ITRS(obstime=tt))
    return np.array([it.x.to(u.km).value,it.y.to(u.km).value,it.z.to(u.km).value]).T
def chapman(h,Nm,hm=300.,H=55.):
    z=(h-hm)/H; return Nm*np.exp(0.5*(1-z-np.exp(-z)))
def stec(A,B,Nm):
    s=np.linspace(0,1,400)[:,None]; P=A[None,:]+s*(B-A)[None,:]; h=np.linalg.norm(P,axis=1)-RE; L=np.linalg.norm(B-A)*1e3
    return np.trapezoid(chapman(h,Nm),dx=1/399)*L
S2=json.load(open('../results/phase2_summary.json'))
keys=sorted({d['file'] for d in S2['detections'] if d['class']=='ISS complex' and d['file'] in path})
NMIN=int(os.environ.get('NMIN',19))
def proc(key):
    md,x,pps=bwav.read(path[key]); fs=md['samplerate']; t0=md['start_us']/1e6
    p=np.array(pps,dtype=float); off0=0.0
    if len(p)>=10:
        kk=p[:,0]; tt=(p[:,1]-md['start_us'])/1e6; (a,b),*_=np.linalg.lstsq(np.vstack([np.ones_like(kk),kk]).T,tt,rcond=None); off0=-a
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
    ts=np.arange(len(x))/fs; bb=sg.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-sg.oaconvolve(bb,np.ones(w)/w,mode='same')
    if not key.endswith('BEOTTI'): bb=maybe_clean(bb,ts,fs)
    sat,age_h=sat_before(t0)
    t=np.arange(-60,360,0.25); X=tle_ecef(sat,t0+t); R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1); fdg=-np.gradient(R1+R2,t)/lam
    fd=np.interp(ts,t,fdg); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    q=int(fs//100); yd=sg.decimate(sg.decimate(y,q//4,ftype='fir'),4,ftype='fir'); fsd=fs/((q//4)*4)
    W=int(2.0*fsd); H=int(0.5*fsd); NF=4096; fmax=0.9*min(fb,fs/2-fb); rows=[]
    for k in range(0,len(yd)-W,H):
        tc=(k+W/2)/fsd; fdc=np.interp(tc,t,fdg)
        if not (40<abs(fdc)<fmax): continue
        seg=yd[k:k+W]*np.hanning(W); Sp=np.abs(np.fft.fftshift(np.fft.fft(seg,NF)))**2; f=np.fft.fftshift(np.fft.fftfreq(NF,1/fsd))
        mm=abs(f)<20; ring=(abs(f)>25)&(abs(f)<45); noise=np.median(Sp[ring]); i=np.argmax(Sp*mm); snr=10*np.log10(Sp[i]/noise)
        rows.append((tc,f[i],snr))
    rows=np.array(rows); g=rows[rows[:,2]>12] if len(rows) else rows
    if len(g)<NMIN: print(key,'few points',len(g),flush=True); return key,dict(n=int(len(g)),fit=False)
    med=np.median(g[:,1]); g=g[abs(g[:,1]-med)<4]
    if len(g)<NMIN: return key,dict(n=int(len(g)),fit=False)
    tc=g[:,0]; tt=np.arange(tc.min()-5,tc.max()+5,0.5); Xt=tle_ecef(sat,t0+tt)
    fT=lambda XX: -np.gradient(np.linalg.norm(XX-TX,axis=1)+np.linalg.norm(XX-RX,axis=1),tt)/lam
    ft=fT(Xt); dfd=np.gradient(ft,tt); Sx=np.array([stec(TX,pp,1.0)+stec(pp,RX,1.0) for pp in Xt]); dS=np.gradient(Sx,tt)
    mt=-np.interp(tc,tt,dfd); obs=g[:,1]-mt*off0; mi=K*np.interp(tc,tt,dS)
    A=np.column_stack([mt,mi]); pr,*_=np.linalg.lstsq(A,obs,rcond=None); r=obs-A@pr; rms=np.sqrt(np.mean(r**2))
    C=np.linalg.inv(A.T@A)*rms**2*len(obs)/max(len(obs)-2,1)
    p1=float(np.dot(mt,obs)/np.dot(mt,mt))                                 # timing term alone, for comparison
    res0=dict(dt_noiono_ms=p1*1e3,rms_noiono=float(np.sqrt(np.mean((obs-mt*p1)**2))),corr=float(np.corrcoef(mt,mi)[0,1]))
    res=dict(n=int(len(obs)),dt_ms=float(pr[0]*1e3),sdt_ms=float(np.sqrt(C[0,0])*1e3),Nm=float(pr[1]),rms=float(rms),pps_ms=float(off0*1e3),tle_age_h=float(age_h),oem={})
    msk=(tt>=tc.min())&(tt<=tc.max()); tp=t0+0.5*(tc.min()+tc.max())
    for lab,cr,T,S in OEMS:
        if t0+tt[0]<T[0] or t0+tt[-1]>T[-1]: continue
        Xo=oem_ecef(T,S,t0+tt); fo=fT(Xo); yy=(fo-ft)[msk]; Aa=-dfd[msk]; d_=float(np.dot(Aa,yy)/np.dot(Aa,Aa))
        res['oem'][lab]=dict(dt_ms=d_*1e3,issued_before=bool(cr<tp),age_h=float((tp-cr)/3600))
    # quasi-definitive reference from the element sets: the first ISS set with epoch after the pass, whose fit span
    # contains it, and the set after that (their difference measures the reference's own scatter)
    jdp=tp/86400+2440587.5; post=[e for e in ISS if e[0]>jdp+1/24][:2]; res['tle_post']={}
    for j,e in enumerate(post):
        Xp=tle_ecef(Satrec.twoline2rv(e[1],e[2]),t0+tt); fp=fT(Xp); yy=(fp-ft)[msk]; Aa=-dfd[msk]
        res['tle_post'][str(j)]=dict(dt_ms=float(np.dot(Aa,yy)/np.dot(Aa,Aa)*1e3),epoch_after_h=float((e[0]-jdp)*24))
    for lab in res['oem']: res['oem'][lab]['prop_age_h']=float((tp-[o for o in OEMS if o[0]==lab][0][2][0])/3600)
    res['t_pass']=float(tp); res['arc_s']=[float(tc.min()),float(tc.max())]; res.update(res0); res['track']=g.tolist()
    print(key,res['n'],'%.1f±%.1f'%(res['dt_ms'],res['sdt_ms']),{k:round(v['dt_ms'],1) for k,v in res['oem'].items()},flush=True)
    return key,res
if __name__=='__main__':
    from multiprocessing import Pool
    with Pool(int(os.environ.get('NPROC',os.cpu_count() or 1))) as pool: out=dict(pool.map(proc,keys,chunksize=1))
    json.dump(out,open('../results/phase2_iss_timing.json','w'),indent=1)
