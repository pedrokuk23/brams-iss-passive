# Varredura do catálogo com busca em frequência (Δf via FFT por bloco), τ=0 e controle τ=+30 s
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import numpy as np, sys, datetime as dt, bwav, json
from scipy import signal
from sgp4.api import Satrec, SatrecArray, jday
from sgp4.propagation import gstime
fn=sys.argv[1]
def ecef(lat,lon,h):
    a=6378.137; f=1/298.257223563; e2=f*(2-f); la,lo=np.radians(lat),np.radians(lon); N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
lam=299792.458/49.97e6
md,x,pps=bwav.read(fn); fs=md['samplerate']; t0=md['start_us']/1e6
TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
n=len(x); ts=np.arange(n)/fs
bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); bb=bb-signal.oaconvolve(bb,np.ones(int(fs))/int(fs),mode='same')
bb=signal.decimate(bb,2,ftype='fir'); fsd=fs/2; tsd=np.arange(len(bb))/fsd
fmax=0.9*min(fb,fs/2-fb)
L=[l.rstrip() for l in open('../data/tle/leo.txt')]
names=[];sats=[];seen=set()
for k in range(0,len(L)-2,3):
    if L[k+1].startswith('1 ') and L[k+2].startswith('2 ') and L[k+1][2:7] not in seen: seen.add(L[k+1][2:7]); names.append(L[k].strip()); sats.append(Satrec.twoline2rv(L[k+1],L[k+2]))
arr=SatrecArray(sats)
tg=np.arange(-40,341,1.0); d0=dt.datetime.fromtimestamp(t0,dt.timezone.utc)
jd0,fr0=jday(d0.year,d0.month,d0.day,d0.hour,d0.minute,d0.second+d0.microsecond/1e6)
jd=np.full(len(tg),jd0); fr=fr0+tg/86400
e,r,v=arr.sgp4(jd,fr); g=np.array([gstime(a+b) for a,b in zip(jd,fr)])
X=np.stack([r[...,0]*np.cos(g)+r[...,1]*np.sin(g),-r[...,0]*np.sin(g)+r[...,1]*np.cos(g),r[...,2]],-1)
R1=np.linalg.norm(X-TX,axis=-1); R2=np.linalg.norm(X-RX,axis=-1)
up=((X-RX)@(RX/np.linalg.norm(RX)))>0
fd=-np.gradient(R1+R2,tg,axis=1)/lam
cand=np.where(((np.abs(fd)>30)&(np.abs(fd)<fmax)&up&(e==0)).sum(1)>=20)[0]
B=int(0.5*fsd); nb=len(bb)//B; fax=np.arange(-40,40.01,1.0); kk=np.abs(fax)<=12
E=np.exp(-2j*np.pi*np.outer(np.arange(B)/fsd,fax)).astype(np.complex64); tb=(np.arange(nb)*B/fsd)
def stat(i,tau):
    f=np.interp(tsd+tau,tg,fd[i]); y=bb*np.exp(-2j*np.pi*np.cumsum(f)/fsd)
    valid=((np.abs(f)>30)&(np.abs(f)<fmax)&(np.interp(tsd+tau,tg,up[i].astype(float))>0.5))[:nb*B].reshape(nb,B).all(1)
    if valid.sum()<20: return None
    Y=y[:nb*B].reshape(nb,B)[valid].astype(np.complex64)
    S=np.abs(Y@E)**2                     # DFT por bloco de 0,5 s em Δf de -40 a 40 Hz, passo 1 Hz
    P=S.mean(0)
    return P[kk].max()/np.median(P), float(fax[kk][np.argmax(P[kk])]), int(valid.sum())
out=[]
for i in cand:
    a=stat(i,0.0); c=stat(i,30.0)
    if a is None or c is None: continue
    rr=float(np.sqrt((R1[i]*R2[i]).min()))
    out.append(dict(file=fn.split('/')[-1][11:31],obj=names[i],peak=float(a[0]),df=a[1],nb=a[2],ctrl=float(c[0]),rr=rr))
json.dump(out,open('../results/catscan/cat3_'+fn.split('/')[-1][11:31]+'.json','w'))
print(fn.split('/')[-1][11:31],len(cand),'candidatos',len(out),'avaliados',flush=True)
