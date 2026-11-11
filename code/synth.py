# Validação do estimador com ecos sintéticos injetados em gravações reais (sem ISS) e distribuição nula do detector.
import numpy as np, json, sys, glob, itertools, os
from scipy import signal as sg
sys.argv=['x','X']
exec(open('iono_fit.py').read().split("passes={")[0])   # ecef, lam, iss_ecef, bwav, K, stec, chapman, cleaned, tr(iono_tracks)
D=os.environ.get('WAVDIR','../data/raw_wav/')
def baseband(x,fs):
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
    ts=np.arange(len(x))/fs; bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-signal.oaconvolve(bb,np.ones(w)/w,mode='same')
    return fb,bb,ts
class Geo:
    def __init__(s,key):
        md,_,_=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); s.t0=md['start_us']/1e6; s.fs=md['samplerate']
        s.TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); s.RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
        s.t=np.arange(-60,360,0.25); s.X=iss_ecef(s.t0+s.t)
        s.R1=np.linalg.norm(s.X-s.TX,axis=1); s.R2=np.linalg.norm(s.X-s.RX,axis=1); s.fd=-np.gradient(s.R1+s.R2,s.t)/lam
        s.fdot=np.gradient(s.fd,s.t)
        tq=np.arange(-60,360,1.0); Xq=iss_ecef(s.t0+tq); S=np.array([stec(s.TX,p,1.0)+stec(p,s.RX,1.0) for p in Xq])
        s.ion_unit=np.interp(s.t,tq,K*np.gradient(S,tq))      # Hz por (el/m^3) de NmF2
def search(bb,ts,fs,fb,g,taus=np.arange(-6,6.01,0.5),dfs=np.arange(-12,12.01,0.5),Tb=0.5):
    fmax=0.9*min(fb,fs/2-fb); B=int(Tb*fs); nb=len(bb)//B; E=np.exp(-2j*np.pi*np.outer(np.arange(B)/fs,dfs)).astype(np.complex64)
    Z=np.full((len(taus),len(dfs)),np.nan)
    for a,tau in enumerate(taus):
        fd=np.interp(ts+tau,g.t,g.fd); y0=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
        valid=(np.abs(fd)>30)&(np.abs(fd)<fmax); vb=valid[:nb*B].reshape(nb,B).all(1)
        if vb.sum()<5: continue
        Y=y0[:nb*B].reshape(nb,B)[vb].astype(np.complex64); Z[a]=(np.abs(Y@E/B)**2).mean(0)
    v=Z[np.isfinite(Z)]; med=np.median(v); z=(Z-med)/(1.4826*np.median(np.abs(v-med)))
    i,j=np.unravel_index(np.nanargmax(z),z.shape); return float(z[i,j]),float(taus[i]),float(dfs[j])
def track(bb,ts,fs,fb,g):
    fd=np.interp(ts,g.t,g.fd); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    q=int(fs//100); yd=sg.decimate(sg.decimate(y,q//4,ftype='fir'),4,ftype='fir'); fsd=fs/((q//4)*4)
    W=int(2.0*fsd); H=int(0.5*fsd); NF=4096; fmax=0.9*min(fb,fs/2-fb); f=np.fft.fftshift(np.fft.fftfreq(NF,1/fsd)); rows=[]
    for k in range(0,len(yd)-W,H):
        tc=(k+W/2)/fsd; fdc=np.interp(tc,g.t,g.fd)
        if not (40<abs(fdc)<fmax): continue
        S=np.abs(np.fft.fftshift(np.fft.fft(yd[k:k+W]*np.hanning(W),NF)))**2
        m=abs(f)<20; noise=np.median(S[(abs(f)>25)&(abs(f)<45)]); i=np.argmax(S*m); snr=10*np.log10(S[i]/noise)
        if snr>12: rows.append((tc,f[i]))
    return np.array(rows)
def fit(rows,g):
    if len(rows)<5: return None
    a=cleaned(rows); tc=a[:,0]
    A=np.column_stack([-np.interp(tc,g.t,g.fdot),np.interp(tc,g.t,g.ion_unit)]); p,*_=np.linalg.lstsq(A,a[:,1],rcond=None)
    return p[0]*1e3, p[1], len(a)
def inject(x,fs,fb,g,dt_true,Nm,amp,win):
    ts=np.arange(len(x))/fs; fd=np.interp(ts-dt_true,g.t,g.fd)+Nm*np.interp(ts,g.t,g.ion_unit)
    env=np.zeros_like(ts); m=(ts>=win[0])&(ts<=win[1]); env[m]=np.sin(np.pi*(ts[m]-win[0])/(win[1]-win[0]))**2
    ph=2*np.pi*np.cumsum(fb+fd)/fs+np.random.uniform(0,2*np.pi)
    return x+amp*env*np.cos(ph)
if __name__=='__main__':
    mode=sys.argv[-1] if len(sys.argv)>1 else ''
