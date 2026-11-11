import numpy as np, sys, json
from scipy import signal as sg
from sgp4.propagation import gstime
args=sys.argv[1:]; sys.argv=['x','X']
base=open('station_run.py').read().split('print(f"== {ST} ==")')[0]
exec(base)
import os
from linefilt import maybe_clean
OUT=os.environ.get('OUT','../results/iono_tracks.json')
out={}
for fn in args:
    md,fs,t0,TX,RX,fb,bb,ts=load(fn)
    if os.environ.get('CLEAN'): bb=maybe_clean(bb,ts,fs)
    t=np.arange(-60,360,0.25); X=iss_ecef(t0+t); R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1); fdg=-np.gradient(R1+R2,t)/lam
    fd=np.interp(ts,t,fdg); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    q=int(fs//100); yd=sg.decimate(sg.decimate(y,q//4 if q>=8 else 1,ftype='fir'),4,ftype='fir') if q>=8 else y
    fsd=fs/((q//4)*4); 
    W=int(2.0*fsd); H=int(0.5*fsd); NF=4096
    fmax=0.9*min(fb,fs/2-fb)
    rows=[]
    for k in range(0,len(yd)-W,H):
        tc=(k+W/2)/fsd; fdc=np.interp(tc,t,fdg)
        if not (40<abs(fdc)<fmax): continue
        seg=yd[k:k+W]*np.hanning(W); S=np.abs(np.fft.fftshift(np.fft.fft(seg,NF)))**2; f=np.fft.fftshift(np.fft.fftfreq(NF,1/fsd))
        m=abs(f)<20; noise=np.median(S[(abs(f)>25)&(abs(f)<45)]) if ((abs(f)>25)&(abs(f)<45)).any() else np.median(S)
        i=np.argmax(S*m); snr=10*np.log10(S[i]/noise)
        el1=np.degrees(np.arcsin(np.dot(np.interp(tc,t,X[:,0])-TX[0],0)+0)) if False else None
        rows.append((tc,f[i],snr,float(np.interp(tc,t,R1)),float(np.interp(tc,t,R2)),float(fdc)))
    rows=np.array(rows); good=rows[rows[:,2]>12]
    key=fn.split('/')[-1][11:31]; out[key]=good.tolist()
    if len(good): print(f"{key}: {len(good)} pontos SNR>12 dB | Δf de {good[:,1].min():+.1f} a {good[:,1].max():+.1f} Hz | t {good[0,0]:.0f}-{good[-1,0]:.0f}s | Δf vs t: " + " ".join(f"{a:.0f}s:{b:+.1f}" for a,b in good[::max(1,len(good)//8),:2]))
    else: print(key,': sem pontos fortes')
json.dump(out,open(OUT,'w'))
