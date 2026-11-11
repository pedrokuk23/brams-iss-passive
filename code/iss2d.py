import numpy as np, sys
from sgp4.propagation import gstime
args=sys.argv[1:]; sys.argv=['x','X']
base=open('station_run.py').read().split('print(f"== {ST} ==")')[0]
exec(base)
for fn in args:
    md,fs,t0,TX,RX,fb,bb,ts=load(fn)
    t=np.arange(-60,360,0.5); X=iss_ecef(t0+t); Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fdg=-np.gradient(Rb,t)/lam
    fmax=0.9*min(fb,fs/2-fb); taus=np.arange(-6,6.01,0.5); dfs=np.arange(-12,12.01,0.5)
    q=2; from scipy import signal as sg
    bbd=sg.decimate(bb,q,ftype='fir'); fsd=fs/q; tsd=np.arange(len(bbd))/fsd; B=int(0.5*fsd); nb=len(bbd)//B
    Z=np.full((len(taus),len(dfs)),np.nan)
    for a,tau in enumerate(taus):
        fd=np.interp(tsd+tau,t,fdg); y0=bbd*np.exp(-2j*np.pi*np.cumsum(fd)/fsd); valid=(np.abs(fd)>30)&(np.abs(fd)<fmax); vb=valid[:nb*B].reshape(nb,B).all(1)
        for b,df in enumerate(dfs):
            y=y0*np.exp(-2j*np.pi*df*tsd); P=np.abs(y[:nb*B].reshape(nb,B).mean(1))**2; Z[a,b]=P[vb].mean()
    z=(Z-np.median(Z))/(1.4826*np.median(abs(Z-np.median(Z))))
    i,j=np.unravel_index(np.nanargmax(z),z.shape)
    print(f"{fn.split('/')[-1][15:31]}: pico z={z[i,j]:.1f} em tau={taus[i]:+.1f}s, Δf={dfs[j]:+.1f} Hz")
