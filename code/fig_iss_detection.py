# ISS detection figure alone (same code as paper_figs2.py), needs only the Ottignies file of 3 Nov 23:55
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import numpy as np, sys, json, glob, datetime as dt
from scipy import signal as sg
from sgp4.propagation import gstime
from sgp4.api import Satrec, jday
sys.argv=['x','X']
base=open('station_run.py').read().split('print(f"== {ST} ==")')[0]
exec(base)
import bwav
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('figstyle.py').read())
D=WAVDIR_
# ---- Fig ISS detection: 2D z-map + dechirped spectrogram (Ottignies 23:55)
fn=D+'RAD_BEDOUR_20261103_2355_BEOTTI_SYS001.wav'
md,fs,t0,TX,RX,fb,bb,ts=load(fn)
t=np.arange(-60,360,0.5); X=iss_ecef(t0+t); Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fdg=-np.gradient(Rb,t)/lam
q=1; bbd=bb; fsd=fs; tsd=np.arange(len(bbd))/fsd; B=int(0.5*fsd); nb=len(bbd)//B
fmax=0.9*min(fb,fs/2-fb); taus=np.arange(-6,6.01,0.5); dfs=np.arange(-12,12.01,0.5); Z=np.zeros((len(taus),len(dfs)))
for a,tau in enumerate(taus):
    fd=np.interp(tsd+tau,t,fdg); y0=bbd*np.exp(-2j*np.pi*np.cumsum(fd)/fsd); vb=((np.abs(fd)>30)&(np.abs(fd)<fmax))[:nb*B].reshape(nb,B).all(1)
    for b,df in enumerate(dfs):
        y=y0*np.exp(-2j*np.pi*df*tsd); Z[a,b]=(np.abs(y[:nb*B].reshape(nb,B).mean(1))**2)[vb].mean()
Z=(Z-np.median(Z))/(1.4826*np.median(abs(Z-np.median(Z))))
fig,ax=plt.subplots(1,2,figsize=(7.16,2.4),gridspec_kw={'width_ratios':[1,1.5]})
im=ax[0].pcolormesh(dfs,taus,Z,shading='auto',cmap='Blues'); cb=plt.colorbar(im,ax=ax[0]); cb.set_label('Detection statistic $z$')
i,j=np.unravel_index(Z.argmax(),Z.shape); ax[0].plot(dfs[j],taus[i],'+',color='#eb6834',ms=10,mew=1.5)
ax[0].set_xlabel(r'Frequency offset $\nu$ (Hz)'); ax[0].set_ylabel(r'Time offset $\tau$ (s)'); ax[0].set_title(f'(a) peak $z$ = {Z.max():.0f} at $\\tau$ = {taus[i]:.1f} s, $\\nu$ = {dfs[j]:.1f} Hz',fontsize=8)
fd=np.interp(tsd,t,fdg); y=bbd*np.exp(-2j*np.pi*(np.cumsum(fd)/fsd))
yd=sg.decimate(sg.decimate(sg.decimate(y,4,ftype="fir"),4,ftype="fir"),2,ftype="fir"); f2s=fsd/32
f2,t2,S2=sg.spectrogram(yd,f2s,nperseg=256,noverlap=240,return_onesided=False); f2=np.fft.fftshift(f2); S2=np.fft.fftshift(S2,0)
k=abs(f2)<25; Zs=10*np.log10(S2[k]/np.median(S2[k],axis=0,keepdims=True))
valid=(np.abs(fd)>30)&(np.abs(fd)<fmax); tv=tsd[valid]
ax[1].pcolormesh(t2,f2[k],Zs,vmin=0,vmax=14,shading='auto',cmap='Greys'); ax[1].set_xlim(tv.min(),tv.max())
ax[1].set_xlabel('Time since 23:55:00 UTC (s)'); ax[1].set_ylabel('Residual frequency (Hz)'); ax[1].set_title('(b) ISS echo after removing the TLE-predicted Doppler',fontsize=8)
plt.tight_layout(); plt.savefig(OUT+'fig_iss_detection.pdf'); plt.savefig(OUT+'fig_iss_detection.png')
print('ok', Z.max(), taus[i], dfs[j])
