import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import numpy as np, json, sys, os
from scipy import signal as sg
exec(open('figstyle.py').read())
os.environ['ADSB']='../data/adsb/adsb_be_nov03.jsonl'
exec(open('air_overlay.py').read().split("files=")[0])   # AC, ecef, lam
import bwav
rows=np.load('../results/air_mf_rows.npy',allow_pickle=True)   # from air_mf.py
ok={'SP-RSW','9H-WAZ','EI-IKW'}
p=WAVDIR_+'RAD_BEDOUR_20261103_1010_BEOTTI_SYS001.wav'
md,x,pps=bwav.read(p); fs=md['samplerate']; t0=md['start_us']/1e6
TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>900)&(ff<1400); fb=ff[m][np.argmax(F[m])]
f,t,S=sg.spectrogram(x,fs,nperseg=8192,noverlap=8192-1024); k=(f>fb-75)&(f<fb+75)
Z=10*np.log10(S[k]); Z-=np.median(Z,axis=1,keepdims=True)
fig,ax=plt.subplots(1,2,figsize=(7.16,2.75),sharey=True)
for a in ax: a.pcolormesh(t,f[k]-fb,Z,vmin=0,vmax=28,shading='auto',cmap='Greys')
cols=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4']; i=0
tt=t0+np.arange(0,300,1.0)
for name,typ,ts,la,lo,al in AC:
    if name not in ok: continue
    sel=(tt>=ts.min())&(tt<=ts.max())
    if sel.sum()<10: continue
    X=ecef(np.interp(tt[sel],ts,la),np.interp(tt[sel],ts,lo),np.interp(tt[sel],ts,al))
    fd=-np.gradient(np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1),tt[sel])/lam
    ax[1].plot(tt[sel]-t0,fd,'--',color=cols[i%5],lw=1.2,label=f'{name} ({typ})'); i+=1
ax[0].set_title('(a) Ottignies, 3 Nov 2026, 10:10:00 UTC',fontsize=8); ax[1].set_title('(b) same, with ADS-B-predicted bistatic Doppler',fontsize=8)
ax[0].set_ylabel('Doppler (Hz)')
for a in ax: a.set_xlabel('Time (s)'); a.set_ylim(-75,75); a.set_xlim(0,300)
ax[1].legend(fontsize=7,loc='upper center',bbox_to_anchor=(0.5,-0.22),ncol=3); ax[1].annotate('no ADS-B match',xy=(262,-5),xytext=(205,-58),fontsize=7,color='#444',arrowprops=dict(arrowstyle='->',color='#666',lw=0.7))
plt.tight_layout(); plt.savefig(OUT+'fig_aircraft_example.pdf'); plt.savefig(OUT+'fig_aircraft_example.png',dpi=330)
print(sorted(ok))
