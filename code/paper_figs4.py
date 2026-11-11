import numpy as np, json, sys
exec(open('iono_fit.py').read().split("passes={")[0])
exec(open('figstyle.py').read())
P={'(a) 3 Nov, 02:20 UTC':['20261103_0220_BEOTTI','20261103_0220_BEOUDS'],'(b) 3 Nov, 23:55 UTC':['20261103_2355_BEOTTI','20261103_2355_BEOUDS'],'(c) 4 Nov, 23:05 UTC':['20261104_2305_BEOTTI','20261104_2305_BEOUDS']}
import json as js
tr=js.load(open('../results/iono_tracks.json'))
fig,ax=plt.subplots(1,3,figsize=(7.16,2.3))
POS=[(0.55,0.04,'bottom'),(0.03,0.97,'top'),(0.30,0.97,'top')]
for a,(pn,keys),pos in zip(ax,P.items(),POS):
    data=[]
    for key in keys:
        if key not in tr or len(tr[key])<5: continue
        md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
        TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
        o=cleaned(tr[key]); t=np.arange(o[:,0].min()-8,o[:,0].max()+8,0.5); X=iss_ecef(t0+t)
        Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
        S=np.array([stec(TX,p,1.0)+stec(p,RX,1.0) for p in X]); dS=np.gradient(S,t)
        k=slice(6,-6); data.append((key.split('_')[2],o,t[k],-dfd[k],K*dS[k]))
    A=np.vstack([np.column_stack([np.interp(o[:,0],tt,mt),np.interp(o[:,0],tt,mi)]) for _,o,tt,mt,mi in data]); y=np.concatenate([o[:,1] for _,o,*_ in data])
    p,*_=np.linalg.lstsq(A,y,rcond=None); rms=np.sqrt(np.mean((y-A@p)**2))
    for st,o,tt,mt,mi in data:
        a.plot(o[:,0],o[:,1],MK[st],color=COL[st],ms=2.2,ls='',label=NAME[st])
        a.plot(tt,p[0]*mt+p[1]*mi,'-',color=COL[st],lw=1.2)
        a.plot(tt,p[0]*mt,':',color=COL[st],lw=0.9)
    a.set_title(pn,fontsize=8)
    a.text(pos[0],pos[1],f"$\\delta t$={p[0]*1e3:+.0f} ms\n$N_mF_2$={p[1]/1e11:.1f}$\\times10^{{11}}$ m$^{{-3}}$\nrms {rms:.2f} Hz",transform=a.transAxes,fontsize=6.3,va=pos[2])
    a.set_xlabel('Time in file (s)')
    for s in ['top','right']: a.spines[s].set_visible(False)
ax[0].set_ylabel('Residual Doppler (Hz)'); ax[1].legend(fontsize=6.3,loc='lower right',markerscale=2)
plt.tight_layout(); plt.savefig(OUT+'fig_residual.pdf'); plt.savefig(OUT+'fig_residual.png')
