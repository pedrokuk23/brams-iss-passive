exec(open('bb2.py').read())
exec(open('figstyle.py').read())
from scipy import signal as sg
from linefilt import maybe_clean
import os
DL=os.environ.get('BBDIR','../data/raw_wav/')
fig,ax=plt.subplots(1,3,figsize=(7.16,2.3))
info={}
for a,st in zip(ax,['BEOTTI','BEOUDS','FRMAUB']):
    f=DL+f'RAD_BEDOUR_20261105_0310_{st}_SYS001.wav'
    md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs); g=GeoSat(f,'SPACEMOBILE-008')
    fd=np.interp(ts-0.5*0,g.t,g.fd); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    q=int(fs//100); yd=sg.decimate(sg.decimate(y,q//4,ftype='fir'),4,ftype='fir'); fsd=fs/((q//4)*4)
    f_,t_,S=sg.spectrogram(yd,fsd,nperseg=int(2*fsd),noverlap=int(1.75*fsd),nfft=1024,return_onesided=False,window='hann')
    f_=np.fft.fftshift(f_); S=np.fft.fftshift(S,axes=0); k=np.abs(f_)<=30
    Z=10*np.log10(S[k]+1e-12); Z-=np.median(Z,axis=0,keepdims=True)
    a.pcolormesh(t_,f_[k],Z,vmin=0,vmax=15,cmap='Greys',shading='auto')
    fdt=np.interp(t_,g.t,g.fd); bad=(np.abs(fdt)<30)|(np.abs(fdt)>0.9*min(fb,fs/2-fb))
    for i0 in np.where(np.diff(bad.astype(int))!=0)[0]: a.axvline(t_[i0],color='#c0392b',lw=0.5,ls=':')
    a.set_xlim(15,105); a.set_title(f'({"abc"[["BEOTTI","BEOUDS","FRMAUB"].index(st)]}) {NAME[st]}',fontsize=8); a.set_xlabel('Time since 03:10:00 UTC (s)')
    i=np.argmin(g.rr); info[st]=dict(t_spec=float(g.t[i]),rr=float(g.rr[i]))
ax[0].set_ylabel('Residual frequency (Hz)')
plt.tight_layout(); plt.savefig(OUT+'fig_bluebird.png',dpi=330); plt.savefig(OUT+'fig_bluebird.pdf')
print(info)
