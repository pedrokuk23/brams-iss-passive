# Figure: dechirped spectrograms of the close Starlink passes at Oudsbergen with z >= 7 (with line removal):
# six echoes that follow the prediction, one short marginal echo and one broadband burst (meteor) for comparison.
import os, json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('bb2.py').read())
from linefilt import maybe_clean
exec(open('glint.py').read().split("def glint_stat")[0].split("WS=[4,8,16,32,64]")[1])
D=os.environ.get('WAVDIR','../data/raw_wav/'); DF=np.arange(-20,20.01,0.5)
C={(r['obj'],r['file']):r for r in json.load(open('../results/close_check.json'))}
sel=[('STARLINK-11196 [DTC]','20261103_0040_BEOUDS','echo'),('STARLINK-11136 [DTC]','20261103_0045_BEOUDS','echo'),
     ('STARLINK-11330 [DTC]','20261103_0045_BEOUDS','echo'),('STARLINK-11317 [DTC]','20261104_0130_BEOUDS','echo'),
     ('STARLINK-11488 [DTC]','20261104_0125_BEOUDS','echo'),('STARLINK-11466 [DTC]','20261104_0130_BEOUDS','echo'),
     ('STARLINK-11293 [DTC]','20261104_0125_BEOUDS','short echo'),('STARLINK-11340 [DTC]','20261103_2350_BEOUDS','burst (meteor)')]
_G=dict(globals()); exec(open('figstyle.py').read()); plt.rcParams.update({'font.size':7})
fig,ax=plt.subplots(4,2,figsize=(6.4,7.6),sharey=True)
for a,(obj,f,kind),L in zip(ax.flat,sel,'abcdefgh'):
    r=C[(obj,f)]; p=D+f'RAD_BEDOUR_{f}_SYS001.wav'; md,x,_=bwav.read(p); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    g=GeoSat(p,obj); s,vb=snr_map(bb,ts,fs,fb,g,r['tau'],dfs=DF); t=(np.arange(len(s))+0.5)*0.5
    w=vb&(abs(t-r['t_echo'])<25)
    a.pcolormesh(t[w]-r['t_echo'],DF,10*np.log10(np.clip(s[w],1e-2,None)).T,vmin=0,vmax=15,cmap='Greys',shading='auto',rasterized=True)
    a.set_title(f"({L}) {obj.replace(' [DTC]','')}\n{int(f[6:8])} Nov {f[9:11]}:{f[11:13]} UTC, $z$={r['z']:.0f}, $\\tau$={r['tau']:+.1f} s".replace('-','\u2212').replace('STARLINK\u2212','STARLINK-'),fontsize=7)
    a.text(-23.5,17,kind,fontsize=7,color='#b03a2e' if kind!='echo' else '#1c5cab',bbox=dict(fc='white',ec='none',pad=0.5))
    a.set_xlim(-25,25)
for a in ax[-1]: a.set_xlabel('Time from echo peak (s)')
for a in ax[:,0]: a.set_ylabel('Residual frequency (Hz)')
fig.subplots_adjust(hspace=0.7,wspace=0.12)
cb=fig.colorbar(ax[0,0].collections[0],ax=ax,shrink=0.5,pad=0.02); cb.set_label('Block power above noise (dB)')
fig.savefig(FIGDIR+'fig_starlink.png',dpi=330,bbox_inches='tight')
