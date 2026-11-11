# Figura: linhas de interferência estreitas no sinal em banda base e efeito da sua remoção na detecção da ISS
exec(open('synth.py').read().split("if __name__")[0])
exec(open('figstyle.py').read())
from linefilt import clean
from scipy import signal as sg
fig,ax=plt.subplots(1,2,figsize=(7.2,2.5),gridspec_kw=dict(width_ratios=[1.6,1],wspace=0.28))
a=ax[0]
for st,off in [('BEHUMA',0),('BEOUDS',-60)]:
    md,x,_=bwav.read(D+f'RAD_BEDOUR_20261104_2305_{st}_SYS001.wav'); fs=md['samplerate']; fb,bb,ts=baseband(x,fs)
    for y,lab in [(clean(bb,ts,fs),'after'),(bb,'before')]:
        f,P=sg.welch(y,fs,nperseg=8192,return_onesided=False); f=np.fft.fftshift(f); L=10*np.log10(np.fft.fftshift(P)); L-=np.median(L)
        if lab=='before': a.plot(f,L+off,'-',lw=0.7,color=COL[st])
        else: a.plot(f,L+off,'-',lw=0.9,color='#222')
a.axvspan(-900,-30,color='#999',alpha=0.12,lw=0); a.axvspan(30,900,color='#999',alpha=0.12,lw=0)
a.text(-880,92,'band searched for the ISS',fontsize=7,color='#555')
a.set_xlim(-1500,1500); a.set_ylim(-75,100); a.set_xlabel('Frequency relative to the beacon (Hz)'); a.set_ylabel('Power spectral density (dB, offset)')
a.annotate('+500 Hz',xy=(500,85),xytext=(650,70),fontsize=7,arrowprops=dict(arrowstyle='-',lw=0.5,color='#555'))
a.annotate('+600 Hz',xy=(600,86),xytext=(950,88),fontsize=7,arrowprops=dict(arrowstyle='-',lw=0.5,color='#555'))
a.text(-1450,-3,'Humain',fontsize=7,color=COL['BEHUMA']); a.text(-1450,-63,'Oudsbergen (−60 dB)',fontsize=7,color=COL['BEOUDS'])
a.set_title('(a) Baseband spectrum (color: before, black: after line removal)',fontsize=7.5)
b=ax[1]; T=json.load(open('../results/iss_detection_table_v2.json'))
for r in T:
    st=r['st']; b.plot(max(r['z_raw'],1),r['z'],MK[st],color=COL[st],ms=4,ls='',label=NAME[st])
h,l=b.get_legend_handles_labels(); u=dict(zip(l,h)); b.legend(u.values(),u.keys(),fontsize=7,loc='lower right')
b.set_xscale('log'); b.set_yscale('log'); b.plot([1,400],[1,400],':',color='#888',lw=0.6)
b.axhline(7,ls='--',color='#888',lw=0.6); b.axvline(7,ls='--',color='#888',lw=0.6)
b.set_xlim(1,400); b.set_ylim(1,400); b.set_xlabel('$z$ without line removal'); b.set_ylabel('$z$ with line removal')
b.set_title('(b) ISS detection statistic',fontsize=7.5)
plt.savefig(OUT+'fig_interference.pdf',bbox_inches='tight'); plt.savefig(OUT+'fig_interference.png',bbox_inches='tight',dpi=200)
