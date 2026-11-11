# Figura: padrão relativo das antenas receptoras obtido com aviões (modelo censurado, dia 04/11)
exec(open('figstyle.py').read())
F=json.load(open('../results/antenna_map_tobit.json')); STS=['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']
fig,ax=plt.subplots(1,4,figsize=(7.4,2.3),gridspec_kw=dict(wspace=0.5))
azc=np.arange(15,360,30)
for i,(k,lab) in enumerate([('low','el < 10°'),('high','el ≥ 10°')]):
    a=ax[i]
    for s in STS:
        y=np.array(F['az'][s][k])+F['off'][s]; n=np.array(F['az'][s]['ndet_'+k]); ok=n>=50
        a.plot(azc,y,'-',color=COL[s],lw=1.1)
        a.plot(azc[ok],y[ok],MK[s],color=COL[s],ms=3.2,label=NAME[s]); a.plot(azc[~ok],y[~ok],MK[s],mfc='white',mec=COL[s],ms=3.2)
    a.set_xlim(0,360); a.set_xticks([0,90,180,270,360]); a.set_xlabel('Azimuth at receiver (°)'); a.set_title(f'({"ab"[i]}) Azimuth, {lab}',fontsize=8)
    a.grid(alpha=0.25,lw=0.4); a.set_ylim(-16,16)
ax[0].set_ylabel('Relative gain (dB)')
ec=0.5*(np.array(F['el_edges'][:-1])+np.array(F['el_edges'][1:]))
for s in STS: ax[2].plot(ec,F['el'][s]['x'],'-'+MK[s],color=COL[s],ms=3,lw=1.1)
ax[2].set_xlabel('Elevation at receiver (°)'); ax[2].set_title('(c) Elevation',fontsize=8); ax[2].grid(alpha=0.25,lw=0.4)
bc=np.arange(5,180,10); ax[3].plot(bc,F['beta'],'-o',color='#444',ms=2.8,lw=1.1,label='bistatic angle')
tc=np.array([2.5,7.5,12.5,17.5,25,37.5,67.5]); ax[3].plot(tc,F['eltx'],'--s',color='#8a8a8a',ms=2.8,lw=1.1,label='elevation from TX')
ax[3].set_xlabel('Angle (°)'); ax[3].set_title('(d) Common terms',fontsize=8); ax[3].set_ylim(-6,15); ax[3].legend(fontsize=7,loc='upper left',handlelength=1.5); ax[3].grid(alpha=0.25,lw=0.4)
h,l=ax[0].get_legend_handles_labels(); fig.legend(h,l,loc='lower center',ncol=5,fontsize=7,bbox_to_anchor=(0.5,-0.2))
plt.savefig(OUT+'fig_antenna_map.pdf',bbox_inches='tight'); plt.savefig(OUT+'fig_antenna_map.png',bbox_inches='tight',dpi=200)
