# Figura: TEC vertical de GNSS (ROB/EPN) vs. perfil de Chapman ajustado ao resíduo Doppler da ISS
exec(open('figstyle.py').read())
R=json.load(open('../results/gnss_tec_cmp.json')); P=list(R)
lab=['3 Nov\n00:40','3 Nov\n02:20','3 Nov\n23:55','4 Nov\n01:30','4 Nov\n23:05']
x=np.arange(len(P))
fig,ax=plt.subplots(1,2,figsize=(7.2,2.6),gridspec_kw=dict(wspace=0.3))
a=ax[0]
V=np.array([R[p]['vtec_gnss'] for p in P]); lo=np.array([R[p]['vtec_range'][0] for p in P]); hi=np.array([R[p]['vtec_range'][1] for p in P])
a.errorbar(x,V,yerr=[V-lo,hi-V],fmt='s',color='#222',ms=5,capsize=2,label='GNSS VTEC (ROB/EPN maps)',zorder=4)
cols={'300':'#9ecae1','350':'#4292c6','420':'#08306b'}   # sequential: higher peak, darker
for i,hm in enumerate(['300','350','420']):
    a.plot(x+(i-1)*0.12,[R[p]['fits'][hm]['free']['vtec_total'] for p in P],'o',color=cols[hm],ms=4,label=f'fitted Chapman, total, $h_mF_2$={hm} km')
a.plot(x,[R[p]['fits']['300']['free']['vtec_below_iss'] for p in P],'_',color='#444',ms=12,mew=1.5,label='fitted, below the ISS (any $h_mF_2$)')
a.set_xticks(x); a.set_xticklabels(lab,fontsize=7); a.set_ylabel('Vertical TEC (TECU)'); a.set_ylim(0,24)
a.set_title('(a) Fitted profile vs. GNSS TEC',fontsize=8); a.legend(fontsize=7,loc='upper left',frameon=False); a.grid(alpha=0.25,lw=0.4)
b=ax[1]
for i,hm in enumerate(['300','350','420']):
    b.plot(x+(i-1)*0.12,[R[p]['fits'][hm]['gnss_fixed']['rms'] for p in P],'o',color=cols[hm],ms=4,label=f'$N_mF_2$ from GNSS, $h_mF_2$={hm} km')
b.plot(x,[R[p]['fits']['300']['free']['rms'] for p in P],'_',color='#444',ms=12,mew=1.5,label='$N_mF_2$ free')
b.set_xticks(x); b.set_xticklabels(lab,fontsize=7); b.set_ylabel('Residual rms (Hz)'); b.set_ylim(0,1.0)
b.set_title('(b) Doppler residual with the GNSS constraint',fontsize=8); b.legend(fontsize=7,loc='upper left',frameon=False); b.grid(alpha=0.25,lw=0.4)
plt.savefig(OUT+'fig_gnss_tec.pdf',bbox_inches='tight'); plt.savefig(OUT+'fig_gnss_tec.png',bbox_inches='tight',dpi=200)
