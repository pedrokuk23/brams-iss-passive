# Figure: (a) detection statistic versus the value predicted from the measured echo energy; (b) minimum detectable cross
# section versus geometric mean range, with the measured objects.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('detect_theory.py').read().split("cn=json.load")[0])
T=json.load(open('../results/detect_theory.json')); cn=json.load(open('../results/detectability_cn0.json'))
G={k:10**(v/10) for k,v in T['Geff_dB'].items()}
_G=dict(globals()); exec(open('figstyle.py').read())          # common style: serif font, receiver colors COL, markers MK, NAME
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False})
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,3.2),gridspec_kw=dict(width_ratios=[1,1.45]))
for st in COL:
    pts=[(r['z_pred'],r['z_obs']) for r in cn if r['label'].split()[-1]==st]
    if not pts: continue
    x,y=np.array(pts).T; a.loglog(x,y,MK[st],ms=4,color=COL[st],label=NAME[st],mec='white',mew=0.5,ls='')
bbp=[(r['z_pred'],r['z_obs']) for r in cn if r['label'].startswith('SPACE')]
x,y=np.array(bbp).T; a.loglog(x,y,'*',ms=9,color='none',mec='#222',mew=0.8,label='SPACEMOBILE-008')
g=np.array([3,800]); a.loglog(g,g,color='#999',lw=1); a.axhline(7,color='#999',lw=0.8,ls=':')
a.text(200,7.6,'$z_\\mathrm{th}=7$',color='#666'); a.set_xlabel('$(E/N_0)/\\sqrt{n_b}$ from the measured echo'); a.set_ylabel('$z$ of the search')
a.set_xlim(3,800); a.set_ylim(3,800); a.legend(frameon=False,fontsize=7,loc='upper left'); a.set_title('(a) statistic versus echo energy',loc='left',fontsize=8)
Rk=np.linspace(200,1000,200); 
for st,ls in [('BEOUDS','-'),('BEOTTI','-'),('FRMAUB','-')]:
    b.semilogy(Rk,[T['EN0_min_window']/50/cn0(1.,r*1e3,G[st]) for r in Rk],ls,color=COL[st],lw=1.5,label=f'{NAME[st]}, 50-s echo')
for Tg,ls in [(16,'--'),(2,':')]:
    b.semilogy(Rk,[T['EN0_min_window']/Tg/cn0(1.,r*1e3,G['BEOUDS']) for r in Rk],ls,color=COL['BEOUDS'],lw=1.2,label=f'Oudsbergen, glint, $T_\\mathrm{{eff}}$ = {Tg} s')
for r in T['iss_sigma']:
    st=r['label'].split()[-1]; rr=[c['rr_at_peak'] for c in cn if c['label']==r['label']][0]
    b.semilogy(rr,r['sigma_m2'],MK[st],ms=4,color=COL[st],mec='white',mew=0.5)
ce=json.load(open('../results/close_echo.json'))['echoes']
clear={'STARLINK-11196 [DTC]','STARLINK-11136 [DTC]','STARLINK-11330 [DTC]','STARLINK-11317 [DTC]','STARLINK-11488 [DTC]','STARLINK-11466 [DTC]'}
for e in ce:
    if e['obj'] in clear: b.semilogy(e['rr_km'],e['sigma_m2'],'^',ms=4.5,color=COL['BEOUDS'],mec='white',mew=0.5)
b.text(300,3.2e3,'Starlink\nDTC',fontsize=7,ha='center')
b.semilogy(538,T['sigma_bluebird_m2'],'*',ms=10,color='none',mec='#222',mew=0.8)
ul=T['envisat_upper_limits']['BEOUDS']
b.annotate('',xy=(ul['R_km'],ul['sigma_ul_window']*0.45),xytext=(ul['R_km'],ul['sigma_ul_window']),arrowprops=dict(arrowstyle='->',color='#222',lw=1))
b.plot(ul['R_km'],ul['sigma_ul_window'],'v',color='#222',ms=5); b.text(ul['R_km']-10,ul['sigma_ul_window']*0.5,'Envisat\n(nominal limit)',fontsize=7,va='top',ha='right')

b.text(445,1.6e5,'ISS',fontsize=7); b.text(548,5.5e4,'BlueBird\narray',fontsize=7)
for k,lab in [('resonant_wire','resonant $\\lambda/2$ wire'),('plate_5m','5-m flat plate'),('plate_15m','15-m flat plate'),('rayleigh_sphere_1m','sphere, 1-m diameter')]:
    v=T['ref_sigma_m2'][k]; b.axhline(v,color='#bbb',lw=0.7,ls=(0,(2,2)))
    if k=='plate_15m': b.text(205,v*1.15,lab,fontsize=7,color='#666',ha='left')      # left end: below the curves only there
    elif k=='resonant_wire': b.text(330,v*0.55,lab,fontsize=7,color='#666',ha='left')   # clear of the curves and of the legend
    else: b.text(995,v*1.15,lab,fontsize=7,color='#666',ha='right')                  # right end: above the curves only there
b.set_xlim(200,1000); b.set_ylim(0.3,5e5); b.set_xlabel('Geometric mean range $\\sqrt{R_1R_2}$ (km)'); b.set_ylabel('Cross section $\\sigma$ (m$^2$)')
b.legend(frameon=False,fontsize=7,loc='lower right',bbox_to_anchor=(1,0.08)); b.set_title('(b) minimum detectable cross section ($z=7$)',loc='left',fontsize=8)
fig.tight_layout(); fig.savefig(FIGDIR+'fig_detectability.pdf'); fig.savefig('../results/fig_detectability.png',dpi=200)
