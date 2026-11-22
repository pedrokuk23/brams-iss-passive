# Graphical abstract (Elsevier: at least 1328 x 531 pixels, read at about 13 x 5 cm): smallest cross section that the
# network detects against range, for a steady echo and for glints at Oudsbergen, with the objects detected (ISS,
# Starlink DTC, BlueBird array) and the upper limit for Envisat; based on Fig. 8b of the article.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('detect_theory.py').read().split("cn=json.load")[0])
T=json.load(open('../results/detect_theory.json')); cn=json.load(open('../results/detectability_cn0.json'))
G={k:10**(v/10) for k,v in T['Geff_dB'].items()}
exec(open('figstyle.py').read())
plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':9})
fig,b=plt.subplots(figsize=(8.0,3.2))
Rk=np.linspace(250,1000,200); c=COL['BEOUDS']
b.fill_between(Rk,[T['EN0_min_window']/50/cn0(1.,r*1e3,G['BEOUDS']) for r in Rk],5e5,color=c,alpha=0.07,lw=0)
b.semilogy(Rk,[T['EN0_min_window']/50/cn0(1.,r*1e3,G['BEOUDS']) for r in Rk],'-',color=c,lw=2,label='Oudsbergen limit, steady 50-s echo')
for Tg,ls in [(16,'--'),(2,':')]:
    b.semilogy(Rk,[T['EN0_min_window']/Tg/cn0(1.,r*1e3,G['BEOUDS']) for r in Rk],ls,color=c,lw=1.5,label=f'Oudsbergen limit, {Tg}-s glint')
iss=[(([x['rr_at_peak'] for x in cn if x['label']==r['label']][0]),r['sigma_m2']) for r in T['iss_sigma']]
x,y=np.array(iss).T; b.semilogy(x,y,'o',ms=5,color='#c0392b',mec='white',mew=0.5,ls='',label='ISS')
ce=json.load(open('../results/close_echo.json'))['echoes']
clear={'STARLINK-11196 [DTC]','STARLINK-11136 [DTC]','STARLINK-11330 [DTC]','STARLINK-11317 [DTC]','STARLINK-11488 [DTC]','STARLINK-11466 [DTC]'}
d=np.array([(e['rr_km'],e['sigma_m2']) for e in ce if e['obj'] in clear]).T
b.semilogy(d[0],d[1],'^',ms=6,color='#e67e22',mec='white',mew=0.5,ls='',label='Starlink DTC')
b.semilogy(538,T['sigma_bluebird_m2'],'*',ms=13,color='#2c3e50',ls='',label='BlueBird array')
ul=T['envisat_upper_limits']['BEOUDS']
b.annotate('',xy=(ul['R_km'],ul['sigma_ul_window']*0.4),xytext=(ul['R_km'],ul['sigma_ul_window']),arrowprops=dict(arrowstyle='->',color='#555',lw=1.2))
b.plot(ul['R_km'],ul['sigma_ul_window'],'v',color='#555',ms=6); b.text(ul['R_km']-12,ul['sigma_ul_window']*0.22,'Envisat, tumbling:\nnot detected',fontsize=8,va='center',ha='right',color='#555')
for k,lab in [('plate_15m','15-m flat plate'),('plate_5m','5-m flat plate'),('resonant_wire','resonant $\\lambda/2$ wire')]:
    v=T['ref_sigma_m2'][k]; b.axhline(v,color='#bbb',lw=0.8,ls=(0,(2,2))); b.text(990,v*1.2,lab,fontsize=8,color='#777',ha='right')
b.set_xlim(250,1000); b.set_ylim(3,5e5); b.set_xlabel('Geometric mean range (km)'); b.set_ylabel('Radar cross section (m$^2$)')
b.set_title('An open meteor radar network (BRAMS, 50 MHz) as a multistatic sensor of objects in low Earth orbit',loc='left',fontsize=9.5)
b.legend(frameon=False,fontsize=8,loc='upper left',bbox_to_anchor=(1.0,1.0))
fig.tight_layout(); fig.savefig(FIGDIR+'graphical_abstract.png',dpi=332,facecolor='white'); from PIL import Image; Image.open(FIGDIR+'graphical_abstract.png').convert('RGB').save(FIGDIR+'graphical_abstract.png'); fig.savefig(FIGDIR+'graphical_abstract.pdf')
