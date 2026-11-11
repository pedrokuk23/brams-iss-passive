# Graphical abstract (Elsevier: about 13 x 5 cm displayed, legible fonts): (a) blind search of the validation campaign,
# detection fraction per class against the pooled false-alarm rate; (b) ISS timing against NASA predictions propagated
# less than one day (both campaigns) and older ones.
import numpy as np, json
exec(open('figstyle.py').read())
P2=json.load(open('../results/phase2_summary.json')); C=json.load(open('../results/iss_oem_combined.json'))
fig,ax=plt.subplots(1,2,figsize=(5.4,2.25),gridspec_kw=dict(width_ratios=[1.25,1],wspace=0.55))
b=ax[0]
cls=[('ISS complex','ISS'),('AST SpaceMobile','BlueBird'),('Starlink DTC','Starlink DTC'),('Starlink V2-mini','V2-mini'),('Starlink V1.x','V1.x')]
for i,(c,lab) in enumerate(cls):
    v=P2['classes'][c]['all']; lo,hi=v['ci']
    b.errorbar(100*v['pd'],i,xerr=[[100*(v['pd']-lo)],[100*(hi-v['pd'])]],fmt='o',color='#2b5d8a',ms=4,lw=1,capsize=2)
b.axvspan(0.12,0.16,color='#c0392b',alpha=0.4,lw=0)
b.annotate('false alarms',(0.18,len(cls)-0.75),ha='left',va='center',fontsize=7,color='#c0392b')
b.set_xscale('log'); b.set_xlim(0.08,100); b.set_ylim(len(cls)-0.4+0.3,-0.6)
b.set_yticks(range(len(cls))); b.set_yticklabels([l for _,l in cls],fontsize=8); b.tick_params(labelsize=7.5)
b.set_xlabel('Passes detected (%)',fontsize=8); b.set_title('Blind search, 24,284 passes',fontsize=8.5)
for s_ in ['top','right']: b.spines[s_].set_visible(False)
c=ax[1]; R=[r for r in C['rows'] if r['oem']!='20261025T0255']
old=[r for r in R if r['age_start_h']>=24]; fr=[r for r in R if r['age_start_h']<24]
c.plot([r['dt_oem_ms'] for r in old],[r['dt_ms'] for r in old],'o',mfc='none',mec='#aaa',ms=3.5,mew=0.7,label='> 1 day')
c.plot([r['dt_oem_ms'] for r in fr],[r['dt_ms'] for r in fr],'o',color='#2b5d8a',ms=4,label='< 1 day')
c.plot([-400,400],[-400,400],':',color='#888',lw=0.6); c.set_xlim(-550,550); c.set_ylim(-300,300)
c.set_xticks([-400,0,400]); c.set_yticks([-200,0,200]); c.tick_params(labelsize=7.5)
c.set_xlabel('NASA prediction (ms)',fontsize=8); c.set_ylabel('Measured $\\delta t$ (ms)',fontsize=8)
s=C['by_start_age_excl_20261025']['0-24 h|both']
c.set_title('ISS timing: %.0f ms rms\n(%d passes, prediction < 1 day)'%(s['rms_ms'],s['n']),fontsize=8.5)
c.legend(fontsize=7,loc='upper left',handletextpad=0.2,borderpad=0.2)
for s_ in ['top','right']: c.spines[s_].set_visible(False)
plt.savefig(OUT+'graphical_abstract.png',dpi=300,bbox_inches='tight',pad_inches=0.05)
