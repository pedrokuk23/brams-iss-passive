# Figure of the blind search: (a) ISS complex, detection probability per recording against the sigma_min of the station
# (Phase 2), with Phase 1; (b) detection probability per class with the satellite in the Earth's shadow or in sunlight,
# Phase 1 and Phase 2, with 95% intervals, the wider of the bootstraps over recordings and over satellites.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import beta
exec(open('figstyle.py').read()); plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False})
def cp(k,n,a=0.05): return (0.0 if k==0 else beta.ppf(a/2,k,n-k+1), 1.0 if k==n else beta.ppf(1-a/2,k+1,n-k))
C=json.load(open('../results/phase2_compare.json'))['ISS complex by station sigma_min']
S=json.load(open('../results/cluster_stats.json'))
fig,(a,b)=plt.subplots(1,2,figsize=(7.2,2.7),gridspec_kw=dict(width_ratios=[1,1.35]))
x=[np.sqrt(r['lo']*min(r['hi'],3e4)) for r in C]; p=[r['k']/r['n'] for r in C]; lo=[r['ci'][0] for r in C]; hi=[r['ci'][1] for r in C]
a.errorbar(x,p,yerr=[np.array(p)-lo,np.array(hi)-p],fmt='s',color='#c0392b',ms=4,capsize=2,lw=0.8,label='Phase 2 (419 recordings)')
for xi,r in zip(x,C): a.text(xi*1.15,r['k']/r['n']+0.02,f"{r['k']}/{r['n']}",ha='left',fontsize=7,color='#444')
v=S['classes']['phase1']['ISS complex']; a.errorbar([501],[v['pd']],yerr=[[v['pd']-v['ci_binomial'][0]],[v['ci_binomial'][1]-v['pd']]],fmt='o',color='#2a78d6',ms=4,capsize=2,lw=0.8,label='Phase 1 (22 recordings, median $\\sigma_\\mathrm{min}$ 500 m$^2$)')
for e in (100,300,1e3,3e3,1e4): a.axvline(e,color='#ddd',lw=0.5,zorder=0)
a.set_xscale('log'); a.set_xlim(80,3e4); a.set_ylim(-0.03,1.35); a.set_xlabel('$\\sigma_\\mathrm{min}$ of the ISS pass (m$^2$)'); a.set_ylabel('Detection probability')
a.set_yticks([0,0.2,0.4,0.6,0.8,1.0]); a.legend(fontsize=7,loc='upper right'); a.set_title('(a) ISS complex',fontsize=8,loc='left')
cl=['Starlink V2-mini','Starlink DTC','AST SpaceMobile','Starlink V1.x']; lab=['V2-mini','DTC','AST\nBlueBird','V1.x\n(control)']
for j,(ph,mk,dx) in enumerate((('phase1','o',-0.2),('phase2','s',0.2))):
    for g,col,ddx in (('plate_shadow','#222',-0.07),('plate_sunlit','#eda100',0.07)):
        for i,c in enumerate(cl):
            r=S['eclipse'][ph][f'{c}|{g}']; y=max(r['pd'],1e-3); l,h=r['ci_twoway']   # two-way bootstrap over satellites and time clusters
            b.errorbar(i+dx+ddx,y,yerr=[[y-max(l,1e-3)],[max(h-y,0)]],fmt=mk,color=col,mfc=col if g=='plate_shadow' else 'white',ms=4,capsize=1.5,lw=0.8,
                       label=(f"{'Phase 1' if ph=='phase1' else 'Phase 2'}, {'Earth shadow' if g=='plate_shadow' else 'sunlit'}" if i==0 else None))
b.set_yscale('log'); b.set_ylim(8e-4,4); b.set_xticks(range(4)); b.set_xticklabels(lab); b.set_xlim(-0.6,3.6)
b.set_ylabel('Detection probability per pass'); b.set_title('(b) Satellites with flat arrays',fontsize=8,loc='left')
b.legend(fontsize=7,loc='upper center',ncol=2)
fig.tight_layout(); fig.savefig(FIGDIR+'fig_blind.pdf'); fig.savefig('../results/fig_blind.png',dpi=200)
