# Figures of the feasibility test: (1) decision flow; (2) detection probability against the margin M = sigma_max/sigma_min,
# Phase 1 (fit) and Phase 2 (prediction test), by group.
import json, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy.stats import beta
exec(open('figstyle.py').read()); plt.rcParams.update({'font.size':8,'axes.spines.top':False,'axes.spines.right':False})
def cp(k,n,a=0.05): return (0.0 if k==0 else beta.ppf(a/2,k,n-k+1), 1.0 if k==n else beta.ppf(1-a/2,k+1,n-k))
# ---- (1) decision flow
fig,ax=plt.subplots(figsize=(6.3,3.1)); ax.axis('off'); ax.set_xlim(-1,101); ax.set_ylim(-4,44)
W=18.4; X=[0,20.4,40.8,61.2,81.6]
B=['1. Geometry\n$30<|f_D|<900$\nHz for 20 s?',
   '2. Sensitivity\n$\\sigma_\\min$, Eq. (10)\n(range, gains,\n$z_\\mathrm{th}$, $T_\\mathrm{eff}$)',
   '3. Target\nplate: $4\\pi A^2/\\lambda^2$\nglint: $T_g/2$\ntumbling: $\\lambda P/4\\pi L$',
   '4. Margin\n$M=\\sigma_\\max/\\sigma_\\min$',
   '5. Group\narray in shadow\narray in sunlight\nlarge structure']
for x,t in zip(X,B):
    ax.add_patch(FancyBboxPatch((x,27),W,14,boxstyle='round,pad=0.3',fc='#eeeeee',ec='k',lw=0.6)); ax.text(x+W/2,34,t,ha='center',va='center',fontsize=7)
for x in X[:-1]: ax.annotate('',xy=(x+W+2.0,34),xytext=(x+W+0.4,34),arrowprops=dict(arrowstyle='->',lw=0.7))
ax.text(X[0]+W/2,24.3,'no: not observable',ha='center',fontsize=7,color='0.3')
# box labels from the fitted curves (feasibility_curve.json) and the M<1 counts (feasibility_validation.json)
_C=json.load(open('../results/feasibility_curve.json'))['groups']; _I=json.load(open('../results/feasibility_validation.json'))['impossible']
def _P(g,M): a,b=_C[g]['a'],_C[g]['b']; return 1/(1+np.exp(-(a+b*np.log10(M))))
def _pc(v): return f"{100*v:.2g}%" if v<0.1 else f"{100*v:.0f}%"
_i1,_i2=(100*_I[p]['k']/_I[p]['n'] for p in ('phase1','phase2'))
D=[(0.5,f'$M<1$: nominally not visible in\nany attitude (detected in {_i1:.2f}%\nand {_i2:.2f}% of passes)','#f4cccc'),
   (34.5,f'$1\\leq M<10$: marginal\nPhase 1 fit: arrays {_pc(min(_P("plate_sunlit",1),_P("plate_shadow",1)))}–{_pc(_P("plate_shadow",10))},\nISS {_pc(_P("structure",1))}–{_pc(_P("structure",10))}','#fff2cc'),
   (68.5,f'$M\\geq10$, Phase 1 fit:\nshadow {_pc(_P("plate_shadow",10))}–{_pc(_P("plate_shadow",100))}, sunlit {_pc(_P("plate_sunlit",10))}–{_pc(_P("plate_sunlit",100))},\nISS >{_pc(_P("structure",10))[:-1]}% (Phase 2 lower)','#d9ead3')]
for x,t,c in D:
    ax.add_patch(FancyBboxPatch((x,6),31,11,boxstyle='round,pad=0.3',fc=c,ec='k',lw=0.6)); ax.text(x+15.5,11.5,t,ha='center',va='center',fontsize=7)
ax.plot([81.6+W/2]*2,[26.6,21],color='k',lw=0.6); ax.plot([16,81.6+W/2],[21,21],color='k',lw=0.6)
for x in (16,50,84): ax.annotate('',xy=(x,17.6),xytext=(x,21),arrowprops=dict(arrowstyle='->',lw=0.7))
ax.text(50,1.5,'6. Expected detections = sum of $P_d(M)$ over the passes of a campaign; worth it if at least a few',ha='center',fontsize=7,style='italic')
fig.savefig(FIGDIR+'fig_feasibility_flow.pdf',bbox_inches='tight'); fig.savefig('../results/fig_feasibility_flow.png',dpi=200,bbox_inches='tight')
# ---- (2) validation
C=json.load(open('../results/feasibility_curve.json')); V=json.load(open('../results/feasibility_validation.json'))
E=[0,1,3,10,30,100,1000]; ctr=[0.6,np.sqrt(3),np.sqrt(30),np.sqrt(300),np.sqrt(3000),300]
fig,axs=plt.subplots(1,3,figsize=(7.2,2.4),sharey=True)
names={'structure':'(a) ISS complex','plate_shadow':'(b) Flat arrays, Earth shadow','plate_sunlit':'(c) Flat arrays, sunlit'}
for ax,g in zip(axs,('structure','plate_shadow','plate_sunlit')):
    gg=C['groups'][g]; xm=np.logspace(np.log10(0.3),3,400)          # logistic in log10 M for M >= 1, false-alarm rate below
    ym=np.where(xm<1,1e-3,1/(1+np.exp(-(gg['a']+gg['b']*np.log10(np.maximum(xm,1))))))
    ax.plot(xm,np.maximum(ym,1e-3),color='0.5',lw=1,label='logistic curve (fit on Phase 1)')
    for rows,mk,col,lab,dx in ((C['groups'][g]['phase1'],'o','#2a78d6','Phase 1',0.9),(V['by_bin'][g],'s','#c0392b','Phase 2',1.1)):
        for i,r in enumerate(rows):
            n=r['n']; k=r.get('k',r.get('observed'))
            if n<5: continue
            lo,hi=cp(k,n); p=max(k/n,1e-3); x=ctr[i]*dx
            ax.errorbar(x,p,yerr=[[p-max(lo,1e-3)],[max(hi-p,0)]],fmt=mk,ms=3.5,color=col,capsize=1.5,lw=0.7,label=lab if i==3 else None)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_ylim(8e-4,1.5); ax.set_xlim(0.3,1000); ax.axvline(1,color='k',ls=':',lw=0.7)
    ax.set_title(names[g],fontsize=8); ax.set_xlabel('margin $M=\\sigma_\\max/\\sigma_\\min$')
axs[0].set_ylabel('Detection probability'); axs[0].legend(fontsize=7,frameon=False,loc='lower right')
fig.tight_layout(); fig.savefig(FIGDIR+'fig_feasibility_valid.pdf'); fig.savefig('../results/fig_feasibility_valid.png',dpi=200)
