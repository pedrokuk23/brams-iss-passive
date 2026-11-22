# Avaliação: distribuição nula, probabilidade de detecção e erro do estimador de δt (ecos sintéticos)
import json, numpy as np, glob
exec(open('figstyle.py').read())
N=json.load(open('../results/null_iss.json')); zn=np.array([r['z'] for r in N])
S=[]; [S.extend(json.load(open(f))) for f in sorted(glob.glob('../results/synth_BEO*.json'))]
z=np.array([r['z'] for r in S]); al=np.array([r['alpha'] for r in S]); tau=np.array([r['tau'] for r in S])
e=np.array([(r['dt_est_ms']-r['dt_true_ms']) if r['dt_est_ms'] is not None and r['npts']>=5 else np.nan for r in S])
npts=np.array([r['npts'] for r in S])
q=lambda p: np.quantile(zn,p)
print(f'null: n={len(zn)} median {np.median(zn):.2f} q95 {q(.95):.2f} q99 {q(.99):.2f} max {zn.max():.2f}')
for th in [4,6,8,10]: print(f'  Pfa(z>={th}) = {np.mean(zn>=th):.3f}  ({np.sum(zn>=th)}/{len(zn)})')
TH=7.0  # acima de todos os máximos nulos com |tau|<=0.5 s
print('threshold:',TH, ' null(|tau|<=0.5) max', max(r['z'] for r in N if abs(r['tau'])<=0.5))
ok=(z>=TH)&(np.abs(tau)<=0.5)
res={'null':dict(n=int(len(zn)),median=float(np.median(zn)),q95=float(q(.95)),q99=float(q(.99)),max=float(zn.max()),pfa={str(t):float(np.mean(zn>=t)) for t in [4,6,8,10]}),'threshold':TH,'by_alpha':[], 'by_z':[]}
for a in sorted(set(al)):
    m=al==a; d=ok[m]; ee=e[m&ok&np.isfinite(e)]
    row=dict(alpha=float(a),n=int(m.sum()),Pd=float(d.mean()),z_med=float(np.median(z[m])),bias_ms=float(np.mean(ee)) if len(ee) else None,std_ms=float(np.std(ee,ddof=1)) if len(ee)>1 else None,rms_ms=float(np.sqrt(np.mean(ee**2))) if len(ee) else None)
    res['by_alpha'].append(row); print(row)
bins=[(TH,20),(20,40),(40,80),(80,1e9)]
for lo,hi in bins:
    m=ok&(z>=lo)&(z<hi)&np.isfinite(e); ee=e[m]
    if len(ee)>2:
        row=dict(z_lo=lo,z_hi=hi,n=int(m.sum()),bias_ms=float(np.mean(ee)),std_ms=float(np.std(ee,ddof=1)),rms_ms=float(np.sqrt(np.mean(ee**2))),p95_abs_ms=float(np.quantile(np.abs(ee),.95)))
        res['by_z'].append(row); print(row)
json.dump(res,open('../results/synthetic_summary.json','w'),indent=1)
# figura
T=json.load(open('../results/iss_detection_table.json'))
fig,ax=plt.subplots(1,3,figsize=(7.16,2.25),gridspec_kw=dict(width_ratios=[1.15,1,1.15]))
b=np.arange(0,62,0.5)
ax[0].hist(zn,bins=b,color='#b5b5b0',label='No target (265)')
for st in ['BEOTTI','BEOUDS','FRMAUB']:
    v=[r['z'] for r in T if r['st']==st and abs(r['tau'])<=0.5]
    ax[0].plot(v,[0.5]*len(v),MK[st],color=COL[st],ms=4,label='ISS, '+NAME[st],ls='none')
ax[0].set_yscale('log'); ax[0].set_ylim(0.3,200); ax[0].set_xlabel(r'Maximum $z$ over the grid'); ax[0].set_ylabel('Searches')
ax[0].axvline(7,color='#c0392b',ls='--',lw=0.8,label='threshold $z=7$'); ax[0].legend(fontsize=7,loc='upper right'); ax[0].set_title('(a) Null and ISS',fontsize=8)
A=[r['alpha'] for r in res['by_alpha']]; P=[r['Pd'] for r in res['by_alpha']]
ax[1].plot(A,P,'o-',color='#2a78d6',ms=3,lw=1)
for r in res['by_alpha']: ax[1].annotate(f"{r['z_med']:.0f}",(r['alpha'],r['Pd']),textcoords='offset points',xytext=((-10,4) if r['alpha']==0.03 else (2,6)) if r['Pd']<0.1 else (4,-9),fontsize=7,color='#555')
ax[1].text(0.98,0.42,'labels: median $z$',transform=ax[1].transAxes,ha='right',fontsize=7,color='#555')
ax[1].set_xlabel(r'Injected amplitude / noise rms'); ax[1].set_ylabel(r'$P_d$ ($z\geq7$, $|\tau|\leq0.5$ s)'); ax[1].set_ylim(-0.03,1.05)
ax[1].set_title('(b) Detection probability',fontsize=8)
m=ok&np.isfinite(e)
ax[2].scatter(z[m],e[m],s=4,color='#2a78d6',alpha=0.6,lw=0); ax[2].axhline(0,color='#555',lw=0.6); ax[2].set_xscale('log')
ax[2].set_ylim(-40,40); ax[2].set_xlabel(r'Detection statistic $z$'); ax[2].set_ylabel(r'$\hat{\delta t}-\delta t$ (ms)'); ax[2].set_title(r'(c) Error of $\hat{\delta t}$',fontsize=8)
for a_ in ax:
    for s_ in ['top','right']: a_.spines[s_].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'fig_synthetic.pdf'); plt.savefig(OUT+'fig_synthetic.png')
