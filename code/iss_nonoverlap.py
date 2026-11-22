# ISS timing fit with non-overlapping analysis windows (Sec. 6). The tracks of iono_tracks_v2.json come from 2-s
# windows every 0.5 s, so neighbouring points share data. Taking every fourth point gives four disjoint sets of
# non-overlapping 2-s windows per receiver (offsets 0, 0.5, 1 and 1.5 s). Each set is fitted with the model of
# iono_fit_v2.py (timing plus Chapman term, time-base correction applied). The spread of delta t over the four
# disjoint sets is compared with the formal error of one set, and the latency model of iss_checks2.py is refitted
# with each set. Writes results/iss_nonoverlap.json.
import json, numpy as np
src=open('iono_fit_v2.py').read().split('passes=sorted(')[0]
exec(src)
keys=sorted(k for k in tr if len(k)>13 and k[14:] in ('BEOTTI','BEOUDS','BEHUMA','BEUCCE') and len(cleaned(tr[k]))>=NMIN)
V5={'BEHUMA','BEUCCE'}
rows=[]
for k in keys:
    obs,mt,mi=cols(k); n=len(obs); sub=[]
    for ph in range(4):
        idx=np.arange(ph,n,4)
        A=np.column_stack([mt[idx],mi[idx]]); p,*_=np.linalg.lstsq(A,obs[idx],rcond=None); r=obs[idx]-A@p
        s2=r@r/max(len(idx)-2,1); C=np.linalg.inv(A.T@A)*s2
        sub.append((p[0]*1e3,np.sqrt(C[0,0])*1e3,len(idx)))
    d=np.array([x[0] for x in sub]); s=np.array([x[1] for x in sub])
    rows.append(dict(key=k,n_full=n,dt_sets_ms=d.tolist(),sdt_sets_ms=s.tolist(),spread_sd_ms=float(d.std(ddof=1)),
                     formal_one_set_ms=float(np.median(s)),dt_mean_ms=float(d.mean()),dt_full_ms=float(F_full[k]) if False else None))
full=json.load(open('../results/iono_fit_v2.json'))
for r in rows: r['dt_full_ms']=full[r['key']]['dt_ms']; r['sdt_full_ms']=full[r['key']]['sdt_ms']
ratio=np.array([r['spread_sd_ms']/r['formal_one_set_ms'] for r in rows])
# latency model per set
def latency(ph):
    passes=sorted({r['key'][:13] for r in rows}); X=np.zeros((len(rows),len(passes)+1)); y=np.array([r['dt_sets_ms'][ph] for r in rows]); s=np.array([r['sdt_sets_ms'][ph] for r in rows])
    for i,r in enumerate(rows): X[i,passes.index(r['key'][:13])]=1; X[i,-1]=1 if r['key'][14:] in V5 else 0
    b=np.linalg.lstsq(X,y,rcond=None)[0]; e=y-X@b; dof=len(y)-X.shape[1]
    return dict(latency_ms=float(b[-1]),resid_sd_ms=float(np.sqrt(e@e/dof)),chi2_red=float(np.sum((e/s)**2)/dof),dof=dof)
L=[latency(ph) for ph in range(4)]
out=dict(n_receiver_passes=len(rows),spread_over_formal_median=float(np.median(ratio)),spread_over_formal_range=[float(ratio.min()),float(ratio.max())],
         formal_one_set_ms=[float(min(r['formal_one_set_ms'] for r in rows)),float(max(r['formal_one_set_ms'] for r in rows))],
         max_abs_mean_minus_full_ms=float(max(abs(r['dt_mean_ms']-r['dt_full_ms']) for r in rows)),
         latency_by_set=L,chi2_red_median=float(np.median([x['chi2_red'] for x in L])),latency_range_ms=[min(x['latency_ms'] for x in L),max(x['latency_ms'] for x in L)],rows=rows)
json.dump(out,open('../results/iss_nonoverlap.json','w'),indent=1)
print({k:v for k,v in out.items() if k!='rows'})
for r in rows: print(r['key'],round(r['dt_full_ms'],1),np.round(r['dt_sets_ms'],1),round(r['spread_sd_ms'],1),round(r['formal_one_set_ms'],1),round(r['sdt_full_ms'],1))
