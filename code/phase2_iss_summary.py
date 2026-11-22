# Summary of the Phase 2 ISS timing (Sec. 6, "Phase 2 passes"): reads results/phase2_iss_timing.json (written by
# phase2_iss_timing.py) and writes results/phase2_iss_summary.json with every number quoted in the text.
#  - recordings of the same pass (start times within 10 min, e.g. 25 Oct 06:10 and 06:15) are merged: per receiver the
#    fit with more points is kept;
#  - latency of the version-5 receivers: one offset per pass plus a common latency of Humain and Uccle (as in Phase 1)
#    and of Humain, Uccle and Maubeuge (used for the correction, since Maubeuge lags by the same amount);
#  - comparison with the NASA OEM: for each pass, the latest file created BEFORE the pass (a prediction); ages are given
#    both from the file's creation and from its START_TIME (initial state), the same convention in both phases;
#  - divergence of successive OEM files, per-file statistics, element sets issued after the pass, pass order in a night.
import os, json, numpy as np
D=json.load(open('../results/phase2_iss_timing.json'))
V5=('BEHUMA','BEUCCE'); V5M=('BEHUMA','BEUCCE','FRMAUB')
fits={k:v for k,v in D.items() if 'dt_ms' in v}
assert fits, 'no fits in phase2_iss_timing.json'
# --- group recordings into passes
starts=sorted({k[:13] for k in fits})
def tsec(s): import datetime as d; return d.datetime.strptime(s,'%Y%m%d_%H%M').timestamp()
groups=[]
for s in starts:
    if groups and tsec(s)-tsec(groups[-1][-1])<=600: groups[-1].append(s)
    else: groups.append([s])
passes={}
for g in groups:
    rx={}
    for s in g:
        for k,v in fits.items():
            if k[:13]==s and (k[-6:] not in rx or v['n']>rx[k[-6:]]['n']): rx[k[-6:]]=v
    passes[g[0]]=rx
# --- receiver offsets (rx_offsets.py): dt = a_pass + b[group], Oudsbergen as reference; version-5 group either Humain
#     and Uccle (as in Phase 1) or Humain, Uccle and Maubeuge (used for the correction); Ottignies has its own offset
from rx_offsets import fit_offsets
def latency(v5):
    M=fit_offsets([(p,r,v['dt_ms'],v.get('sdt_ms')) for p,rx in passes.items() for r,v in rx.items()],v5=v5)
    return dict(L_ms=M['offsets_ms']['v5'],L_se_ms=M['offsets_se_ms']['v5'],ottignies_ms=M['offsets_ms']['ottignies'],
                ottignies_se_ms=M['offsets_se_ms']['ottignies'],resid_sd_ms=M['resid_sd_ms'],resid_rms_ms=M['resid_rms_ms'],
                n_fits=M['n'],dof=M['dof'])
lat=latency(V5); latM=latency(V5M)
# per-receiver offset relative to Oudsbergen in the same pass
rel={}
for r in ('BEHUMA','BEUCCE','FRMAUB','BEOTTI','BEOUDS'):
    d=[]
    for p,rx in passes.items():
        ref=[rx[q]['dt_ms'] for q in ('BEOUDS',) if q in rx and q!=r]
        if r in rx and ref: d.append(rx[r]['dt_ms']-np.mean(ref))
    if d: rel[r]=dict(mean_ms=float(np.mean(d)),se_ms=float(np.std(d,ddof=1)/np.sqrt(len(d))) if len(d)>1 else None,n=len(d))
# correction applied below: the latency fitted with all three version-5 receivers (Maubeuge lags like Humain and Uccle)
L=latM['L_ms']; LO=latM['ottignies_ms']
def corr(r,v): return v['dt_ms']-(L if r in V5M else LO if r=='BEOTTI' else 0.0)
# --- per pass: mean corrected dt, OEM dt (mean over receivers) per file
OEMSTART={}
per=[]
for p,rx in sorted(passes.items()):
    dt=np.mean([corr(r,v) for r,v in rx.items()])
    labs=sorted(set().union(*[v['oem'].keys() for v in rx.values()]))   # sorted: byte-stable output
    oem={}
    for lab in labs:
        vv=[v['oem'][lab] for v in rx.values() if lab in v['oem']]
        oem[lab]=dict(dt_ms=float(np.mean([w['dt_ms'] for w in vv])),issued_before=vv[0]['issued_before'],
                      age_h=float(np.mean([w['age_h'] for w in vv])),prop_age_h=float(np.mean([w['prop_age_h'] for w in vv])))
    post=[v['tle_post']['0'] for v in rx.values() if '0' in v.get('tle_post',{})]
    per.append(dict(pass_=p,files=[s for g in groups if g[0]==p for s in g],receivers=sorted(rx),dt_ms=float(dt),
                    t_pass=float(np.mean([v['t_pass'] for v in rx.values()])),oem=oem,
                    tle_age_h=float(np.mean([v['tle_age_h'] for v in rx.values()])),
                    tle_post_ms=float(np.mean([w['dt_ms'] for w in post])) if post else None,
                    tle_post_lag_h=float(np.mean([w['epoch_after_h'] for w in post])) if post else None))
# pass order within a night (passes within 6 h of each other)
night=0; k=0
for i,q in enumerate(per):
    if i and q['t_pass']-per[i-1]['t_pass']>6*3600: night+=1; k=0
    elif i: k+=1
    q['night']=night; q['order']=k
# --- prediction comparison: latest OEM created before the pass
def stats(d):
    d=np.array(d); return dict(n=len(d),mean_ms=float(d.mean()),sd_ms=float(d.std(ddof=1)) if len(d)>1 else None,
                               rms_ms=float(np.sqrt(np.mean(d**2))),se_ms=float(d.std(ddof=1)/np.sqrt(len(d))) if len(d)>1 else None)
pred=[]
for q in per:
    c=[(v['age_h'],lab,v) for lab,v in q['oem'].items() if v['issued_before']]
    if not c: continue
    a,lab,v=min(c); pred.append(dict(pass_=q['pass_'],oem=lab,age_h=a,prop_age_h=v['prop_age_h'],dt_ms=q['dt_ms'],
                                     dt_oem_ms=v['dt_ms'],diff_ms=q['dt_ms']-v['dt_ms']))
bins=[(0,24),(24,48),(48,72),(72,1e9)]
by_age={'%d-%s h'%(a,'' if b>1e8 else b):stats([x['diff_ms'] for x in pred if a<=x['age_h']<b]) for a,b in bins
        if any(a<=x['age_h']<b for x in pred)}
by_prop={'%d-%s h'%(a,'' if b>1e8 else b):stats([x['diff_ms'] for x in pred if a<=x['prop_age_h']<b]) for a,b in bins
         if any(a<=x['prop_age_h']<b for x in pred)}
young=[x for x in pred if x['age_h']<24]
reg=None
if len(young)>2:
    X=np.array([x['dt_oem_ms'] for x in young]); Y=np.array([x['dt_ms'] for x in young])
    sl,ic=np.polyfit(X,Y,1); r=np.corrcoef(X,Y)[0,1]; res=Y-(sl*X+ic)
    se=np.sqrt(res@res/(len(X)-2)/np.sum((X-X.mean())**2)); reg=dict(slope=float(sl),slope_se=float(se),intercept_ms=float(ic),r=float(r))
fresh=[x for x in pred if x['prop_age_h']<24]
X=np.array([x['dt_oem_ms'] for x in fresh]); Y=np.array([x['dt_ms'] for x in fresh]); sl,ic=np.polyfit(X,Y,1); res=Y-(sl*X+ic)
reg_fresh=dict(n=len(fresh),slope=float(sl),slope_se=float(np.sqrt(res@res/(len(X)-2)/np.sum((X-X.mean())**2))),intercept_ms=float(ic),
               r=float(np.corrcoef(X,Y)[0,1]),stats=stats([x['diff_ms'] for x in fresh]))
pf={}
for x in fresh: pf.setdefault(x['oem'],[]).append(x['diff_ms'])
fm=np.array([np.mean(v) for v in pf.values()])
reg_fresh['per_file_mean_ms']={k:float(np.mean(v)) for k,v in pf.items()}
reg_fresh['n_files']=len(pf); reg_fresh['mean_of_file_means_ms']=float(fm.mean()); reg_fresh['se_file_level_ms']=float(fm.std(ddof=1)/np.sqrt(len(fm)))
excl={'%d-%s h'%(a,'' if b>1e8 else b):stats([x['diff_ms'] for x in pred if a<=x['prop_age_h']<b and x['oem']!='20261025T0255'])
      for a,b in bins if any(a<=x['prop_age_h']<b and x['oem']!='20261025T0255' for x in pred)}
per_file={}
for x in young: per_file.setdefault(x['oem'],[]).append(x['diff_ms'])
per_file={k:stats(v) for k,v in per_file.items()}
# --- divergence of successive OEM files at the passes they both cover
labs=sorted({l for q in per for l in q['oem']})
div=[]
for a,b in zip(labs[:-1],labs[1:]):
    for q in per:
        if a in q['oem'] and b in q['oem'] and q['oem'][b]['prop_age_h']>=0:
            div.append(dict(older=a,newer=b,pass_=q['pass_'],t_since_newer_start_h=q['oem'][b]['prop_age_h'],
                            diff_ms=q['oem'][a]['dt_ms']-q['oem'][b]['dt_ms']))
# --- element sets issued after the pass
tp=[q for q in per if q['tle_post_ms'] is not None]
post=dict(stats([q['dt_ms']-q['tle_post_ms'] for q in tp]),lag_h_min=float(min(q['tle_post_lag_h'] for q in tp)),
          lag_h_max=float(max(q['tle_post_lag_h'] for q in tp)))
# --- pass order within a night
order={}
for k in sorted({q['order'] for q in per}):
    qq=[q for q in per if q['order']==k]
    order[str(k)]=dict(n=len(qq),dt_mean_ms=float(np.mean([q['dt_ms'] for q in qq])),
                       tle_age_h_mean=float(np.mean([q['tle_age_h'] for q in qq])),
                       post_mean_ms=float(np.mean([q['dt_ms']-q['tle_post_ms'] for q in qq if q['tle_post_ms'] is not None])))
out=dict(n_recordings=len({k[:13] for k in fits}),n_passes=len(per),n_fits=len(fits),
         n_fits_used=sum(len(q['receivers']) for q in per),latency_v5_HU=lat,latency_v5_HUM=latM,receiver_offsets=rel,
         prediction=pred,by_creation_age=by_age,by_start_age=by_prop,young_regression=reg,fresh_start_lt24h=reg_fresh,by_start_age_excl_20261025=excl,young_per_file=per_file,
         oem_divergence=div,tle_post=post,pass_order=order,passes=per)
json.dump(out,open('../results/phase2_iss_summary.json','w'),indent=1)
print(json.dumps({k:out[k] for k in ('n_recordings','n_passes','n_fits','n_fits_used','latency_v5_HU','latency_v5_HUM',
      'receiver_offsets','by_creation_age','by_start_age','young_regression','fresh_start_lt24h','by_start_age_excl_20261025','tle_post','pass_order')},indent=1))
