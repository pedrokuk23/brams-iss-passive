# ISS timing against the NASA OEM predictions, both campaigns with one rule (Sec. 6 and Supplementary S10).
# For every pass: the Doppler-equivalent offset delta t (mean over the receivers of the pass, version-5 receivers
# corrected for their fitted latency: Phase 1 Humain and Uccle, iss_checks2.json; Phase 2 Humain, Uccle and
# Maubeuge, phase2_iss_summary.json) and delta t_OEM of the latest OEM file created before the pass (mean over the same
# receivers). Two ages of that prediction are given: from its initial state (START_TIME, the propagation span) and
# from its creation. Reads iono_fit_v2.json, oem_cmp_v2.json, iss_checks2.json, phase2_iss_summary.json;
# writes results/iss_oem_combined.json.
import json, datetime as dt, numpy as np
R='../results/'
def ts(s): return dt.datetime.strptime(s,'%Y-%m-%dT%H:%M').replace(tzinfo=dt.timezone.utc).timestamp()
F=json.load(open(R+'iono_fit_v2.json')); O=json.load(open(R+'oem_cmp_v2.json')); D=json.load(open(R+'iss_checks2.json'))['latency_model']
Q=json.load(open(R+'phase2_iss_summary.json'))
V5_1=('BEHUMA','BEUCCE'); L1=D['latency_ms']   # joint model of Sec. 6 (iss_checks2.py)
# Phase 1 OEM files (data/oem_archive): creation and initial state
P1FILES={'1102':dict(created=ts('2026-11-02T15:17'),start=ts('2026-11-02T12:00'),name='20261102T1517'),
         '1104':dict(created=ts('2026-11-04T16:10'),start=ts('2026-11-04T12:00'),name='20261104T1610')}
rows=[]
p1=sorted(k for k in F if len(k)==13)
for k in p1:
    rx=F[k]['receivers']; t=dt.datetime.strptime(k,'%Y%m%d_%H%M').replace(tzinfo=dt.timezone.utc).timestamp()+150.
    lab=max((l for l,v in P1FILES.items() if v['created']<t), key=lambda l:P1FILES[l]['created'])
    d=[F[k+'_'+r]['dt_ms']-(L1 if r in V5_1 else 0.) for r in rx]
    o=[O[k+'_'+r+'_'+lab]['dt_ms'] for r in rx]
    rows.append(dict(campaign=1,pass_=k,receivers=rx,dt_ms=float(np.mean(d)),oem=P1FILES[lab]['name'],dt_oem_ms=float(np.mean(o)),
                     age_created_h=(t-P1FILES[lab]['created'])/3600.,age_start_h=(t-P1FILES[lab]['start'])/3600.,
                     night='p1_'+{'20261103_0040':'a','20261103_0220':'a','20261103_2355':'b','20261104_0130':'b','20261104_2305':'c'}[k]))
P={q['pass_']:q for q in Q['passes']}
for x in Q['prediction']:
    q=P[x['pass_']]
    rows.append(dict(campaign=2,pass_=x['pass_'],receivers=q['receivers'],dt_ms=x['dt_ms'],oem=x['oem'],dt_oem_ms=x['dt_oem_ms'],
                     age_created_h=x['age_h'],age_start_h=x['prop_age_h'],night='p2_%d'%q['night'],order=q['order'],
                     tle_post_ms=q['tle_post_ms']))
for r in rows: r['diff_ms']=r['dt_ms']-r['dt_oem_ms']
def st(v):
    v=np.array(v,float)
    if not len(v): return None
    return dict(n=int(len(v)),mean_ms=float(v.mean()),sd_ms=float(v.std(ddof=1)) if len(v)>1 else None,rms_ms=float(np.sqrt(np.mean(v**2))))
BAD='20261025T0255'
bins=[(0,24),(24,48),(48,1e9)]
def binned(key,excl=False):
    out={}
    for a,b in bins:
        for c in (1,2,'both'):
            v=[r['diff_ms'] for r in rows if a<=r[key]<b and (c=='both' or r['campaign']==c) and not (excl and r['oem']==BAD)]
            s=st(v)
            if s: out['%d-%s h|%s'%(a,'' if b>1e8 else b,c)]=s
    return out
fresh=[r for r in rows if r['age_start_h']<24]
X=np.array([r['dt_oem_ms'] for r in fresh]); Y=np.array([r['dt_ms'] for r in fresh]); sl,ic=np.polyfit(X,Y,1); res=Y-(sl*X+ic)
reg=dict(n=len(fresh),slope=float(sl),slope_se=float(np.sqrt(res@res/(len(X)-2)/np.sum((X-X.mean())**2))),intercept_ms=float(ic),r=float(np.corrcoef(X,Y)[0,1]))
# scatter within a night (around the night mean) and between nights, fresh predictions; in this data set a
# prediction file and a night coincide, so the two cannot be separated
g={}
for r in fresh: g.setdefault(r['night'],[]).append(r['diff_ms'])
within=[v-np.mean(vs) for vs in g.values() if len(vs)>1 for v in vs]; dof=sum(len(vs)-1 for vs in g.values() if len(vs)>1)
means=np.array([np.mean(vs) for vs in g.values()])
scatter=dict(n_nights=len(g),within_sd_ms=float(np.sqrt(np.sum(np.square(within))/dof)),within_dof=int(dof),
             between_sd_of_means_ms=float(means.std(ddof=1)),mean_of_night_means_ms=float(means.mean()),
             se_night_level_ms=float(means.std(ddof=1)/np.sqrt(len(means))),
             nights_one_file=all(len({r['oem'] for r in fresh if r['night']==n})==1 for n in g))
out=dict(rule='latest OEM file created before the pass; age from its START_TIME (propagation span) and from creation',
         latency_ms=dict(phase1_HU=L1,phase2_HUM=Q['latency_v5_HUM']['L_ms']),
         by_start_age=binned('age_start_h'),by_start_age_excl_20261025=binned('age_start_h',True),
         by_creation_age=binned('age_created_h'),fresh_regression=reg,fresh_scatter=scatter,rows=rows)
json.dump(out,open(R+'iss_oem_combined.json','w'),indent=1)
for k in ('by_start_age_excl_20261025','by_creation_age'):
    print(k); [print('  ',a,{kk:round(vv,1) if isinstance(vv,float) else vv for kk,vv in b.items()}) for a,b in out[k].items()]
print(reg); print(scatter)
for r in rows[:5]: print(r['pass_'],r['oem'],round(r['age_created_h'],1),round(r['age_start_h'],1),round(r['diff_ms'],1))
