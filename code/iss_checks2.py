# Further checks of the ISS timing results:
# (1) delta t is a Doppler-equivalent time offset: compare, for the NASA ephemeris (OEM) against the element set, the
#     offset delta t_OEM fitted from the predicted Doppler with the purely geometric along-track difference of the two
#     orbits (oem_cmp_v2.json); the difference is what radial and cross-track differences add to delta t.
# (2) receiver offsets: delta t[pass, receiver] = a[pass] + b[group] (version-5 receivers; Ottignies), Oudsbergen as
#     reference, by least squares (rx_offsets.py); scatter about the model and reduced chi-square.
# (3) light time: the echo received at t left the target at t - R2/c, which delays the observed Doppler history by R2/c;
#     both predictions ignore it, so the measured delta t carries a bias of +R2/c that delta t_OEM (a difference of two
#     predictions) does not: delta t - delta t_OEM is biased by +R2/c (1.3-1.7 ms).
import json, numpy as np
c=299792.458
F=json.load(open('../results/iono_fit_v2.json')); O=json.load(open('../results/oem_cmp_v2.json'))
out={}
# (1)
rows=[]
for k,v in O.items():
    if not k.endswith('_1102') and not k.startswith('20261104_2305'): continue
    if k.startswith('20261104_2305') and not k.endswith('_1104'): continue
    rows.append(dict(key=k,dt_oem_ms=v['dt_ms'],along_geom_ms=v['along_ms'],diff_ms=v['dt_ms']-v['along_ms'],dpos_m=v['dpos_km']*1e3,
                     v_ground_kms=v['along_km']/v['along_ms']*1e3 if v['along_ms'] else None))
d=np.array([r['diff_ms'] for r in rows]); out['doppler_vs_geometric']=dict(rows=rows,rms_ms=float(np.sqrt(np.mean(d**2))),max_abs_ms=float(np.max(abs(d))),
    per_pass={p:float(np.mean([r['diff_ms'] for r in rows if r['key'].startswith(p)])) for p in sorted({r['key'][:13] for r in rows})},
    v_ground_kms=float(np.median([r['v_ground_kms'] for r in rows if r['v_ground_kms']])))
# (2) receiver offsets (rx_offsets.py): delta t[pass, receiver] = a[pass] + b[group], Oudsbergen as reference; Phase 1
#     has Humain and Uccle as version-5 receivers. 'latency_ms' is the offset of the version-5 receivers.
from rx_offsets import fit_offsets
V5=('BEHUMA','BEUCCE')
recs=[(k[:13],k[14:],v['dt_ms'],v['sdt_ms']) for k,v in F.items() if len(k)>13 and k[14:] in ('BEOTTI','BEOUDS','BEHUMA','BEUCCE')]
M=fit_offsets(recs,v5=V5)
P0=sorted({r[0] for r in recs}); pm={p:np.mean([r[2] for r in recs if r[0]==p]) for p in P0}
e0=np.array([r[2]-pm[r[0]] for r in recs])
out['latency_model']=dict(n=M['n'],n_pass=M['n_pass'],latency_ms=M['offsets_ms']['v5'],latency_se_ms=M['offsets_se_ms']['v5'],
    ottignies_offset_ms=M['offsets_ms']['ottignies'],ottignies_offset_se_ms=M['offsets_se_ms']['ottignies'],dof=M['dof'],
    resid_rms_N_ms=M['resid_rms_ms'],resid_sd_dof_ms=M['resid_sd_ms'],chi2_red=M['chi2_red'],formal_scale=float(np.sqrt(M['chi2_red'])),
    scatter_no_offsets_rms_N_ms=float(np.sqrt(np.mean(e0**2))),scatter_no_offsets_sd_dof_ms=float(np.sqrt(e0@e0/(len(e0)-len(P0)))),
    formal_err_ms=[float(min(r[3] for r in recs)),float(max(r[3] for r in recs))])
loo=[]
for p in P0:   # offsets fitted without a pass, scatter of that pass after correction
    Mo=fit_offsets([r for r in recs if r[0]!=p],v5=V5)
    cc=np.array([r[2]-(Mo['offsets_ms']['v5'] if r[1] in V5 else Mo['offsets_ms']['ottignies'] if r[1]=='BEOTTI' else 0) for r in recs if r[0]==p])
    loo+=list(cc-cc.mean())
out['latency_model']['loo_max_abs_ms']=float(np.max(np.abs(loo)))
# (3)
out['light_time']=dict(R2_km=[400,520],bias_ms=[400/c*1e3,520/c*1e3],note='in the measured delta t only; biases delta t - delta t_OEM by +R2/c')
json.dump(out,open('../results/iss_checks2.json','w'),indent=1)
print(json.dumps({k:(v if k!='doppler_vs_geometric' else {kk:vv for kk,vv in v.items() if kk!='rows'}) for k,v in out.items()},indent=1))
