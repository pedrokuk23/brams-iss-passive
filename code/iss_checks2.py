# Further checks of the ISS timing results:
# (1) delta t is a Doppler-equivalent time offset: compare, for the NASA ephemeris (OEM) against the element set, the
#     offset delta t_OEM fitted from the predicted Doppler with the purely geometric along-track difference of the two
#     orbits (oem_cmp_v2.json); the difference is what radial and cross-track differences add to delta t.
# (2) receiver latency: linear model delta t[pass, receiver] = a[pass] + b * (version-5 receiver) + e, fitted by least
#     squares; standard error of b and scatter about the model with the degrees of freedom accounted for; formal errors
#     rescaled by sqrt(chi2_red).
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
# (2)
V5={'BEHUMA','BEUCCE','FRMAUB'}
recs=[(k[:13],k[14:],v['dt_ms'],v['sdt_ms']) for k,v in F.items() if len(k)>13 and k[14:] in ('BEOTTI','BEOUDS','BEHUMA','BEUCCE')]
passes=sorted({r[0] for r in recs}); X=np.zeros((len(recs),len(passes)+1)); y=np.array([r[2] for r in recs]); s=np.array([r[3] for r in recs])
for i,r in enumerate(recs): X[i,passes.index(r[0])]=1; X[i,-1]=1 if r[1] in V5 else 0
beta,res,rk,_=np.linalg.lstsq(X,y,rcond=None); e=y-X@beta; dof=len(y)-X.shape[1]; s2=e@e/dof
cov=s2*np.linalg.inv(X.T@X); chi2r=float(np.sum((e/s)**2)/dof)
# scatter about the pass mean without the latency term (what the text calls 4.3 ms), with N and with N - p
X0=X[:,:-1]; b0=np.linalg.lstsq(X0,y,rcond=None)[0]; e0=y-X0@b0
out['latency_model']=dict(n=len(y),n_pass=len(passes),latency_ms=float(beta[-1]),latency_se_ms=float(np.sqrt(cov[-1,-1])),dof=dof,
    resid_rms_N_ms=float(np.sqrt(np.mean(e**2))),resid_sd_dof_ms=float(np.sqrt(s2)),chi2_red=chi2r,formal_scale=float(np.sqrt(chi2r)),
    scatter_no_latency_rms_N_ms=float(np.sqrt(np.mean(e0**2))),scatter_no_latency_sd_dof_ms=float(np.sqrt(e0@e0/(len(y)-X0.shape[1]))),
    formal_err_ms=[float(s.min()),float(s.max())])
# (3)
out['light_time']=dict(R2_km=[400,520],bias_ms=[400/c*1e3,520/c*1e3],note='in the measured delta t only; biases delta t - delta t_OEM by +R2/c')
json.dump(out,open('../results/iss_checks2.json','w'),indent=1)
print(json.dumps({k:(v if k!='doppler_vs_geometric' else {kk:vv for kk,vv in v.items() if kk!='rows'}) for k,v in out.items()},indent=1))
