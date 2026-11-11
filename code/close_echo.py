# Echo properties of the close-pass candidates: duration of the echo along the best (tau, nu), measured C/N0 and E/N0,
# and the effective cross section from the radar equation with the receiver gain calibrated on SPACEMOBILE-008.
# An echo is accepted when it lasts at least 5 s along the prediction (a meteor burst lasts ~1 s and spans all frequencies).
import os, json, numpy as np
exec(open('detectability.py').read().split("out=[]")[0])
T=json.load(open('../results/detect_theory.json'))
def cn0(sigma,R,Gst):   # same radar equation as detect_theory.py (kept local to avoid name clashes)
    return T['Pt_W']*Gst*T['lam_m']**2*sigma/((4*np.pi)**3*R**4*1.380649e-23*T['Tsky_K'])
G=10**(T['Geff_dB']['BEOUDS']/10)
cs=json.load(open('../results/close_scan.json')); out=[]
for o in [o for o in cs if o['case']=='pass' and o['z']>=7]:
    st=o['file'][-6:]; p=D+f"RAD_BEDOUR_{o['file']}_SYS001.wav"
    r=blocks(p,o['obj'],o['tau'],o['df'])
    if r['nb']<8: print('no valid blocks',o['obj'],o['file']); continue
    s=summarize(r)
    Gst=G if st=='BEOUDS' else 10**(T['Geff_dB'][st]/10)
    sig=10**(s['cn0_peak_dBHz']/10)/cn0(1.,s['rr_at_peak']*1e3,Gst)
    e=dict(obj=o['obj'],file=o['file'],z=o['z'],tau=o['tau'],nu=o['df'],dtc='[DTC]' in o['obj'],duration_s=s['duration_s'],
           cn0_dBHz=s['cn0_peak_dBHz'],E_N0_dB=s['E_N0_dB'],rr_km=s['rr_at_peak'],sigma_m2=float(sig),echo=s['duration_s']>=5)
    out.append(e); print(json.dumps({k:(round(v,1) if isinstance(v,float) else v) for k,v in e.items()}),flush=True)
P=[o for o in cs if o['case']=='pass']; N=[o for o in cs if o['case']=='null']
summ=dict(n_pass=len(P),n_null=len(N),null_max_z_any_tau=float(max(o['z'] for o in N)),null_max_z_tau05=float(max(o['z'] for o in N if abs(o['tau'])<=0.5)),
  by_type={})
for st in ['BEOUDS','FRMAUB']:
    for dtc in [True,False]:
        pp=[o for o in P if o['file'].endswith(st) and ('[DTC]' in o['obj'])==dtc]
        det=[e for e in out if e['file'].endswith(st) and e['dtc']==dtc and e['echo']]
        summ['by_type'][f"{st}_{'DTC' if dtc else 'other'}"]=dict(n=len(pp),detected=len(det),rr=[float(min(o['rrmin'] for o in pp)),float(max(o['rrmin'] for o in pp))] if pp else None)
json.dump(dict(echoes=out,summary=summ),open('../results/close_echo.json','w'),indent=1)
print(json.dumps(summ,indent=1))
