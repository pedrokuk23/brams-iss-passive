# Glint-matched search (glint.py) applied to the close passes (< 400 km) at Oudsbergen, with the timing grid of the scan
# (tau in [-6, 6] s), for the DTC satellites (reference: echoes known) and for the other objects (older Starlink satellites
# and CubeSats, not detected by the window statistic). Null: each geometry applied to four other recordings of Oudsbergen.
import os, json, glob, numpy as np
exec(open('glint.py').read().split("if __name__")[0])
D=os.environ.get('WAVDIR','../data/raw_wav/'); TAUS=np.arange(-6,6.01,0.5)
cs=[r for r in json.load(open('../results/catalog_scan.json')) if not r['issc'] and r['rr']<400 and r['file'].endswith('BEOUDS')]
others=sorted(glob.glob(D+'RAD_BEDOUR_*_BEOUDS_SYS001.wav')); rng=np.random.default_rng(21); out=[]
for r in cs:
    p=D+f"RAD_BEDOUR_{r['file']}_SYS001.wav"; g=GeoSat(p,r['obj'])
    z,info=glint_stat(p,g,taus=TAUS)
    if info is None: continue
    o=dict(case='pass',obj=r['obj'],dtc='[DTC]' in r['obj'],file=r['file'],rrmin=float(g.rr.min()),zg=z,**info); out.append(o); print(json.dumps(o),flush=True)
    for f in rng.choice([f for f in others if r['file'][:13] not in f],4,replace=False):
        z,info=glint_stat(f,g,taus=TAUS)
        if info: out.append(dict(case='null',obj=r['obj'],file=os.path.basename(f),geom=r['file'],zg=z,**info))
    json.dump(out,open('../results/glint_close.json','w'),indent=1)
N=np.array([o['zg'] for o in out if o['case']=='null'])
print('null n',len(N),'p50/p90/p99/max',np.percentile(N,[50,90,99]).round(1),N.max().round(1))
for dtc in (True,False):
    P=[o for o in out if o['case']=='pass' and o['dtc']==dtc]
    print('DTC' if dtc else 'other',len(P),'above null max:',[(o['obj'],round(o['zg'],1),o['W'],o['tau']) for o in P if o['zg']>N.max()])
