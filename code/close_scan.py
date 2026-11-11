# Repeats the catalog scan for the closest passes (geometric mean range below 400 km, mostly Starlink) at the two SDR
# receivers of the scan, now with line removal and the full ISS search grid; null: same geometry on two other recordings.
import os, json, glob, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
D=os.environ.get('WAVDIR','../data/raw_wav/'); TA=np.arange(-6,6.01,0.5); DF=np.arange(-12,12.01,0.5)
cs=[r for r in json.load(open('../results/catalog_scan.json')) if not r['issc'] and r['rr']<400 and r['file'][-6:] in ('BEOUDS','FRMAUB')]
def run(path,g):
    md,x,_=bwav.read(path); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    try: return search(bb,ts,fs,fb,g,TA,DF)
    except ValueError: return None
rng=np.random.default_rng(11); out=[]
for r in cs:
    p=D+f"RAD_BEDOUR_{r['file']}_SYS001.wav"; g=GeoSat(p,r['obj']); res=run(p,g)
    if res is None: continue
    z,tau,df=res; o=dict(case='pass',file=r['file'],obj=r['obj'],rrmin=float(g.rr.min()),z=z,tau=tau,df=df); out.append(o); print(json.dumps(o),flush=True)
    st=r['file'][-6:]
    for f in rng.choice([f for f in sorted(glob.glob(D+f'RAD_BEDOUR_*_{st}_SYS001.wav')) if r['file'][:13] not in f],2,replace=False):
        res=run(f,g)
        if res: out.append(dict(case='null',file=os.path.basename(f),geom=r['file'],obj=r['obj'],z=res[0],tau=res[1],df=res[2]))
json.dump(out,open('../results/close_scan.json','w'),indent=1)
P=[o for o in out if o['case']=='pass']; N=[o for o in out if o['case']=='null']
crit=lambda o:o['z']>=7 and abs(o['tau'])<=0.5 and abs(o['df'])<12
print('pass',len(P),'detections',[(o['obj'],o['file'],round(o['z'],1)) for o in P if crit(o)],'| null',len(N),'false',sum(crit(o) for o in N))
