# Validation of search_fast against search() on known echoes and on null cases.
import os, json, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
import fastsearch as F
D=os.environ.get('WAVDIR','../data/raw_wav/')
TLE['ISS']=tuple(l.strip() for l in open('../data/tle/iss.txt').read().split('\n')[1:3])
cases=[]
for d in json.load(open('../results/iss_detection_table_v2.json')):
    cases.append((f"RAD_BEDOUR_{d['pass_']}_{d['st']}_SYS001.wav",'ISS',None))
for o in json.load(open('../results/close_scan.json')):
    if o['case']=='pass' and o['z']>=7: cases.append((f"RAD_BEDOUR_{o['file']}_SYS001.wav",o['obj'],None))
nl=[o for o in json.load(open('../results/close_scan.json')) if o['case']=='null'][:30]
for o in nl: cases.append((o['file'],o['obj'],f"RAD_BEDOUR_{o['geom']}_SYS001.wav"))
out=[]
for f,obj,geomf in cases:
    p=D+f; gp=D+(geomf or f)
    md,x,_=bwav.read(p); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    g=GeoSat(gp,obj); ex=search(bb,ts,fs,fb,g,F.TAUS,F.DFS)
    fa=F.search_fast(bb,ts,fs,fb,g.t,g.fd)
    out.append(dict(file=f,obj=obj,null=geomf is not None,exact=ex,fast=fa)); print(json.dumps(out[-1]),flush=True)
json.dump(out,open('../results/fast_validate.json','w'),indent=1)
