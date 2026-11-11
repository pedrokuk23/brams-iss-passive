# Phase 2 search: the frozen statistic (fastsearch.py with burst blanking and limiter, PROTOCOL items 2-2c) and the
# sign-reversed null for every selected pass (phase2_plan.py, phase2_select.py) in every retrieved recording.
import os, sys, json, glob, time, numpy as np
shard,nsh=int(sys.argv[1]),int(sys.argv[2])
import bwav, fastsearch as F
from sgp4.api import Satrec
exec(open('synth.py').read().split("def search")[0])          # baseband()
from linefilt import maybe_clean
LIM={'ISS complex':1e9,'AST SpaceMobile':3000,'Starlink DTC':1000,'Starlink V2-mini':300,'Starlink V1.x':300}
WD=os.environ.get('P2DIR','../data/raw_wav_phase2/'); OUT=os.environ.get('OUTTAG','p2_shard')
path={os.path.basename(p)[11:31]:p for d in WD.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
if not path: sys.exit('phase2_search: no WAV files found in '+WD)
sel=json.load(open('../results/phase2_files.json'))
byfile={}
if os.path.exists('../results/phase2_byfile.json'):      # same selection, precomputed (phase2_plan_passes.json is 484 MB)
    byfile={k:v for k,v in json.load(open('../results/phase2_byfile.json')).items() if k in path}
else:
    for p in json.load(open('../results/phase2_plan_passes.json')):
        if p['file'] in path and p['sigma_min']<=min(LIM[p['cls']],p['sigma_max']): byfile.setdefault(p['file'],[]).append(p)
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '):
            nm=L[k][2:].strip()
            sets.setdefault(nm,[]).append((L[k+1],L[k+2]))
EP={}
def satrec(nm,jd0):
    if nm not in EP: EP[nm]=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets[nm]])
    return Satrec.twoline2rv(*sets[nm][int(np.argmin(abs(EP[nm]-jd0)))])
files=[f for f in sorted(byfile)][shard::nsh]
if os.environ.get('REVERSE'): files=files[::-1]
if os.environ.get('SUB'): a,b=map(int,os.environ['SUB'].split('/')); files=files[a::b]
of=f'../results/{OUT}{shard}.json'; out=json.load(open(of)) if os.path.exists(of) else []
done={r['file'] for g in glob.glob('../results/'+OUT+'*.json') for r in json.load(open(g)) if r['case']=='filedone'}; t00=time.time(); BUDGET=float(os.environ.get('BUDGET','1e9'))
for k in files:
    if k in done: continue
    if time.time()-t00>BUDGET: break
    md,x,_=bwav.read(path[k]); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    t0=md['start_us']/1e6; jd0=t0/86400+2440587.5
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    for p in byfile[k]:
        tg,fdg,rr=F.doppler(satrec(p['obj'],jd0),t0,TX,RX)
        r=F.search_fast(bb,ts,fs,fb,tg,fdg,full=True,blank=2.0,clip=20)
        if r: out.append(dict(case='pass',file=k,obj=p['obj'],cls=p['cls'],rrmin=p['rrmin'],n_valid_s=p['n_valid_s'],sigma_min=p['sigma_min'],z=r[0],tau=r[1],nu=r[2],t_pk=r[3],f_raw=r[4]))
        r=F.search_fast(bb,ts,fs,fb,tg,-fdg,full=True,blank=2.0,clip=20)
        if r: out.append(dict(case='null',file=k,obj=p['obj'],cls=p['cls'],z=r[0],tau=r[1],nu=r[2],t_pk=r[3],f_raw=r[4]))
    out.append(dict(case='filedone',file=k)); json.dump(out,open(of,'w')); print(k,len(out),flush=True)
