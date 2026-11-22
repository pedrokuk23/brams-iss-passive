# Unconditional test for unintended emission. For every Phase 2 pass of a Starlink V2-mini or
# DTC satellite in the recordings at hand, the frozen statistic (fastsearch.py, blanking and limiter as in the blind
# search) is computed along the ONE-WAY Doppler history of the satellite-receiver leg, -(1/lambda) dR2/dt, as an
# emission near the beacon frequency would follow, and along its sign-reversed version as a null. Thresholds are set
# on these one-way nulls per receiver with the rule of the blind search (quantile 1e-3). Unlike uemr_test.py, this
# search does not start from the bistatic detections, so it can find emission lines that the bistatic search misses.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os, sys, json, glob, numpy as np
ARGS=sys.argv[1:]
import bwav, fastsearch as F
from sgp4.api import Satrec
from sgp4.propagation import gstime
exec(open('synth.py').read().split("def search")[0])          # baseband()
from linefilt import maybe_clean
shard,nsh=(int(ARGS[0]),int(ARGS[1])) if len(ARGS)>1 else (0,1)
WD=os.environ.get('P2DIR','../data/raw_wav_phase2/')
path={os.path.basename(p)[11:31]:p for d in WD.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
if not path: sys.exit(f'{__file__}: no WAV files found in '+WD)

# The published run used the Phase 2 recordings of 19-26 October that were at hand (every third one; frozen list below); SUBSET=none uses
# every file found.
_SUB=os.environ.get('SUBSET','../data/raw_wav_phase2/SUBSET_UEMR_SEARCH.txt')
if _SUB!='none':
    _keep={l.strip()[11:31] for l in open(_SUB) if l.strip()}
    path={k:v for k,v in path.items() if k in _keep}
    if len(path)<len(_keep): print(f'warning: {len(_keep)-len(path)} files of {_SUB} are missing; results will differ')
BY=json.load(open('../data/planning/phase2_byfile.json'))
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '): sets.setdefault(L[k][2:].strip(),[]).append((L[k+1],L[k+2]))
EP={}
def satrec(nm,jd0):
    if nm not in EP: EP[nm]=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets[nm]])
    return Satrec.twoline2rv(*sets[nm][int(np.argmin(abs(EP[nm]-jd0)))])
def oneway(sat,t0,RX,tg):
    jd=(t0+tg)/86400.0+2440587.5; ji=np.floor(jd-0.5)+0.5; fr=jd-ji; e,r,v=sat.sgp4_array(ji,fr); g=np.array([gstime(a+b) for a,b in zip(ji,fr)])
    X=np.stack([r[:,0]*np.cos(g)+r[:,1]*np.sin(g),-r[:,0]*np.sin(g)+r[:,1]*np.cos(g),r[:,2]],1)
    return -np.gradient(np.linalg.norm(X-RX,axis=1),tg)*1e3/F.LAM
files=sorted(k for k in path if k in BY)[::int(os.environ.get('EVERY','1'))][shard::nsh]   # with SUBSET=none, EVERY=3 takes every third recording
of=f"../results/uemr_search{os.environ.get('OUTTAG','')}_{shard}.json"; out=json.load(open(of)) if os.path.exists(of) else []
done={r['file'] for r in out if r.get('case')=='filedone'}
for k in files:
    if k in done: continue
    md,x,_=bwav.read(path[k]); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    t0=md['start_us']/1e6; jd0=t0/86400+2440587.5
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    for p in BY[k]:
        if p['cls'] not in ('Starlink V2-mini','Starlink DTC') or p['obj'] not in sets: continue
        sat=satrec(p['obj'],jd0); tg,fdg,rr=F.doppler(sat,t0,TX,RX); f1=oneway(sat,t0,RX,tg)
        for case,h in (('pass',f1),('null',-f1)):
            r=F.search_fast(bb,ts,fs,fb,tg,h,full=True,blank=2.0,clip=20)
            if r: out.append(dict(case=case,file=k,obj=p['obj'],cls=p['cls'],z=r[0],tau=r[1],nu=r[2],t_pk=r[3],f_raw=r[4]))
    out.append(dict(case='filedone',file=k)); json.dump(out,open(of,'w')); print(k,len(out),flush=True)
