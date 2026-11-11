# Are the Starlink detections reflections of the beacon or unintended emission of the satellites (UEMR, reported for
# V2-mini at 40-70 MHz by Bassa et al. 2024)? A reflection of the beacon follows the bistatic Doppler of both legs,
# f_D = -(1/lambda) d(R1 + R2)/dt, relative to the direct beacon; an emission from the satellite follows the one-way
# Doppler of the satellite-receiver leg, -(1/lambda) dR2/dt, about half as steep in this near-backscatter geometry, at a
# frequency of its own. For every Phase 2 detection of a Starlink satellite whose recording is at hand, the frozen
# statistic (fastsearch.py, blanking and limiter as in the search, same tau and nu grid) is computed along the bistatic
# history and along the one-way history shifted in frequency so that the two coincide at the detected peak (time t_pk,
# offset tau): if the line were an emission, the one-way history would follow it over the arc and give the higher z.
import sys
import os, json, glob, numpy as np
import bwav, fastsearch as F
from sgp4.api import Satrec
exec(open('synth.py').read().split("def search")[0])          # baseband()
from linefilt import maybe_clean
WD=os.environ.get('P2DIR','../data/raw_wav_phase2/')
path={os.path.basename(p)[11:31]:p for d in WD.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
if not path: sys.exit(f'{__file__}: no WAV files found in '+WD)

# The published run used the Phase 2 recordings of 19-26 October that were at hand (frozen list below); SUBSET=none uses
# every file found.
_SUB=os.environ.get('SUBSET','../data/raw_wav_phase2/SUBSET_CHECKS.txt')
if _SUB!='none':
    _keep={l.strip()[11:31] for l in open(_SUB) if l.strip()}
    path={k:v for k,v in path.items() if k in _keep}
    if len(path)<len(_keep): print(f'warning: {len(_keep)-len(path)} files of {_SUB} are missing; results will differ')
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '): sets.setdefault(L[k][2:].strip(),[]).append((L[k+1],L[k+2]))
EP={}
def satrec(nm,jd0):
    if nm not in EP: EP[nm]=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets[nm]])
    return Satrec.twoline2rv(*sets[nm][int(np.argmin(abs(EP[nm]-jd0)))])
from sgp4.propagation import gstime
def oneway(sat,t0,RX,tg):
    jd=(t0+tg)/86400.0+2440587.5; ji=np.floor(jd-0.5)+0.5; fr=jd-ji; e,r,v=sat.sgp4_array(ji,fr); g=np.array([gstime(a+b) for a,b in zip(ji,fr)])
    X=np.stack([r[:,0]*np.cos(g)+r[:,1]*np.sin(g),-r[:,0]*np.sin(g)+r[:,1]*np.cos(g),r[:,2]],1)
    return -np.gradient(np.linalg.norm(X-RX,axis=1),tg)*1e3/F.LAM
S=json.load(open('../results/phase2_summary.json'))
D=[d for d in S['detections'] if d['class'] in ('Starlink V2-mini','Starlink DTC') and d['file'] in path]
rng=np.random.default_rng(3); D=[D[i] for i in sorted(rng.choice(len(D),min(len(D),int(os.environ.get('NMAX','200'))),replace=False))]
out=[]; cache={}
for d in D:
    k=d['file']
    if k not in cache:
        md,x,_=bwav.read(path[k]); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
        cache={k:(md,fs,fb,bb,ts)}
    md,fs,fb,bb,ts=cache[k]; t0=md['start_us']/1e6; jd0=t0/86400+2440587.5
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    sat=satrec(d['obj'],jd0); tg,fdg,rr=F.doppler(sat,t0,TX,RX); f1=oneway(sat,t0,RX,tg)
    tp=d['t_pk']+d['tau']; f1s=f1+(np.interp(tp,tg,fdg)-np.interp(tp,tg,f1))
    zb=F.search_fast(bb,ts,fs,fb,tg,fdg,blank=2.0,clip=20); z1=F.search_fast(bb,ts,fs,fb,tg,f1s,blank=2.0,clip=20)
    if zb and z1: out.append(dict(file=k,obj=d['obj'],cls=d['class'],z_search=d['z'],nu_search=d['nu'],z_bistatic=zb[0],z_oneway=z1[0],nu_oneway=z1[2]))
zb=np.array([o['z_bistatic'] for o in out]); z1=np.array([o['z_oneway'] for o in out])
summ=dict(n=len(out),frac_bistatic_higher=float(np.mean(zb>z1)),median_z_bistatic=float(np.median(zb)),median_z_oneway=float(np.median(z1)),
          median_ratio=float(np.median(z1/zb)),frac_oneway_ge_z_bistatic=float(np.mean(z1>=zb)),
          nu_search_abs_median_Hz=float(np.median([abs(o['nu_search']) for o in out])))
json.dump(dict(summary=summ,rows=out),open('../results/uemr_test.json','w'),indent=1)
print(json.dumps(summ,indent=1))
