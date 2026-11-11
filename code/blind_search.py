# Blind search over all candidate passes (PROTOCOL.md items 1-2) and the null searches (item 4), sharded by recording.
import os, sys, json, glob, numpy as np
shard,nsh=int(sys.argv[1]),int(sys.argv[2])          # read before the exec() calls below, which reset sys.argv
import bwav, fastsearch as F
from sgp4.api import Satrec
exec(open('synth.py').read().split("def search")[0])          # baseband()
from linefilt import maybe_clean
exec(open('blind_enumerate.py').read().split("names=[x[0]")[0]) # read_sets, catalog, special element sets
TLEc={nm:(l1,l2) for nm,l1,l2 in cat}
def satrec(obj,jd0):
    if obj in special:
        sets=special[obj]; ep=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets])
        return Satrec.twoline2rv(*sets[int(np.argmin(abs(ep-jd0)))])
    return Satrec.twoline2rv(*TLEc[obj])
NNULL=30
BLANK=float(os.environ.get('BLANK','2.0'))   # burst blanking (PROTOCOL item 2b)
CLIP=float(os.environ.get('CLIP','20'))      # limiter (PROTOCOL item 2c)
D=os.environ.get('WAVDIR','../data/raw_wav/'); EV=os.environ.get('EVDIR','../data/raw_wav_envisat/')
path={os.path.basename(p)[11:31]:p for d in (D+':'+EV).split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
P=[p for p in json.load(open('../results/blind_passes.json')) if p['n_valid_s']>=20 and p['rrmin_valid']<=600]
byfile={}; 
for p in P: byfile.setdefault(p['file'],[]).append(p)
ALL=sorted({p['file'] for p in json.load(open('../results/blind_passes.json'))}); files=[f for f in ALL[shard::nsh] if f in path]
if os.environ.get('REVERSE'): files=files[::-1]
if os.environ.get('ONLY'): files=[f for f in ALL if any(f.startswith(p) for p in os.environ['ONLY'].split(','))]
OUT=os.environ.get('OUTTAG','blind_shard')
rng=np.random.default_rng(1000+shard); out=[]
def tmin(k): return int(k[:8])*1440+int(k[9:11])*60+int(k[11:13])
for k in files:
    md,x,_=bwav.read(path[k]); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    t0=md['start_us']/1e6; jd0=t0/86400+2440587.5
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    for p in byfile.get(k,[]):
        tg,fdg,rr=F.doppler(satrec(p['obj'],jd0),t0,TX,RX)
        r=F.search_fast(bb,ts,fs,fb,tg,fdg,full=True,blank=BLANK,clip=CLIP)
        if r: out.append(dict(case='pass',file=k,obj=p['obj'],rrmin=p['rrmin_valid'],n_valid_s=p['n_valid_s'],z=r[0],tau=r[1],nu=r[2],t_pk=r[3],f_raw=r[4]))
        # null: the same Doppler history with its sign reversed, in the same recording; no object produces it
        r=F.search_fast(bb,ts,fs,fb,tg,-fdg,full=True,blank=BLANK,clip=CLIP)
        if r: out.append(dict(case='null',file=k,obj=p['obj'],rrmin=p['rrmin_valid'],n_valid_s=p['n_valid_s'],z=r[0],tau=r[1],nu=r[2],t_pk=r[3],f_raw=r[4]))
    print(k,len(out),flush=True)
    json.dump(out,open(f'../results/{OUT}{shard}.json','w'))
