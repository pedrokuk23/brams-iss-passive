# Envisat (NORAD 27386, inactive, tumbling): search of the five best passes of 3-5 Nov 2026 with the ISS statistic and grid,
# with line removal at the SDR receivers. Element set: the one with epoch closest to each pass (data/tle/envisat_gp_history.txt).
# Null: the same pass geometry applied to recordings of the same receiver at other times (no Envisat), 8 per receiver and pass.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os, glob, json, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
EV=os.environ.get('EVDIR','../data/raw_wav_envisat/'); D=os.environ.get('WAVDIR','../data/raw_wav/')
L=[l.strip() for l in open('../data/tle/envisat_gp_history.txt') if l.strip()]
SETS=[(L[k+1],L[k+2]) for k in range(0,len(L),3)]
EP=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in SETS])
PASSES={'20261103_0241':['20261103_0235','20261103_0240'],'20261104_0203':['20261104_0200'],'20261104_1322':['20261104_1320'],
        '20261105_0304':['20261105_0300','20261105_0305'],'20261105_1244':['20261105_1240','20261105_1245']}
STS=['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']
TA=np.arange(-6,6.01,0.5); DF=np.arange(-12,12.01,0.5)
def geo(path):
    md,_,_=bwav.read(path); jd=md['start_us']/1e6/86400+2440587.5
    TLE['ENVISAT']=SETS[int(np.argmin(abs(EP-jd)))]; return GeoSat(path,'ENVISAT')
def run(path,g,st):
    md,x,_=bwav.read(path); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    try: return search(bb,ts,fs,fb,g,TA,DF)
    except ValueError: return None      # Doppler outside the usable band during the whole file
rng=np.random.default_rng(3); out=[]
for p,keys in PASSES.items():
    for st in STS:
        for k in keys:
            f=EV+f'RAD_BEDOUR_{k}_{st}_SYS001.wav'
            if not os.path.exists(f): continue
            g=geo(f); res=run(f,g,st)
            if res is None: print('skip',os.path.basename(f)); continue
            z,tau,df=res
            r=dict(case='pass',pass_=p,file=os.path.basename(f),st=st,rrmin=float(g.rr.min()),z=z,tau=tau,df=df); out.append(r); print(json.dumps(r),flush=True)
            others=sorted(glob.glob(D+f'RAD_BEDOUR_*_{st}_SYS001.wav'))
            for o in rng.choice(others,4,replace=False):
                res=run(o,g,st)
                if res is None: continue
                z,tau,df=res
                r=dict(case='null',pass_=p,file=os.path.basename(o),geom=os.path.basename(f),st=st,z=z,tau=tau,df=df); out.append(r)
json.dump(out,open('../results/envisat.json','w'),indent=1)
nul=[r for r in out if r['case']=='null']; print('null n',len(nul),'max z |tau|<=0.5:',max([r['z'] for r in nul if abs(r['tau'])<=0.5]+[0]),'n z>=7 & |tau|<=.5:',sum(r['z']>=7 and abs(r['tau'])<=0.5 for r in nul))
