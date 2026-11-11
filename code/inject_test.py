# End-to-end injection test of the blind-search statistic. Synthetic echoes that follow the
# predicted bistatic Doppler of a real Phase 2 pass are added to the baseband of real Phase 2 recordings (Oudsbergen
# and Humain, the receivers of most Phase 2 passes) before line removal, and the frozen search (fastsearch.py with
# burst blanking and limiter, as in phase2_search.py) is run. Two envelopes: a steady echo over the usable arc, and a
# glint with a sinc^2 lobe of full width 20 s (first null to first null) centred in the arc. The echo energy E/N0 is set
# relative to the threshold of Eq. (zE), z_th sqrt(n_b), with the receiver's Phase 2 threshold and n_b the number of
# usable 0.5-s blocks. Output: detection probability (z >= z_th at |tau| <= 0.5 s and |nu - 0| <= 1 Hz) and z / z_pred.
import sys
import os, json, glob, numpy as np
import bwav, fastsearch as F
from sgp4.api import Satrec
exec(open('synth.py').read().split("def search")[0])          # baseband()
from linefilt import maybe_clean
WD=os.environ.get('P2DIR','../data/raw_wav_phase2/')
path={os.path.basename(p)[11:31]:p for d in WD.split(':') if d for p in glob.glob(os.path.join(d,'*.wav'))}
if not path: sys.exit(f'{__file__}: no WAV files found in '+WD)
BY=json.load(open('../results/phase2_byfile.json')); S=json.load(open('../results/phase2_summary.json'))
ZTH={k:v['z_th'] for k,v in S['thresholds'].items()}; DET={(d['file'],d['obj']) for d in S['detections']}
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '): sets.setdefault(L[k][2:].strip(),[]).append((L[k+1],L[k+2]))
def satrec(nm,jd0):
    E=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets[nm]])
    return Satrec.twoline2rv(*sets[nm][int(np.argmin(abs(E-jd0)))])
rng=np.random.default_rng(11)
files=sorted(k for k in path if k[-6:] in ('BEOUDS','BEHUMA') and k in BY)
files=[files[i] for i in sorted(rng.choice(len(files),int(os.environ.get('NFILE','24')),replace=False))]
LEVELS=[0.5,0.75,1.0,1.25,1.5,2.0,3.0]; TG=20.0; out=[]
for k in files:
    md,x,_=bwav.read(path[k]); fs=md['samplerate']; fb,bb0,ts=baseband(x,fs)
    t0=md['start_us']/1e6; jd0=t0/86400+2440587.5
    TX=F.ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=F.ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    cand=[p for p in BY[k] if (k,p['obj']) not in DET and p['obj'] in sets]
    if not cand: continue
    p=cand[int(rng.integers(len(cand)))]; tg,fdg,rr=F.doppler(satrec(p['obj'],jd0),t0,TX,RX)
    fmax=0.9*min(fb,fs/2-fb); fd=np.interp(ts,tg,fdg); band=(abs(fd)>30)&(abs(fd)<fmax)
    if band.sum()<20*fs: continue
    nb=int(band.sum()/fs/0.5); tb=ts[band]; tc=0.5*(tb[0]+tb[-1])
    # noise density along the echo track, block by block: the recording is dechirped along the prediction and the noise
    # per 0.5-s block is the median power 13-40 Hz from the track (as the limiter of the search measures it); the
    # receiver noise is colored (it falls by 10 dB toward the edge of the audio band), so E/N0 is the sum of the
    # per-block signal-to-noise ratios with the local noise, as in Eq. (zE)
    cn=maybe_clean(bb0.copy(),ts,fs); B=int(0.5*fs); nbk=len(cn)//B; ph=2*np.pi*np.cumsum(fd)/fs
    Yb=np.abs(np.fft.fft((cn*np.exp(-1j*ph))[:nbk*B].reshape(nbk,B),4*B,axis=1)/B)**2; ffb=np.fft.fftfreq(4*B,1/fs)
    ring=(abs(ffb)>=13)&(abs(ffb)<=40); sig2b=B*np.median(Yb[:,ring],axis=1)/np.log(2)     # per-sample noise power, per block
    sig2=np.repeat(sig2b,B); sig2=np.r_[sig2,np.full(len(ts)-len(sig2),sig2b[-1])]
    for env in ('steady','glint'):
        w=band.astype(float) if env=='steady' else band*np.sinc((ts-tc)/(TG/2))**2
        for lev in LEVELS:
            EN0=lev*ZTH[k[-6:]]*np.sqrt(nb)               # target E/N0
            A2=EN0/np.sum(w**2/sig2)                       # E/N0 = A^2 sum(w^2 / sigma2_local) (per-sample SNR summed)
            y=bb0+np.sqrt(A2)*w*np.exp(1j*(ph+rng.uniform(0,2*np.pi)))
            y=maybe_clean(y,ts,fs)
            r=F.search_fast(y,ts,fs,fb,tg,fdg,blank=2.0,clip=20)
            if not r: continue
            z,tau,nu=r; zp=EN0/np.sqrt(nb)
            out.append(dict(file=k,obj=p['obj'],env=env,level=lev,EN0=float(EN0),nb=nb,z=z,tau=tau,nu=nu,z_pred=float(zp),
                            det=bool(z>=ZTH[k[-6:]] and abs(tau)<=0.5 and abs(nu)<=1.0)))
    print(k,p['obj'],nb,flush=True)
summ={}
for env in ('steady','glint'):
    summ[env]={}
    for lev in LEVELS:
        q=[o for o in out if o['env']==env and o['level']==lev]
        if q: summ[env][str(lev)]=dict(n=len(q),pd=float(np.mean([o['det'] for o in q])),median_z_over_pred=float(np.median([o['z']/o['z_pred'] for o in q])))
json.dump(dict(summary=summ,glint_full_width_s=TG,rows=out),open('../results/inject_test.json','w'),indent=1)
print(json.dumps(summ,indent=1))
