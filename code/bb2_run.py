import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
exec(open('bb2.py').read())
names=[n for n in TLE if n.startswith('SPACEMOBILE') or n.startswith('BLUEWALKER')]
files=[WAVDIR_+f for f in json.load(open('bluebird_files.json'))['search']]   # fixed list of the BlueBird-pass recordings
out=[]
for f in files:
    md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs)
    for n in names:
        g=GeoSat(f,n); m=(g.t>=0)&(g.t<=300)
        if g.rr[m].min()>750: continue
        z,tau,df=search(bb,ts,fs,fb,g,TAUS,DFS)
        r=dict(file=f.split('/')[-1],sat=n,rrmin=float(g.rr[m].min()),tmin=float(g.t[m][np.argmin(g.rr[m])]),z=z,tau=tau,df=df); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/bb2_results.json','w'),indent=1)
