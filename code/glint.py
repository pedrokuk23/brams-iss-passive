# Search matched to short glints (tumbling or briefly specular objects). After dechirping along the element-set Doppler,
# the power of each 0.5-s block in each frequency cell is divided by the block's noise level (median over 14-40 Hz offsets),
# giving s = SNR_bin; window sums over W blocks give z_W = (sum_W s - W)/sqrt(W) under Gaussian noise; the statistic is the
# maximum over start block, frequency cell (|nu| <= 12 Hz), timing offset tau and W in {4, 8, 16, 32, 64} (2-32 s).
# The threshold comes from the same statistic computed with the same geometry on recordings without the object.
import os, json, glob, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
WS=[4,8,16,32,64]
def snr_map(bb,ts,fs,fb,g,tau,Tb=0.5,dfs=np.arange(-12,12.01,0.5)):
    fmax=0.9*min(fb,fs/2-fb); B=int(Tb*fs); nb=len(bb)//B
    fd=np.interp(ts+tau,g.t,g.fd); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    valid=(np.abs(fd)>30)&(np.abs(fd)<fmax); vb=valid[:nb*B].reshape(nb,B).all(1)
    Y=y[:nb*B].reshape(nb,B)
    E=np.exp(-2j*np.pi*np.outer(np.arange(B)/fs,dfs)).astype(np.complex64); P=np.abs(Y@E/B)**2
    F=np.abs(np.fft.fft(Y,axis=1)/B)**2; f=np.fft.fftfreq(B,1/fs); nm=(abs(f)>14)&(abs(f)<40)
    noise=np.median(F[:,nm],axis=1)/np.log(2)
    s=P/noise[:,None]; s[~vb]=np.nan
    return s, vb
def glint_stat(path,g,taus=np.arange(-1,1.01,0.5)):
    md,x,_=bwav.read(path); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    best=(-np.inf,None)
    for tau in taus:
        s,vb=snr_map(bb,ts,fs,fb,g,tau)
        if vb.sum()<8: continue
        for W in WS:
            if W>vb.sum(): continue
            c=np.nancumsum(np.nan_to_num(s-1),axis=0); n=np.cumsum(vb)
            Sw=c[W:]-c[:-W]; nw=n[W:]-n[:-W]; ok=nw==W
            z=np.where(ok[:,None],Sw/np.sqrt(W),-np.inf)
            i,j=np.unravel_index(np.argmax(z),z.shape)
            if z[i,j]>best[0]: best=(float(z[i,j]),dict(tau=float(tau),W=W,t0=float((i+1)*0.5),nu=float(np.arange(-12,12.01,0.5)[j])))
    return best
if __name__=='__main__':
    D=os.environ.get('WAVDIR','../data/raw_wav/'); EV='../data/raw_wav_envisat/'
    exec(open('envisat.py').read().split("rng=np.random")[0].split("import os, glob, json, numpy as np")[1])
    out=[]
    # sanity check on the known glint of SPACEMOBILE-008
    for st in ['BEOUDS','FRMAUB','BEOTTI']:
        p=D+f'RAD_BEDOUR_20261105_0310_{st}_SYS001.wav'; g=GeoSat(p,'SPACEMOBILE-008'); z,info=glint_stat(p,g,taus=np.array([-0.5]))
        out.append(dict(case='bluebird',file=os.path.basename(p),zg=z,**info)); print(json.dumps(out[-1]),flush=True)
    rng=np.random.default_rng(5)
    for pa,keys in PASSES.items():
        for st in STS:
            for k in keys:
                f=EV+f'RAD_BEDOUR_{k}_{st}_SYS001.wav'
                if not os.path.exists(f): continue
                g=geo(f); z,info=glint_stat(f,g)
                if info is None: continue
                out.append(dict(case='pass',pass_=pa,file=os.path.basename(f),zg=z,**info)); print(json.dumps(out[-1]),flush=True)
                for o in rng.choice(sorted(glob.glob(D+f'RAD_BEDOUR_*_{st}_SYS001.wav')),6,replace=False):
                    z,info=glint_stat(o,g)
                    if info is None: continue
                    out.append(dict(case='null',pass_=pa,file=os.path.basename(o),geom=os.path.basename(f),zg=z,**info))
        json.dump(out,open('../results/glint_envisat.json','w'),indent=1)
    nz=np.array([r['zg'] for r in out if r['case']=='null']); print('null n',len(nz),'pcts 50/90/99/max',np.percentile(nz,[50,90,99]).round(1),nz.max().round(1))
