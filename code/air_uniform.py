# Análise uniforme de aviões: mesmos parâmetros para todas as estações e dias
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import json,numpy as np,datetime as dt,bwav,glob,sys,os
from scipy import signal
ST,DAY=sys.argv[1],sys.argv[2]
os.environ['ADSB']={'20261103':'../data/adsb/adsb_be_nov03.jsonl','20261104':'../data/adsb/adsb_be_nov04.jsonl'}[DAY]
exec(open('air_overlay.py').read().split("files=")[0])
RMAX=70.0; TAUS=np.arange(-10,10.001,0.25); BLK=1.0; FMIN=3.0
out=[]
for p in sorted(glob.glob(f'{WAVDIR_}RAD_BEDOUR_{DAY}_10[0-2]?_{ST}_SYS001.wav')):
    md,x,pps=bwav.read(p); fs=md['samplerate']; t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
    n=len(x); ts=np.arange(n)/fs
    bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-signal.oaconvolve(bb,np.ones(w)/w,mode='same')
    bb=signal.decimate(signal.decimate(bb,4,ftype='fir'),4,ftype='fir'); fsd=fs/16; tsd=np.arange(len(bb))/fsd
    B=int(BLK*fsd); nb=len(bb)//B
    for name,typ,tsa,la,lo,al in AC:
        tg=np.arange(-20,320,0.5)+t0
        if tsa.min()>tg[0] or tsa.max()<tg[-1]: continue
        X=ecef(np.interp(tg,tsa,la),np.interp(tg,tsa,lo),np.interp(tg,tsa,al))
        R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1); rr=np.sqrt((R1*R2).min())
        if rr>RMAX: continue
        fdg=-np.gradient(R1+R2,tg)/lam
        sc=[]
        for tau in TAUS:
            fd=np.interp(tsd+tau+t0,tg,fdg); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fsd)
            P=np.abs(y[:nb*B].reshape(nb,B).mean(1))**2; vb=(np.abs(fd)>FMIN)[:nb*B].reshape(nb,B).all(1)
            sc.append(P[vb].mean() if vb.sum()>20 else np.nan)
        sc=np.array(sc)
        if np.isnan(sc).all(): continue
        z=(sc-np.nanmedian(sc))/(1.4826*np.nanmedian(abs(sc-np.nanmedian(sc)))+1e-30)
        i=int(np.nanargmax(z)); tr=TAUS[i]
        if 0<i<len(TAUS)-1 and np.all(np.isfinite(sc[i-1:i+2])):   # refinamento parabólico
            a,b,c=sc[i-1],sc[i],sc[i+1]; d=a-2*b+c
            if d<0: tr=TAUS[i]+0.25*0.5*(a-c)/d
        fdm=-np.interp(tsd+t0,tg,fdg); ym=bb*np.exp(-2j*np.pi*np.cumsum(fdm)/fsd)
        Pm=np.abs(ym[:nb*B].reshape(nb,B).mean(1))**2; vbm=(np.abs(fdm)>FMIN)[:nb*B].reshape(nb,B).all(1)
        g=10*np.log10(sc[i]/np.median(Pm[vbm])) if vbm.any() else np.nan
        out.append(dict(day=DAY,st=ST,file=os.path.basename(p).split('_')[3],ac=name,type=typ,tau=float(TAUS[i]),tau_ref=float(tr),z=float(z[i]),gain_db=float(g),rr_km=float(rr)))
json.dump(out,open(f'../results/uni_{ST}_{DAY}.json','w'))
print(ST,DAY,len(out),'associações')
