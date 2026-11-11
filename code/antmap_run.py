# Mapa de antena com aviões: SNR do eco bloco a bloco (1 s) ao longo de cada trajetória ADS-B,
# junto com a geometria (az/el vistos do receptor e do transmissor, distâncias, rumo do avião).
# Uso: python3 antmap_run.py <ESTAÇÃO> <ADSB.jsonl> <dir_wav> <saida.json>
import json,numpy as np,glob,sys,os,bwav
from scipy import signal
ST,ADSB,WD,OUT=sys.argv[1:5]
os.environ['ADSB']=ADSB
exec(open('air_overlay.py').read().split("files=")[0])
RMAX=80.0; TAUS=np.arange(-2,2.001,0.25); FMIN=3.0
def enu(X,ref,lat,lon):
    la,lo=np.radians(lat),np.radians(lon); d=X-ref
    E=-np.sin(lo)*d[...,0]+np.cos(lo)*d[...,1]
    N=-np.sin(la)*np.cos(lo)*d[...,0]-np.sin(la)*np.sin(lo)*d[...,1]+np.cos(la)*d[...,2]
    U=np.cos(la)*np.cos(lo)*d[...,0]+np.cos(la)*np.sin(lo)*d[...,1]+np.sin(la)*d[...,2]
    return np.degrees(np.arctan2(E,N))%360,np.degrees(np.arctan2(U,np.hypot(E,N)))
out=[]
for p in sorted(glob.glob(f'{WD}/RAD_BEDOUR_20261104_*_{ST}_SYS001.wav')):
    hh=os.path.basename(p).split('_')[3]   # file name only, so the repository path may contain underscores
    if hh[2:] not in ('00','30') or hh=='2130': continue
    md,x,pps=bwav.read(p); fs=md['samplerate']; t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>700)&(ff<1500); fb=ff[m][np.argmax(F[m])]
    ts=np.arange(len(x))/fs
    bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-signal.oaconvolve(bb,np.ones(w)/w,mode='same')
    bb=signal.decimate(signal.decimate(bb,4,ftype='fir'),4,ftype='fir'); fsd=fs/16; tsd=np.arange(len(bb))/fsd
    B=int(fsd); nb=len(bb)//B; tb=t0+(np.arange(nb)+0.5)
    nac=0
    for name,typ,tsa,la,lo,al in AC:
        tg=np.arange(-5,305,0.5)+t0
        sel=(tg>=tsa.min())&(tg<=tsa.max())
        if sel.sum()<200: continue
        tg=tg[sel]
        X=ecef(np.interp(tg,tsa,la),np.interp(tg,tsa,lo),np.interp(tg,tsa,al))
        R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1)
        if np.sqrt((R1*R2).min())>RMAX: continue
        fdg=-np.gradient(R1+R2,tg)/lam
        best=None
        for tau in TAUS:
            fd=np.interp(tsd+tau+t0,tg,fdg,left=np.nan,right=np.nan)
            ok=np.isfinite(fd); fd0=np.where(ok,fd,0)
            y=bb*np.exp(-2j*np.pi*np.cumsum(fd0)/fsd)
            Y=np.fft.fft(y[:nb*B].reshape(nb,B),axis=1)          # bins de 1 Hz
            Pw=np.abs(Y)**2
            sig=Pw[:,0]; far=np.r_[8:60,B-60:B-8]
            noise=np.median(Pw[:,far],axis=1)
            vb=ok[:nb*B].reshape(nb,B).all(1)&(np.abs(fd0)>FMIN)[:nb*B].reshape(nb,B).all(1)
            if vb.sum()<20: break
            snr=sig/noise
            score=np.mean(snr[vb])
            if best is None or score>best[0]: best=(score,tau,snr,vb,fd0[:nb*B].reshape(nb,B).mean(1))
        if best is None: continue
        score,tau,snr,vb,fdb=best
        k=np.where(vb)[0]; tt=tb[k]
        Xb=ecef(np.interp(tt,tsa,la),np.interp(tt,tsa,lo),np.interp(tt,tsa,al))
        r1=np.linalg.norm(Xb-TX,axis=1); r2=np.linalg.norm(Xb-RX,axis=1)
        azr,elr=enu(Xb,RX,md['ant_lat'],md['ant_long']); azt,elt=enu(Xb,TX,md['beacon_lat'],md['beacon_long'])
        hd=np.degrees(np.arctan2(np.gradient(np.interp(tt,tsa,lo))*np.cos(np.radians(np.interp(tt,tsa,la))),np.gradient(np.interp(tt,tsa,la))))%360
        u1=(Xb-TX)/r1[:,None]; u2=(Xb-RX)/r2[:,None]; beta=np.degrees(np.arccos(np.clip((u1*u2).sum(1),-1,1)))
        out.append(dict(st=ST,file=hh,ac=name,type=typ,tau=float(tau),score=float(score),
            t=np.round(tt-t0,1).tolist(),snr=np.round(snr[k],3).tolist(),fd=np.round(fdb[k],2).tolist(),
            r1=np.round(r1,2).tolist(),r2=np.round(r2,2).tolist(),az_rx=np.round(azr,1).tolist(),el_rx=np.round(elr,1).tolist(),
            az_tx=np.round(azt,1).tolist(),el_tx=np.round(elt,1).tolist(),hdg=np.round(hd,1).tolist(),beta=np.round(beta,1).tolist()))
        nac+=1
    print(ST,hh,nac,flush=True)
json.dump(dict(st=ST,ant=[md['ant_lat'],md['ant_long'],md['ant_alt']],tx=[md['beacon_lat'],md['beacon_long']],tracks=out),open(OUT,'w'))
print('FIM',ST,len(out))
