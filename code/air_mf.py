# Matched filter of every ADS-B aircraft near the link in the Ottignies files of 3 Nov 2026, 10:00-10:25 UTC (Fig. 2):
# best delay tau and gain over the mirrored-Doppler control; writes ../results/air_mf_rows.npy, read by paper_figs3.py.
import json,numpy as np,datetime as dt,bwav,glob,os
os.environ.setdefault('ADSB','../data/adsb/adsb_be_nov03.jsonl')
from scipy import signal
exec(open('air_overlay.py').read().split("files=")[0])
rows=[]
for p in sorted(glob.glob(os.environ.get('WAVDIR','../data/raw_wav/')+'RAD_BEDOUR_20261103_10[0-2]?_BEOTTI_SYS001.wav')):
    md,x,pps=bwav.read(p); fs=md['samplerate']; t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    F=np.abs(np.fft.rfft(x))**2; ff=np.fft.rfftfreq(len(x),1/fs); m=(ff>900)&(ff<1400); fb=ff[m][np.argmax(F[m])]
    n=len(x); ts=np.arange(n)/fs
    bb=signal.hilbert(x)*np.exp(-2j*np.pi*fb*ts); w=int(fs); bb=bb-signal.oaconvolve(bb,np.ones(w)/w,mode='same')
    bb=signal.decimate(signal.decimate(bb,4,ftype='fir'),4,ftype='fir'); fsd=fs/16; tsd=np.arange(len(bb))/fsd
    for name,typ,tsa,la,lo,al in AC:
        tg=np.arange(-30,330,0.5)+t0
        if tsa.min()>tg[0] or tsa.max()<tg[-1]: continue
        X=ecef(np.interp(tg,tsa,la),np.interp(tg,tsa,lo),np.interp(tg,tsa,al))
        R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1)
        if (R1*R2).min()>150**2: continue
        fdg=-np.gradient(R1+R2,tg)/lam
        taus=np.arange(-15,15.01,0.25); sc=[]
        for tau in taus:
            fd=np.interp(tsd+tau+t0,tg,fdg)
            y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fsd); valid=np.abs(fd)>3
            B=int(1.0*fsd); nb=len(y)//B; P=np.abs(y[:nb*B].reshape(nb,B).mean(1))**2; vb=valid[:nb*B].reshape(nb,B).all(1)
            sc.append(P[vb].mean() if vb.sum()>20 else np.nan)
        sc=np.array(sc)
        if np.isnan(sc).all(): continue
        # ruído: mesma métrica com Doppler espelhado (-fd), não físico
        fdm=-np.interp(tsd+t0,tg,fdg); ym=bb*np.exp(-2j*np.pi*np.cumsum(fdm)/fsd); vm=np.abs(fdm)>3
        Pm=np.abs(ym[:nb*B].reshape(nb,B).mean(1))**2; vbm=vm[:nb*B].reshape(nb,B).all(1)
        i=np.nanargmax(sc); snr=10*np.log10(sc[i]/np.median(Pm[vbm])) if vbm.any() else np.nan
        rows.append((os.path.basename(p).split('_')[3],name,typ,taus[i],snr,np.sqrt((R1*R2).min()),sc))
rows.sort(key=lambda r:-r[4])
print("arquivo avião tipo  tau_melhor(s)  ganho_vs_espelho(dB)  sqrt(R1R2)min km")
for r in rows[:25]: print(f"{r[0]} {str(r[1]):8} {str(r[2]):5} {r[3]:+6.2f}  {r[4]:6.1f}  {r[5]:6.1f}")
np.save('../results/air_mf_rows.npy',np.array([(r[0],r[1],r[2],r[3],r[4],r[5]) for r in rows],dtype=object),allow_pickle=True)
