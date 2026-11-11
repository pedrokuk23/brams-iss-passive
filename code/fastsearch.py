# Fast form of the detection statistic of search() (synth.py), for searching many objects per recording.
# The recording is read, converted to baseband and cleaned once. For each object, it is dechirped along the predicted
# Doppler shifted by a few reference offsets tau0; the power of each 0.5-s block is taken from a 4x zero-padded FFT, whose
# 0.5-Hz bins coincide with the frequency grid of search(); a nearby offset tau is obtained by shifting each block's spectrum
# by fd(t_b+tau) - fd(t_b+tau0) at the block centre. Mean over valid blocks, median/MAD over the grid and maximum as in search().
import numpy as np
from sgp4.api import Satrec, jday
from sgp4.propagation import gstime
import datetime as dtm
LAM=299792458./49.97e6
TAUS=np.arange(-6,6.01,0.5); DFS=np.arange(-12,12.01,0.5); TAU0=np.array([-4.,0.,4.])
def ecef(lat,lon,h):
    a=6378.137;f=1/298.257223563;e2=f*(2-f);la,lo=np.radians(lat),np.radians(lon);N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
def doppler(sat,t0,TX,RX,tgrid=np.arange(-60,360,0.25)):
    jd=(t0+tgrid)/86400.0+2440587.5; jdi=np.floor(jd-0.5)+0.5; fr=jd-jdi
    e,r,v=sat.sgp4_array(jdi,fr); g=np.array([gstime(a+b) for a,b in zip(jdi,fr)])
    X=np.stack([r[:,0]*np.cos(g)+r[:,1]*np.sin(g),-r[:,0]*np.sin(g)+r[:,1]*np.cos(g),r[:,2]],1); R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1)
    return tgrid, -np.gradient(R1+R2,tgrid)*1e3/LAM, np.sqrt(R1*R2)
RING=np.r_[np.arange(-40,-12.9,0.5),np.arange(13,40.01,0.5)]     # residual frequencies just outside the search band
def search_fast(bb,ts,fs,fb,tg,fdg,taus=TAUS,dfs=DFS,Tb=0.5,full=False,blank=None,clip=None):
    fmax=0.9*min(fb,fs/2-fb); B=int(Tb*fs); nb=len(bb)//B; nfft=int(round(fs/0.5)); assert abs(nfft-4*B)<=4 or True
    df=fs/nfft; tb0=np.arange(nb)*B/fs; tb1=tb0+(B-1)/fs; tbc=tb0+0.5*Tb
    Z=np.full((len(taus),len(dfs)),np.nan); PB={}
    for t0 in TAU0:
        sel=np.abs(taus-t0)<=2.0+1e-9
        if not sel.any(): continue
        fd=np.interp(ts+t0,tg,fdg); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
        # blocks that can be valid for some tau in this group (others are never used)
        fA=np.interp(tb0+t0-2.5,tg,fdg); fB=np.interp(tb1+t0+2.5,tg,fdg); fC=np.interp(tbc+t0,tg,fdg)
        use=(np.maximum.reduce([abs(fA),abs(fB),abs(fC)])>30)&(np.minimum.reduce([abs(fA),abs(fB),abs(fC)])<fmax)
        if use.sum()<5: continue
        S=np.zeros((nb,nfft)); S[use]=np.abs(np.fft.fft(y[:nb*B].reshape(nb,B)[use],nfft,axis=1)/B)**2   # bins k*df, df = 0.5 Hz
        fref=np.interp(tbc+t0,tg,fdg)
        for a in np.flatnonzero(sel):
            tau=taus[a]
            fa=np.interp(tb0+tau,tg,fdg); fb_=np.interp(tb1+tau,tg,fdg)
            vb=(np.minimum(abs(fa),abs(fb_))>30)&(np.maximum(abs(fa),abs(fb_))<fmax)&(np.sign(fa)==np.sign(fb_))&use
            if vb.sum()<5: continue
            shift=np.interp(tbc+tau,tg,fdg)-fref                       # extra residual frequency at offset tau
            if blank or clip:
                kr=np.rint((RING[None,:]+shift[vb,None])/df).astype(int)%nfft
                b=np.median(np.take_along_axis(S[vb],kr,axis=1),axis=1)
            if blank:   # burst blanking: drop blocks whose power just outside the search band is high (broadband bursts)
                ok=b<=blank*np.median(b); idx=np.flatnonzero(vb); vb=vb.copy(); vb[idx[~ok]]=False; b=b[ok]
                if vb.sum()<5: continue
            k=np.rint((dfs[None,:]+shift[vb,None])/df).astype(int)%nfft
            Pb=np.take_along_axis(S[vb],k,axis=1)
            if clip:    # limiter: each block normalized by its own noise level and capped, so that no single block dominates
                Pb=np.minimum(Pb/(b[:,None]/np.log(2)),clip)
            Z[a]=Pb.mean(0); PB[a]=(np.flatnonzero(vb),Pb)
    v=Z[np.isfinite(Z)]
    if v.size==0: return None
    med=np.median(v); z=(Z-med)/(1.4826*np.median(np.abs(v-med)))
    i,j=np.unravel_index(np.nanargmax(z),z.shape)
    if not full: return float(z[i,j]),float(taus[i]),float(dfs[j])
    idx,Pb=PB[i]; c=np.convolve(Pb[:,j],np.ones(4),'same'); kb=idx[int(np.argmax(c))]   # 2-s window of maximum power
    tpk=tbc[kb]; fraw=float(np.interp(tpk+taus[i],tg,fdg)+dfs[j])
    return float(z[i,j]),float(taus[i]),float(dfs[j]),float(tpk),fraw
