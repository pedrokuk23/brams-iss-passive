# Measured carrier-to-noise density of the echoes along the best (tau, nu) of the search, and the relation between the
# detection statistic z and the echo energy. For a 0.5-s block and the |Y E/B|^2 normalization of search(), the power in the
# echo bin relative to the noise power per bin is SNR_bin = (C/N0) * Tb, so C/N0 = SNR_bin / Tb. Summed over the valid blocks,
# sum(SNR_bin) = E/N0 (echo energy over noise density), and for Gaussian noise the statistic of search() is approximately
# z ~ (E/N0) / sqrt(nb), with nb the number of valid blocks.
import os, json, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
D=os.environ.get('WAVDIR','../data/raw_wav/')
TLE['ISS']=tuple(l.strip() for l in open('../data/tle/iss.txt').read().split('\n')[1:3])
def blocks(path,name,tau,nu,Tb=0.5):
    md,x,_=bwav.read(path); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs); g=GeoSat(path,name)
    fmax=0.9*min(fb,fs/2-fb); B=int(Tb*fs); nb=len(bb)//B
    fd=np.interp(ts+tau,g.t,g.fd); y=bb*np.exp(-2j*np.pi*np.cumsum(fd)/fs)
    valid=(np.abs(fd)>30)&(np.abs(fd)<fmax); vb=valid[:nb*B].reshape(nb,B).all(1)
    Y=y[:nb*B].reshape(nb,B); F=np.fft.fftshift(np.fft.fft(Y,axis=1),axes=1)/B; f=np.fft.fftshift(np.fft.fftfreq(B,1/fs))
    P=np.abs(F)**2; i0=np.argmin(abs(f-nu)); noise_cols=(abs(f)>14)&(abs(f)<40)
    noise=np.median(P[:,noise_cols],axis=1)/np.log(2)          # median of exponential -> mean
    snr=P[:,i0]/noise-1; tb=(np.arange(nb)+0.5)*Tb
    k=np.argmin(g.rr); tca=g.t[k]
    return dict(t=tb[vb],snr=snr[vb],nb=int(vb.sum()),rr=np.interp(tb[vb],g.t,g.rr),tca=float(tca),rrmin=float(g.rr.min()))
def summarize(r,Tb=0.5,label=''):
    s=r['snr']; sm=np.convolve(s,np.ones(4)/4,'same')     # 2-s running mean for the peak
    k=int(np.argmax(sm)); cn0=sm[k]/Tb; E=s.sum()
    # glint duration: span where the 2-s mean exceeds 1/4 of its peak (about -6 dB)
    on=r['t'][sm>sm[k]/4]; dur=float(on.max()-on.min()) if len(on) else 0
    return dict(label=label,cn0_peak_dBHz=float(10*np.log10(max(cn0,1e-3))),t_peak=float(r['t'][k]),rr_at_peak=float(r['rr'][k]),
        E_N0_dB=float(10*np.log10(max(E,1e-3))),nb=r['nb'],z_pred=float(E/np.sqrt(r['nb'])),duration_s=dur)
out=[]
iss=json.load(open('../results/iss_detection_table_v2.json'))
for d in iss:
    if d['z']<7 or abs(d['tau'])>0.5: continue
    p=D+f"RAD_BEDOUR_{d['pass_']}_{d['st']}_SYS001.wav"
    r=blocks(p,'ISS',d['tau'],d['df']); s=summarize(r,label=f"ISS {d['pass_']} {d['st']}"); s['z_obs']=d['z']; out.append(s); print(json.dumps(s),flush=True)
# z of the detection: SDR receivers with line removal (bb_clean.json), Ottignies without (bb2_results.json)
for st,z,nu in [('BEOUDS',152.3,-10.5),('FRMAUB',8.2,-10.5),('BEOTTI',21.9,-11.0)]:
    p=D+f'RAD_BEDOUR_20261105_0310_{st}_SYS001.wav'
    r=blocks(p,'SPACEMOBILE-008',-0.5,nu); s=summarize(r,label=f'SPACEMOBILE-008 {st}'); s['z_obs']=z; out.append(s); print(json.dumps(s),flush=True)
json.dump(out,open('../results/detectability_cn0.json','w'),indent=1)
