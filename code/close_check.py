# For each candidate of the close-pass scan: time of the echo in the map dechirped along the candidate's Doppler, and the
# frequency difference at that time between the candidate's predicted echo and the predicted ISS echo in the same file.
# A difference of a few hertz means the candidate's search window contains the ISS echo (ambiguity); a large one excludes it.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os, json, numpy as np
exec(open('bb2.py').read())
from linefilt import maybe_clean
exec(open('glint.py').read().split("def glint_stat")[0].split("WS=[4,8,16,32,64]")[1])
D=os.environ.get('WAVDIR','../data/raw_wav/')
TLE['ISS']=tuple(l.strip() for l in open('../data/tle/iss.txt').read().split('\n')[1:3])
cs=json.load(open('../results/close_scan.json'))
cand=[o for o in cs if o['case']=='pass' and o['z']>=7]
out=[]
for o in cand:
    p=D+f"RAD_BEDOUR_{o['file']}_SYS001.wav"; md,x,_=bwav.read(p); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
    g=GeoSat(p,o['obj']); gi=GeoSat(p,'ISS')
    s,vb=snr_map(bb,ts,fs,fb,g,o['tau'])
    j=int(np.argmin(abs(np.arange(-12,12.01,0.5)-o['df']))); col=np.nan_to_num(s[:,j]-1)
    sm=np.convolve(col,np.ones(4),'same'); k=int(np.argmax(sm)); tk=(k+0.5)*0.5
    f_c=np.interp(tk+o['tau'],g.t,g.fd)+o['df']; f_iss=np.interp(tk,gi.t,gi.fd)
    # same test for the ISS echo itself: power of the ISS-dechirped map at that time near nu=0
    r=dict(obj=o['obj'],file=o['file'],z=o['z'],tau=o['tau'],nu=o['df'],t_echo=tk,f_cand=float(f_c),f_iss=float(f_iss),
           df_iss=float(f_c-f_iss),rr_cand=float(np.interp(tk,g.t,g.rr)),rr_iss=float(np.interp(tk,gi.t,gi.rr)))
    out.append(r); print(json.dumps({k2:(round(v,1) if isinstance(v,float) else v) for k2,v in r.items()}),flush=True)
json.dump(out,open('../results/close_check.json','w'),indent=1)
