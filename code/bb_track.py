# Rastreio do eco do SPACEMOBILE-008 (05/11/2026 03:10) e ajuste do termo de tempo (δt); correção da base de tempo em pps_audit.json
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
exec(open('bb2.py').read())
import json
A={r['file']:r for r in json.load(open('../results/pps_audit.json'))}
res={}
for st in ['BEOTTI','BEOUDS','FRMAUB']:
    f=f'{WAVDIR_}RAD_BEDOUR_20261105_0310_{st}_SYS001.wav'
    md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); g=GeoSat(f,'SPACEMOBILE-008'); g.fdot=np.gradient(g.fd,g.t)
    rows=track(bb,ts,fs,fb,g); a=rows[(rows[:,0]>35)&(rows[:,0]<85)]
    a=a[np.abs(a[:,1]-np.median(a[:,1]))<4]; tc=a[:,0]; A_=-np.interp(tc,g.t,g.fdot)
    p=np.dot(A_,a[:,1])/np.dot(A_,A_); r=a[:,1]-A_*p; off=A.get(f.split('/')[-1],{}).get('off0_ms',np.nan)
    res[st]=dict(n=len(a),df_med=float(np.median(a[:,1])),dt_ms=float(p*1e3),dt_pps_ms=float(p*1e3-off),rms=float(np.sqrt(np.mean(r**2))),pps_off_ms=off)
    print(st,res[st])
json.dump(res,open('../results/bluebird_0310.json','w'),indent=1)
