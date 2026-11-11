# Auditoria da base de tempo: compara a escala de tempo usada (start_us + k/fs_cabeçalho) com os registros PPS/tempo do bloco BRA2
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import bwav, glob, numpy as np, json, sys
files=sorted(glob.glob(WAVDIR_+'*.wav'))
rows=[]
for f in files:
    md,x,p=bwav.read(f); p=np.array(p,dtype=float); n=len(x)
    if len(p)<10: continue
    k=p[:,0]; t=(p[:,1]-md['start_us'])/1e6
    A=np.vstack([np.ones_like(k),k]).T; (a,b),*_=np.linalg.lstsq(A,t,rcond=None); res=t-(a+b*k)
    kk=np.array([0,n-1]); used=kk/md['samplerate']; true=a+b*kk
    st=f.split('_')[-2]
    rows.append(dict(file=f.split('/')[-1],st=st,n_pps=len(p),fs_hdr=md['samplerate'],fs_pps=1/b,ppm=(1/b/md['samplerate']-1)*1e6,
                     off0_ms=-a*1e3, offend_ms=(used[1]-true[1])*1e3, jitter_rms_us=res.std()*1e6, jitter_max_us=np.abs(res).max()*1e6))
json.dump(rows,open('../results/pps_audit.json','w'),indent=0)
import collections
for st in sorted(set(r['st'] for r in rows)):
    R=[r for r in rows if r['st']==st]
    f=lambda k: (min(r[k] for r in R),np.median([r[k] for r in R]),max(r[k] for r in R))
    print(st,len(R),'files | start offset ms (min/med/max) %.3f %.3f %.3f'%f('off0_ms'),'| end offset ms %.3f %.3f %.3f'%f('offend_ms'),'| fs ppm %.2f %.2f %.2f'%f('ppm'),'| jitter rms us %.1f %.1f %.1f'%f('jitter_rms_us'),'| max %.0f'%max(r['jitter_max_us'] for r in R))
