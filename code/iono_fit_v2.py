# Ajuste do resíduo Doppler com todos os receptores (SDR com remoção de linhas CW), correção da base de tempo (PPS)
# e comparação com a efeméride da NASA. Modelo: Δf = -ḟ_D·(δt + off_r) + K·d/dt[STEC] (Chapman, hmF2=300 km, H=55 km)
import numpy as np, json, os
exec(open('iono_fit.py').read().split("passes={")[0])
tr=json.load(open('../results/iono_tracks_v2.json')); D=os.environ.get('WAVDIR','../data/raw_wav/')
PPS={r['file'][11:31]:r['off0_ms'] for r in json.load(open('../results/pps_audit.json'))}
OEM=json.load(open('../results/oem_cmp_v2.json'))
NMIN=int(os.environ.get('NMIN',19))   # pontos mínimos por receptor
cache={}
def cols(key):
    if key in cache: return cache[key]
    md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    a=cleaned(tr[key]); tc=a[:,0]; t=np.arange(tc.min()-5,tc.max()+5,0.5); X=iss_ecef(t0+t)
    Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
    S=np.array([stec(TX,p,1.0)+stec(p,RX,1.0) for p in X]); dS=np.gradient(S,t)
    mt=-np.interp(tc,t,dfd); off=PPS[key]*1e-3
    cache[key]=(a[:,1]-mt*off,mt,K*np.interp(tc,t,dS)); return cache[key]   # observação já corrigida do erro de relógio
def fit(keys,off=None):
    # off: receiver offsets in ms (rx_offsets.py), removed from each record before a joint fit
    obs=np.concatenate([cols(k)[0]-(cols(k)[1]*(off.get(k[14:],0.0)*1e-3) if off else 0.0) for k in keys]); mt=np.concatenate([cols(k)[1] for k in keys]); mi=np.concatenate([cols(k)[2] for k in keys])
    A=np.column_stack([mt,mi]); p,*_=np.linalg.lstsq(A,obs,rcond=None); r=obs-A@p; rms=np.sqrt(np.mean(r**2))
    C=np.linalg.inv(A.T@A)*rms**2*len(obs)/max(len(obs)-2,1)
    return dict(n=int(len(obs)),dt_ms=p[0]*1e3,sdt_ms=np.sqrt(C[0,0])*1e3,Nm=p[1],rms=rms)
passes=sorted({k[:13] for k in tr}); res={}
def oem(k): 
    lab='1104' if k.startswith('20261104_2305') else '1102'; return OEM[k+'_'+lab]['dt_ms']
for p in passes:
    keys=[k for k in tr if k.startswith(p) and len(cleaned(tr[k]))>=NMIN]
    for k in keys:
        r=fit([k]); r['oem_ms']=oem(k); r['pps_ms']=PPS[k]; res[k]=r
# joint fit of each pass after removing the receiver offsets relative to Oudsbergen (Ottignies; Humain and Uccle)
from rx_offsets import fit_offsets, correction
OFF=fit_offsets([(k[:13],k[14:],v['dt_ms'],v['sdt_ms']) for k,v in res.items()],v5=('BEHUMA','BEUCCE'))
offms={rx:correction(rx,OFF,v5=('BEHUMA','BEUCCE')) for rx in ('BEOTTI','BEOUDS','BEHUMA','BEUCCE')}
for p in passes:
    keys=[k for k in tr if k.startswith(p) and len(cleaned(tr[k]))>=NMIN]
    if keys:
        r=fit(keys,offms); r['receivers']=[k[14:] for k in keys]; r['oem_ms']=float(np.mean([oem(k) for k in keys])); r['offsets_removed_ms']=offms; res[p]=r
for k,v in res.items():
    print(f"{k:22s} n={v['n']:3d} δt={v['dt_ms']:+7.1f}±{v['sdt_ms']:4.1f} ms  OEM {v['oem_ms']:+7.1f}  dif {v['dt_ms']-v['oem_ms']:+6.1f}  Nm={v['Nm']:.2e} rms={v['rms']:.2f}")
json.dump(res,open('../results/iono_fit_v2.json','w'),indent=1)
