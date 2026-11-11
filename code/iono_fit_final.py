# Ajuste final do resíduo Doppler (Tabela IV): por passagem (todos os receptores) e por receptor.
# Modelo: Δf = -ḟ_D·δt + K·d/dt[STEC_TS+STEC_SR] com camada de Chapman (hmF2=300 km, H=55 km).
import numpy as np, json
exec(open('iono_fit.py').read().split("passes={")[0])
passes={'0220 03/11':['20261103_0220_BEOTTI','20261103_0220_BEOUDS'],'2355 03/11':['20261103_2355_BEOTTI','20261103_2355_BEOUDS'],
        '0040 03/11':['20261103_0040_BEOTTI'],'0130 04/11':['20261104_0130_BEOTTI','20261104_0130_BEOUDS'],
        '2305 04/11':['20261104_2305_BEOTTI','20261104_2305_BEOUDS']}
cache={}
def cols(key):
    if key in cache: return cache[key]
    md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    a=cleaned(tr[key]); tc=a[:,0]; t=np.arange(tc.min()-5,tc.max()+5,0.5); X=iss_ecef(t0+t)
    Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
    S=np.array([stec(TX,p,1.0)+stec(p,RX,1.0) for p in X]); dS=np.gradient(S,t)
    cache[key]=(a[:,1],-np.interp(tc,t,dfd),K*np.interp(tc,t,dS)); return cache[key]
def fit(keys):
    obs=np.concatenate([cols(k)[0] for k in keys]); mt=np.concatenate([cols(k)[1] for k in keys]); mi=np.concatenate([cols(k)[2] for k in keys])
    def lsq(c):
        A=np.column_stack(c); p,*_=np.linalg.lstsq(A,obs,rcond=None); return p,np.sqrt(np.mean((obs-A@p)**2))
    _,rc=lsq([np.ones_like(obs)]); _,rt=lsq([mt]); p,r=lsq([mt,mi])
    return dict(n=len(obs),dt_ms=p[0]*1e3,Nm=p[1],rms=r,rms_const=rc,rms_time=rt)
res={}
for pn,keys in passes.items():
    keys=[k for k in keys if k in tr and len(tr[k])>=5]
    res[pn]=fit(keys)
    for k in keys: res[pn+'|'+k.split('_')[2]]=fit([k])
for k,v in res.items(): print(k, {a:round(b,3) for a,b in v.items()})
json.dump(res,open('../results/iono_fit_final.json','w'),indent=1)
