# Ecos sintéticos em gravações reais sem ISS, com remoção de linhas CW (receptores SDR). Uso: ST=BEHUMA python3 synth_v2.py
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
exec(open('synth.py').read().split("if __name__")[0])
from linefilt import maybe_clean
ST=os.environ['ST']; tr2=json.load(open('../results/iono_tracks_v2.json')); np.random.seed(3)
passes=[k for k in tr2 if k.endswith(ST) and len(cleaned(tr2[k]))>=20]
noise={'BEHUMA':['20261104_2300','20261104_2310','20261104_0125'],'BEOUDS':['20261103_0045','20261104_0125','20261104_2310']}[ST]
alphas=[0.03,0.04,0.05,0.06,0.08,0.1,0.15]; dts=[-0.1,0.1]; Nm=1.3e11; out=[]
for nk in noise:
    md,x,_=bwav.read(D+f'RAD_BEDOUR_{nk}_{ST}_SYS001.wav'); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bbc=maybe_clean(bb,ts,fs)
    sig=np.sqrt(np.median(np.abs(bbc)**2)/np.log(2))          # ruído após a remoção das linhas
    for pk in passes:
        g=Geo(pk); a=np.array(tr2[pk]); win=(a[:,0].min(),a[:,0].max())
        for al in alphas:
            for dtt in dts:
                xi=inject(x,fs,fb,g,dtt,Nm,al*sig,win); fb2,bb2,_=baseband(xi,fs); bb2=maybe_clean(bb2,ts,fs)
                z,tau,df=search(bb2,ts,fs,fb2,g); rows=track(bb2,ts,fs,fb2,g); f=fit(rows,g)
                r=dict(st=ST,noise=nk,geo=pk,alpha=al,dt_true_ms=dtt*1e3,z=z,tau=tau,df=df,dt_est_ms=(f[0] if f else None),Nm_est=(f[1] if f else None),npts=(f[2] if f else 0))
                out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open(f'../results/synth_v2_{ST}.json','w'))
