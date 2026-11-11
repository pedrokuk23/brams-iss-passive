# Uso: python3 synth_run.py <ESTACAO>  -> ecos sintéticos injetados em gravações reais sem ISS
exec(open('synth.py').read().split("if __name__")[0])
import sys
ST=sys.argv[1] if len(sys.argv)>1 and sys.argv[1].startswith(('BE','FR')) else None
ST=os.environ.get('ST',ST)
np.random.seed(1)
passes=[k for k in tr if k.endswith(ST) and len(tr[k])>=5]
noise={'BEOTTI':['20261103_0045','20261103_0215','20261103_2350'],'BEOUDS':['20261103_0045','20261104_0125','20261104_2310']}[ST]
alphas=[0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2]; dts=[-0.15,-0.05,0.05,0.15]; Nm=1.3e11
out=[]
for nk in noise:
    md,x,_=bwav.read(D+f'RAD_BEDOUR_{nk}_{ST}_SYS001.wav'); fs=md['samplerate']; fb,bb,ts=baseband(x,fs)
    sig=np.sqrt(np.median(np.abs(bb)**2)/np.log(2))
    for pk in passes:
        g=Geo(pk); a=np.array(tr[pk]); win=(a[:,0].min(),a[:,0].max())
        for al in alphas:
            for dtt in dts:
                xi=inject(x,fs,fb,g,dtt,Nm,al*sig,win); fb2,bb2,_=baseband(xi,fs)
                z,tau,df=search(bb2,ts,fs,fb2,g); rows=track(bb2,ts,fs,fb2,g); f=fit(rows,g)
                r=dict(st=ST,noise=nk,geo=pk,alpha=al,dt_true_ms=dtt*1e3,z=z,tau=tau,df=df,dt_est_ms=(f[0] if f else None),Nm_est=(f[1] if f else None),npts=(f[2] if f else 0))
                out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open(f'../results/synth_{ST}.json','w'))
