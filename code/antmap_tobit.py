# Regressão censurada (Tobit) para o padrão das antenas receptoras com aviões.
# y = 10log10(SNR-1) - (-20log10(R1R2))  observado se SNR>T; senão só se sabe y < T + 20log(R1R2)
# y = trilha + f(beta) + g(el_tx) + estação×az(30°)×classe de elevação(<10°, >=10°) + estação×el + ruído N(0,s)
import json,numpy as np,sys
from scipy.sparse import coo_matrix
from scipy.optimize import minimize
from scipy.special import log_ndtr
from scipy.stats import norm
D=sys.argv[1] if len(sys.argv)>1 else '../results/antmap'; TDB=float(sys.argv[2]) if len(sys.argv)>2 else 6
d=json.load(open(f'{D}/blocks.json')); C=d['cols']; ci={k:i for i,k in enumerate(C)}
R=[r for r in d['rows'] if r[-1]==0 and r[ci['el_rx']]>0]
STS=['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']
col=lambda k:np.array([r[ci[k]] for r in R],float)
snr=col('snr'); PL=20*np.log10(col('r1')*col('r2'))
det=10*np.log10(np.maximum(snr-1,1e-6))>TDB
y=np.where(det,10*np.log10(np.maximum(snr-1,1e-6)),TDB)+PL      # valor ou limite superior
st=np.array([STS.index(r[0]) for r in R]); trk=[(r[2],r[1]) for r in R]
tid={t:i for i,t in enumerate(sorted(set(trk)))}; ti=np.array([tid[t] for t in trk])
beta=col('beta'); eltx=col('el_tx'); az=col('az_rx'); el=col('el_rx')
EB=np.array([0,3,6,10,15,20,30,45,90.])
bb=np.clip((beta//10).astype(int),0,17); gb=np.clip(np.digitize(eltx,[0,5,10,15,20,30,45,90])-1,0,6)
ab=(az//30).astype(int)%12; ec=(el>=10).astype(int); eb=np.clip(np.digitize(el,EB)-1,0,len(EB)-2)
NT=len(tid); o=[0,NT,NT+18,NT+25,NT+25+120,NT+25+120+40,NT+25+120+40+5]
idx=np.stack([ti,o[1]+bb,o[2]+gb,o[3]+st*24+ec*12+ab,o[4]+st*8+eb,o[5]+st],1)
n=len(y); P=o[6]+2   # +mu0 (média das trilhas) + log s
A=coo_matrix((np.ones(6*n),(np.repeat(np.arange(n),6),idx.ravel())),shape=(n,P)).tocsr()
groups=[range(o[1],o[2]),range(o[2],o[3])]+[range(o[3]+s*12,o[3]+s*12+12) for s in range(10)]+[range(o[4]+s*8,o[4]+s*8+8) for s in range(5)]+[range(o[5],o[5]+5)]
ST_T=6.0
def f(p):
    mu=A@p; s=np.exp(p[-1]); m0=p[-2]; z=(y-mu)/s
    ll=np.where(det,-0.5*z**2-np.log(s),log_ndtr(z))
    g_mu=np.where(det,z/s,np.exp(norm.logpdf(z)-log_ndtr(z))*(-1/s))   # d ll/d mu
    g_s=np.where(det,z**2-1,np.exp(norm.logpdf(z)-log_ndtr(z))*(-z))   # d ll/d log s
    J=-ll.sum(); G=-(A.T@g_mu); G[-1]=-g_s.sum()
    # priori das trilhas N(m0,ST_T) e vínculo soma-zero nos grupos
    tr=p[:NT]-m0; J+=0.5*(tr**2).sum()/ST_T**2; G[:NT]+=tr/ST_T**2; G[-2]+=-(tr).sum()/ST_T**2
    for g in groups:
        g=list(g); sm=p[g].sum(); J+=50*sm**2; G[g]+=100*sm
    return J,G
p0=np.zeros(P); p0[:NT]=np.median(y[det]) if det.any() else 0; p0[-2]=p0[0]; p0[-1]=np.log(4)
import os
WS=f'{D}/tobit_p_T{int(TDB)}.npy'
if os.path.exists(WS) and len(np.load(WS))==P: p0=np.load(WS)
r=minimize(f,p0,jac=True,method='L-BFGS-B',options=dict(maxiter=int(os.environ.get('MAXIT',3000)),maxfun=20000,gtol=1e-6,ftol=1e-12))
np.save(WS,r.x)
p=r.x; print(r.message,'nll',round(r.fun,1),'s',round(np.exp(p[-1]),2),'n',n,'det',det.sum(),'trilhas',NT)
# incerteza aproximada por bloco: s/sqrt(n_eff) (somente indicativa)
cnt=lambda ii,m:np.bincount(ii,minlength=m)
outd=dict(T_db=TDB,n=int(n),ndet=int(det.sum()),ntrk=NT,s=float(np.exp(p[-1])),
  beta=p[o[1]:o[2]].tolist(),beta_n=cnt(bb,18).tolist(),eltx=p[o[2]:o[3]].tolist(),eltx_n=cnt(gb,7).tolist(),el_edges=EB.tolist(),
  az={STS[k]:{'low':p[o[3]+k*24:o[3]+k*24+12].tolist(),'high':p[o[3]+k*24+12:o[3]+k*24+24].tolist(),
       'n_low':cnt(ab[(st==k)&(ec==0)],12).tolist(),'n_high':cnt(ab[(st==k)&(ec==1)],12).tolist(),
       'ndet_low':cnt(ab[(st==k)&(ec==0)&det],12).tolist(),'ndet_high':cnt(ab[(st==k)&(ec==1)&det],12).tolist()} for k in range(5)},
  off={STS[k]:float(p[o[5]+k]) for k in range(5)},
  el={STS[k]:{'x':p[o[4]+k*8:o[4]+k*8+8].tolist(),'n':cnt(eb[st==k],8).tolist(),'ndet':cnt(eb[(st==k)&det],8).tolist()} for k in range(5)})
json.dump(outd,open('../results/antenna_map_tobit.json','w'))
np.set_printoptions(precision=1,suppress=True,linewidth=160)
print('beta',np.array(outd['beta'])); print('eltx',np.array(outd['eltx']))
for k in STS:
    a=outd['az'][k]; print(k,'az<10',np.array(a['low']),a['ndet_low']); print(k,'az>10',np.array(a['high']),a['ndet_high'])
print('offset',outd['off'])
for k in STS: print(k,'el',np.array(outd['el'][k]['x']),outd['el'][k]['ndet'])
