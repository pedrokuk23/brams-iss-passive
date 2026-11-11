# Incerteza do mapa de ganho por bootstrap agrupado: reamostra trilhas de aviões (com reposição) e reajusta o modelo censurado.
# Blocos de uma mesma trilha não são independentes; a unidade de reamostragem é a trilha (avião x arquivo).
import os,sys,json,numpy as np
exec(open('antmap_tobit.py').read().split("p0=np.zeros(P)")[0])
NB=int(os.environ.get('NBOOT',20)); rng=np.random.default_rng(7)
WS=f'{D}/tobit_p_T{int(TDB)}.npy'
def fit_w(W,p0,maxit):
    def fw(p):
        mu=A@p; s=np.exp(p[-1]); m0=p[-2]; z=(y-mu)/s
        ll=np.where(det,-0.5*z**2-np.log(s),log_ndtr(z))*W
        g_mu=np.where(det,z/s,np.exp(norm.logpdf(z)-log_ndtr(z))*(-1/s))*W
        g_s=np.where(det,z**2-1,np.exp(norm.logpdf(z)-log_ndtr(z))*(-z))*W
        J=-ll.sum(); G=-(A.T@g_mu); G[-1]=-g_s.sum()
        tr=p[:NT]-m0; J+=0.5*(tr**2).sum()/ST_T**2; G[:NT]+=tr/ST_T**2; G[-2]+=-(tr).sum()/ST_T**2
        for g in groups:
            g=list(g); sm=p[g].sum(); J+=50*sm**2; G[g]+=100*sm
        return J,G
    return minimize(fw,p0,jac=True,method='L-BFGS-B',options=dict(maxiter=maxit,maxfun=4*maxit,gtol=1e-6,ftol=1e-12)).x
if os.path.exists(WS) and len(np.load(WS))==P: pfull=np.load(WS)
else:
    pfull=np.zeros(P); pfull[:NT]=np.median(y[det]); pfull[-2]=pfull[0]; pfull[-1]=np.log(4)
    pfull=fit_w(np.ones(n),pfull,12000); np.save(WS,pfull)
reps=[]
for b in range(NB):
    M=rng.multinomial(NT,np.ones(NT)/NT); W=M[ti].astype(float)
    p=fit_w(W,pfull,int(os.environ.get('MAXIT',1500))); reps.append(p[NT:]); print('rep',b,flush=True)
R=np.array(reps); sd=R.std(axis=0,ddof=1); full=pfull[NT:]
def sl(a,b): return sd[a-NT:b-NT].tolist()
out=dict(nboot=NB,unit='aircraft track',
  sd_az={STS[k]:{'low':sl(o[3]+k*24,o[3]+k*24+12),'high':sl(o[3]+k*24+12,o[3]+k*24+24)} for k in range(5)},
  sd_el={STS[k]:sl(o[4]+k*8,o[4]+k*8+8) for k in range(5)},
  sd_off={STS[k]:float(sd[o[5]+k-NT]) for k in range(5)},
  sd_beta=sl(o[1],o[2]),sd_eltx=sl(o[2],o[3]))
json.dump(out,open('../results/antenna_map_bootstrap.json','w'),indent=1)
np.set_printoptions(precision=1,suppress=True,linewidth=160)
for k in STS: print(k,'sd az>10',np.array(out['sd_az'][k]['high']),'sd el',np.array(out['sd_el'][k]),'sd off %.2f'%out['sd_off'][k])
