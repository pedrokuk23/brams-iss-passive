# Cluster bootstrap of the antenna gain map with many replicates (one-step Newton, i.e. linearized, bootstrap).
# Each replicate resamples aircraft tracks with replacement (weights M_c) and moves the full-sample estimate by one
# Newton step of the weighted censored likelihood, using the exact weighted Hessian. Also reports the gradient norm at
# the full-sample optimum (convergence check) and the cluster-robust sandwich covariance as a cross-check.
import os,sys,json,time,numpy as np
from scipy.sparse import diags,csr_matrix
exec(open('antmap_tobit.py').read().split("p0=np.zeros(P)")[0])
NB=int(os.environ.get('NBOOT',1000)); rng=np.random.default_rng(int(os.environ.get('SEED',7)))
pf=np.load(f'{D}/tobit_p_T{int(TDB)}.npy')
J0,G0=f(pf); print('nll',round(J0,1),'|grad|max',float(np.abs(G0).max()),'|grad|2',float(np.linalg.norm(G0)),flush=True)
mu=A@pf; s=np.exp(pf[-1]); z=(y-mu)/s
lam=np.exp(norm.logpdf(z)-log_ndtr(z))
# first derivatives of the per-block log-likelihood
gmu=np.where(det,z/s,-lam/s); gsg=np.where(det,z**2-1,-lam*z)
# second derivatives (observed information = minus these)
hmumu=np.where(det,-1/s**2,-lam*(z+lam)/s**2)
hmusg=np.where(det,-2*z/s,(lam-lam*z*(z+lam))/s)
hsgsg=np.where(det,-2*z**2,lam*z-lam*z**2*(z+lam))
# prior Hessian (track prior and sum-to-zero penalties)
Hp=np.zeros((P,P)); Hp[np.arange(NT),np.arange(NT)]+=1/ST_T**2; Hp[:NT,P-2]-=1/ST_T**2; Hp[P-2,:NT]-=1/ST_T**2; Hp[P-2,P-2]+=NT/ST_T**2
for g in groups:
    g=np.array(list(g)); Hp[np.ix_(g,g)]+=100
def hess(w):
    H=(A.T@diags(-w*hmumu)@A).toarray()
    c=A.T@(-w*hmusg); H[:,-1]+=c; H[-1,:]+=c; H[-1,-1]+=(-w*hsgsg).sum()
    return H+Hp
# per-block score rows (gradient of +loglik wrt p): A_i*gmu_i, last entry gsg_i; aggregated by track
Ci=csr_matrix((np.ones(n),(ti,np.arange(n))),shape=(NT,n))
Sc=(Ci@diags(gmu)@A).toarray(); Sc[:,-1]=Ci@gsg      # NT x P cluster scores
H1=hess(np.ones(n)); H1i=np.linalg.inv(H1)
V=H1i@(Sc.T@Sc)@H1i*NT/(NT-1); sd_sand=np.sqrt(np.diag(V))
reps=[]; t0=time.time()
for b in range(NB):
    M=rng.multinomial(NT,np.ones(NT)/NT); w=M[ti].astype(float)
    # gradient of weighted negative objective at pf: prior gradient minus weighted data score;
    # at the optimum the prior gradient equals the unweighted data score, so G_w = -(M-1)^T Sc  (+ G0 residual)
    Gw=-(M-1)@Sc+G0
    reps.append((pf-np.linalg.solve(hess(w),Gw))[NT:])
    if b%50==0: print('rep',b,round(time.time()-t0),flush=True)
R=np.array(reps); sd=R.std(axis=0,ddof=1)
# convergence of the bootstrap SD with replicate count
conv={int(k):float(np.median(R[:k].std(axis=0,ddof=1)/sd)) for k in [50,100,200,500,NB] if k<=NB}
def pack(v):
    sl=lambda a,b:v[a-NT:b-NT].tolist()
    return dict(sd_az={STS[k]:{'low':sl(o[3]+k*24,o[3]+k*24+12),'high':sl(o[3]+k*24+12,o[3]+k*24+24)} for k in range(5)},
      sd_el={STS[k]:sl(o[4]+k*8,o[4]+k*8+8) for k in range(5)},sd_off={STS[k]:float(v[o[5]+k-NT]) for k in range(5)},
      sd_beta=sl(o[1],o[2]),sd_eltx=sl(o[2],o[3]))
out=dict(nboot=NB,unit='aircraft track',method='one-step Newton cluster bootstrap (exact weighted Hessian)',
  grad_max_at_optimum=float(np.abs(G0).max()),sd_ratio_convergence=conv,**pack(sd),sandwich=pack(sd_sand[NT:]))
json.dump(out,open('../results/antenna_map_bootstrap_1s.json','w'),indent=1)
old=json.load(open('../results/antenna_map_bootstrap.json'))
print('conv',conv)
for k in STS: print(k,'off sd 1s %.2f sand %.2f full10 %.2f'%(out['sd_off'][k],out['sandwich']['sd_off'][k],old['sd_off'][k]),
   'az-high max 1s %.1f full10 %.1f'%(max(out['sd_az'][k]['high']),max(old['sd_az'][k]['high'])))
