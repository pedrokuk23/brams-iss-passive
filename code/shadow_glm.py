# Sunlight effect of the Starlink V2-mini passes of Phase 2 in one regression (Sec. 10.3): logistic model of detection on
# shadow/sunlight with covariates for receiver, margin M (bins), receiver elevation (bins), orbital altitude (bins) and
# 3-h UTC window (altitude: below or above 450 km), with standard errors clustered by satellite, by five-minute slot and by
# both (two-way), and a ridge-penalized refit as a check against separation.
# Writes results/shadow_glm.json.
import json, numpy as np, pandas as pd
import statsmodels.formula.api as smf
exec(open('stats_cluster.py').read().split("out=dict(n_boot")[0])       # U2m, sat2
RE=6378.137; MU=398600.4418
def alt_km(s): n=s.no_kozai/60.0; return (MU/n**2)**(1/3)-RE
rows=[]
for u in U2m:
    if u['cls']!='Starlink V2-mini': continue
    try: a=alt_km(sat2(u['r']['obj'],u['file']))
    except KeyError: continue
    M=u.get('M') or 0.0; el=u['r'].get('el_rx') or 0.0
    rows.append(dict(det=int(bool(u['det'])),shadow=int(u['grp']=='plate_shadow'),rx=u['file'][-6:],sat=u['r']['obj'],slot=u['file'][:13],
                     Mb=int(np.digitize(M,[1,3,10,30])),elb=int(np.digitize(el,[30,45,60,75])),altb=int(np.digitize(a,[450])),
                     utc=int(u['file'][9:11])//3))
D=pd.DataFrame(rows)
out=dict(n=len(D),n_sat=int(D.sat.nunique()),n_slot=int(D.slot.nunique()))
f='det ~ shadow + C(rx) + C(Mb) + C(elb) + C(altb) + C(utc)'
for name,grp in (('satellite','sat'),('recording','slot')):
    m=smf.logit(f,D).fit(disp=0,maxiter=500,method='bfgs',cov_type='cluster',cov_kwds={'groups':pd.factorize(D[grp])[0]})
    b=m.params['shadow']; se=m.bse['shadow']
    out['cluster_'+name]=dict(or_shadow=float(np.exp(b)),ci95=[float(np.exp(b-1.96*se)),float(np.exp(b+1.96*se))],p=float(m.pvalues['shadow']))
# two-way clustering by satellite and slot (Cameron, Gelbach and Miller): V = V_sat + V_slot - V_sat x slot
V={}
for name,g in (('sat',D['sat']),('slot',D['slot']),('both',D['sat']+'|'+D['slot'])):
    m=smf.logit(f,D).fit(disp=0,maxiter=500,method='bfgs',cov_type='cluster',cov_kwds={'groups':pd.factorize(g)[0]})
    V[name]=m.cov_params().loc['shadow','shadow']; b2=m.params['shadow']
se2=float(np.sqrt(max(V['sat']+V['slot']-V['both'],max(V['sat'],V['slot']))))
out['cluster_two_way']=dict(or_shadow=float(np.exp(b2)),ci95=[float(np.exp(b2-1.96*se2)),float(np.exp(b2+1.96*se2))])
# separation check: ridge-penalized fit (small L2 penalty on all coefficients except the intercept)
X=pd.get_dummies(D[['rx','Mb','elb','altb','utc']].astype(str),drop_first=True).astype(float); X.insert(0,'shadow',D['shadow'].astype(float)); X.insert(0,'const',1.0)
from scipy.optimize import minimize
Xv=X.values; y=D['det'].values.astype(float); lam=1.0
def nll(w): z=Xv@w; return np.sum(np.logaddexp(0,z)-y*z)+0.5*lam*np.sum(w[1:]**2)
def grad(w): z=Xv@w; g=Xv.T@(1/(1+np.exp(-z))-y); g[1:]+=lam*w[1:]; return g
w=minimize(nll,np.zeros(Xv.shape[1]),jac=grad,method='L-BFGS-B').x
out['ridge_lambda1']=dict(or_shadow=float(np.exp(w[1])))
# interaction with altitude: is the effect present below 450 km?
m=smf.logit('det ~ shadow*C(altb) + C(rx) + C(Mb) + C(elb) + C(utc)',D).fit(disp=0,maxiter=500,method='bfgs',cov_type='cluster',cov_kwds={'groups':pd.factorize(D['sat'])[0]})
out['shadow_by_altitude_bin']={}
for k in sorted(D.altb.unique()):
    nm='shadow' if k==0 else 'shadow:C(altb)[T.%d]'%k
    if k==0: b=m.params['shadow']; v=m.cov_params().loc['shadow','shadow']
    else:
        if nm not in m.params: continue
        b=m.params['shadow']+m.params[nm]; C_=m.cov_params(); v=C_.loc['shadow','shadow']+C_.loc[nm,nm]+2*C_.loc['shadow',nm]
    out['shadow_by_altitude_bin'][str(k)]=dict(n=int((D.altb==k).sum()),or_shadow=float(np.exp(b)),ci95=[float(np.exp(b-1.96*np.sqrt(v))),float(np.exp(b+1.96*np.sqrt(v)))])
json.dump(out,open('../results/shadow_glm.json','w'),indent=1); print(json.dumps(out,indent=1))
