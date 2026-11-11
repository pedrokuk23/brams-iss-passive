# Test of brams_feasibility.py against the blind searches. For every searched pass of Phase 1 and Phase 2 the margin
# M = sigma_max / sigma_min is recomputed with the feasibility model (energy of the sinc^2 glint lobe for flat plates,
# energy-equivalent cross section for the ISS, n_b = 2 T blocks and the Phase 1 threshold z_th of the receiver), and the
# detection probability P_d(M) is fitted on Phase 1 only, as a logistic function of log10 M for M >= 1 (false-alarm rate
# below), for three groups: 'structure' (ISS complex, one unit per recording), 'plate_shadow' and 'plate_sunlit'
# (Starlink, AST: Earth-facing flat arrays, satellite in the Earth's shadow or in sunlight at closest approach). The
# grouping by sunlight was introduced after the effect was seen in Phase 2 (PROTOCOL item 8), so the prediction is not
# blind in that respect; 'no_sunlight_split' repeats it with one group for all flat arrays. Receiver gain corrections:
# Humain and Uccle were not calibrated on SPACEMOBILE-008 and take the gains of Oudsbergen and Ottignies; a correction
# per receiver (1-dB steps, Oudsbergen fixed) is then fitted on Phase 1 by maximum likelihood, alternating with the
# curve, and the prediction is also given without these corrections ('no_gain_correction').
import sys, json, glob, collections, numpy as np
from scipy.stats import beta, poisson
import brams_feasibility as BF
exec(open('eclipse.py').read().split("out={}")[0])          # sunlit(), picker(), sets_from(), jd_of()
def cp(k,n,a=0.05): return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)), 1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
AREA={'Starlink V2-mini':116,'Starlink DTC':116,'Starlink V1.x':29}
def area(r):
    if r['class'] in AREA: return AREA[r['class']]
    if r['class'].startswith('AST'): return 64 if (r['obj'].startswith('BLUEWALKER') or r['obj'][-3:] in ('001','002','003','004','005')) else 223
    return None
MODEL='zenith'      # 'zenith': gains of the SPACEMOBILE-008 calibration (paper); 'elev': with the elevation dependence
import os
DG={}               # receiver gain corrections (dB) fitted on Phase 1
ZTH={k:v['z_th'] for k,v in json.load(open('../results/blind_summary.json'))['thresholds'].items()}   # Phase 1 thresholds (used in both phases)
def smin(r,Te):
    st=r['file'][-6:]; T=r['n_valid_s']
    s0=BF.sigma_min(r['rrmin'],Te,st,T_arc=T,z_th=ZTH[st]) if MODEL=='zenith' else BF.sigma_min(r['rrmin'],Te,st,r['el_rx'],r['el_tx'],T_arc=T,z_th=ZTH[st])
    return s0/10**(DG.get(st,0)/10)
SISS=BF.TH['iss_sigma_energy_50s_m2']['median']     # energy-equivalent cross section of the ISS echoes over a 50-s arc
def margin(r):
    A=area(r)
    if A is None: return SISS*min(r['n_valid_s'],50)/50/smin(r,min(r['n_valid_s'],50))   # ISS: energy as measured
    Te=BF.glint_energy_s(r['n_valid_s'],BF.glint_s(r['rrmin'],np.sqrt(A)))
    return BF.plate_sigma(A)/smin(r,Te)
def geom(s,f):
    jd0=jd_of(f); best=None; rx=RX[f[-6:]]
    for dt in np.arange(0,300,5):
        jd=jd0+dt/86400; e,p,v=s.sgp4(int(jd),jd-int(jd))
        if e: continue
        p=np.array(p); gm=np.radians((280.46061837+360.98564736629*(jd-2451545))%360)
        pe=np.array([[np.cos(gm),np.sin(gm),0],[-np.sin(gm),np.cos(gm),0],[0,0,1]])@p
        rr=np.linalg.norm(pe-TX)*np.linalg.norm(pe-rx)
        if best is None or rr<best[0]: best=(rr,p,pe,jd)
    _,p,pe,jd=best; su=sunvec(jd); x=p@su
    return bool(x>0 or np.linalg.norm(p-x*su)>6371.0), BF.elev(f[-6:],pe), BF.elev('TX',pe)
# --- the two phases
g2={}; exec(open('phase2_analyze.py').read().split("# echo track of an owner")[0],g2); P2=g2['P']
D2={(d['file'],d['obj']) for d in json.load(open('../results/phase2_summary.json'))['frozen_phase1_thresholds']['detections']}   # thresholds of Phase 1, so that both phases use the same z_th
_tle2=sorted(glob.glob('../data/tle/phase2/*.txt'))
if not _tle2: sys.exit('feasibility_validate: no Phase 2 element sets in data/tle/phase2 (unzip auxiliary_data_phase2.zip)')
sat2=picker(sets_from(_tle2))
g1={}; exec(open('blind_analyze.py').read().split("# echo track of an owner")[0],g1); P1=g1['P']
D1={(d['file'],d['obj']) for d in json.load(open('../results/blind_summary.json'))['detections']}
sat1=picker(sets_from(['../data/tle/leo.txt']+sorted(glob.glob('../data/tle/*.txt'))))
CL=('ISS complex','AST SpaceMobile','Starlink DTC','Starlink V2-mini','Starlink V1.x')
def units(P,D,sat):
    U=[]; iss={}
    for r in P:
        if r['class'] not in CL: continue
        det=(r['file'],r['obj']) in D
        if r['class']=='ISS complex' and not r['obj'].startswith('ISS'):
            iss.setdefault(r['file'],dict(cls='ISS complex',grp='structure',file=r['file'],det=False,M=None))['det']|=det; continue
        try: lit,r['el_rx'],r['el_tx']=geom(sat(r['obj'],r['file']),r['file'])
        except KeyError: continue
        if r['class']=='ISS complex':
            u=iss.setdefault(r['file'],dict(cls='ISS complex',grp='structure',file=r['file'],det=False,M=None))
            u['det']|=det
            u['r']=r
            continue
        U.append(dict(cls=r['class'],grp='plate_sunlit' if lit else 'plate_shadow',file=r['file'],det=det,r=r))
    U+=[u for u in iss.values() if 'r' in u]
    return U
U1=units(P1,D1,sat1); U2=units(P2,D2,sat2)
assert U1 and U2, 'feasibility_validate: no passes could be evaluated (element sets or results missing)'
EDGES=[0,1,3,10,30,100,1e9]; GRP=('structure','plate_shadow','plate_sunlit'); STS=('BEOUDS','BEHUMA','BEOTTI','FRMAUB','BEUCCE')
def setM(U):
    for u in U: u['M']=margin(u['r'])
from scipy.optimize import minimize
def fit_curve(U,groups=None):
    """P_d(M) per group: false-alarm rate (1e-3) for M < 1; for M >= 1 a logistic function of log10 M fitted by maximum
    likelihood to the units of the group (two parameters, weak prior on the slope). Binned counts are kept for the figures."""
    C={'edges':EDGES[:-1],'groups':{},'model':'logistic'}
    for g in (groups or GRP):
        q=[u for u in U if u['grp']==g and u['M']>=1]; x=np.log10([u['M'] for u in q]); y=np.array([u['det'] for u in q],float)
        def nll(p):     # weak Gaussian prior on the slope (mean 3, s.d. 3 per decade of M), which matters only for the
                        # 22 ISS units of Phase 1, all but one detected
            P=np.clip(1/(1+np.exp(-(p[0]+p[1]*x))),1e-6,1-1e-6); return -np.sum(y*np.log(P)+(1-y)*np.log(1-P))+(p[1]-3)**2/18
        p=minimize(nll,[-3.,1.],method='Nelder-Mead',options=dict(xatol=1e-6,fatol=1e-8,maxiter=4000)).x
        rows=[]
        for a,b in zip(EDGES[:-1],EDGES[1:]):
            qq=[u for u in U if u['grp']==g and a<=u['M']<b]; rows.append(dict(lo=a,hi=b,n=len(qq),k=int(sum(u['det'] for u in qq))))
        C['groups'][g]=dict(a=float(p[0]),b=float(p[1]),phase1=rows,n=len(q))
    return C
def pd_of(C,u,g=None): return BF.curve_pd(C,g or u['grp'],u['M'])
def loglik(C,U): return sum(np.log(p) if u['det'] else np.log1p(-min(p,1-1e-9)) for u in U for p in [pd_of(C,u)])
# 1. receiver gain corrections on Phase 1 (Oudsbergen, the calibration receiver, fixed at 0 dB), alternating with the curve
setM(U1); C=fit_curve(U1)
for it in range(3):
    for st in STS[1:]:
        best=None
        for d in np.arange(-12,6.1,1.0):
            DG[st]=d; Us=[u for u in U1 if u['file'].endswith(st)]; setM(Us); L=loglik(C,Us)
            if best is None or L>best[0]: best=(L,d)
        DG[st]=best[1]
    setM(U1); C=fit_curve(U1)
print('receiver gain corrections fitted on Phase 1 (dB):',DG)
json.dump(dict(C,gain_correction_dB=DG,gain_model=MODEL,z_th_phase1=ZTH),open('../results/feasibility_curve.json','w'),indent=1)
# 2. prediction of Phase 2
setM(U2)
def cmp(sel):
    k=sum(u['det'] for u in sel); e=sum(pd_of(C,u) for u in sel)
    return dict(n=len(sel),observed=int(k),predicted=round(e,1),pred_95=(int(poisson.ppf(0.025,e)),int(poisson.ppf(0.975,e))))
res=dict(gain_correction_dB=DG,by_class={c:cmp([u for u in U2 if u['cls']==c]) for c in CL},by_rx={st:cmp([u for u in U2 if u['file'].endswith(st)]) for st in STS},
    by_bin={g:[dict(lo=a,hi=b,**cmp([u for u in U2 if u['grp']==g and a<=u['M']<b])) for a,b in zip(EDGES[:-1],EDGES[1:])] for g in GRP},all=cmp(U2),
    impossible={ph:dict(n=len(q),k=sum(u['det'] for u in q),ci=cp(sum(u['det'] for u in q),len(q)) if q else None) for ph,U in (('phase1',U1),('phase2',U2)) for q in [[u for u in U if u['M']<1]]})
res['by_class_plates_only']=cmp([u for u in U2 if u['grp']!='structure'])
res['table']=[dict(cls=c,rx=st,**cmp([u for u in U2 if u['cls']==c and u['file'].endswith(st)])) for c in CL for st in STS if any(u['cls']==c and u['file'].endswith(st) for u in U2)]
# variants: (a) without the sunlight split (one curve for all flat arrays), (b) without the receiver gain corrections
for u in U1+U2: u['grp0']='structure' if u['grp']=='structure' else 'plate'
C0=fit_curve([dict(u,grp=u['grp0']) for u in U1],groups=('structure','plate'))
res['no_sunlight_split']=dict(observed=int(sum(u['det'] for u in U2)),predicted=round(float(sum(BF.curve_pd(C0,u['grp0'],u['M']) for u in U2)),1))
DG0=dict(DG); DG.clear(); setM(U1); Cn=fit_curve(U1); setM(U2)
res['no_gain_correction']=dict(observed=int(sum(u['det'] for u in U2)),predicted=round(float(sum(pd_of(Cn,u) for u in U2)),1),
    by_rx={st:dict(observed=int(sum(u['det'] for u in U2 if u['file'].endswith(st))),predicted=round(float(sum(pd_of(Cn,u) for u in U2 if u['file'].endswith(st))),1)) for st in STS})
DG.update(DG0); setM(U1); setM(U2)
json.dump(res,open('../results/feasibility_validation.json','w'),indent=1)
print('Phase 1 curve:'); [print('  %-13s'%g,' '.join('%g:%d/%d'%(r['lo'],r['k'],r['n']) for r in v['phase1']),' a=%.2f b=%.2f'%(v['a'],v['b'])) for g,v in C['groups'].items()]
print('  no sunlight split',res['no_sunlight_split'],' no gain correction',res['no_gain_correction'])
for key in ('by_class','by_rx'):
    for c,v in res[key].items(): print('  %-17s n=%5d obs %4d pred %6.1f [%d,%d]'%(c,v['n'],v['observed'],v['predicted'],*v['pred_95']))
for g,rows in res['by_bin'].items():
    for v in rows:
        if v['n']: print('  %-13s M %4g-%-5g n=%5d obs %4d pred %6.1f'%(g,v['lo'],v['hi'],v['n'],v['observed'],v['predicted']))
print('  impossible',res['impossible']); print('  all',res['all'],' plates',res['by_class_plates_only'])
