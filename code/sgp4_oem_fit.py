# Is the pass-order pattern of delta t an error of the SGP4 theory? (Sec. 6, Supplementary S10)
# For each archived NASA OEM file, an SGP4 element set (seven elements including B*) is fitted by least squares to the
# first three days of the ephemeris (EME2000 converted to TEME with astropy). The along-track residual of the best SGP4
# fit, expressed as a time offset, is what SGP4 cannot represent however its elements are chosen. It is compared, at
# the Phase 2 ISS passes inside each window, with the measured delta t, by order of the pass within the night, and its
# spectrum is computed. Writes results/sgp4_oem_fit.json.
import glob, json, os, numpy as np
from sgp4.api import Satrec, WGS72
from scipy.optimize import least_squares
from astropy.time import Time
from astropy.coordinates import GCRS, TEME, CartesianRepresentation, CartesianDifferential
import astropy.units as u
D='../data/'; Q=json.load(open('../results/phase2_iss_summary.json'))
TLE=[]
for f in sorted(glob.glob(D+'tle/phase2/*.txt')):
    L=open(f).read().splitlines()
    for i,l in enumerate(L):
        if l.startswith('1 25544') and i+1<len(L): TLE.append((l,L[i+1]))
if not TLE: raise SystemExit('sgp4_oem_fit: no ISS element sets in data/tle/phase2')
def readoem(fn):
    T=[];X=[];start=None
    for l in open(fn):
        if l.startswith('START_TIME') and start is None: start=l.split('=')[1].strip()
        p=l.split()
        if len(p)==7 and p[0][:2]=='20' and 'T' in p[0]: T.append(p[0]); X.append([float(v) for v in p[1:]])
    return Time(T,scale='utc'),np.array(X),Time(start,scale='utc')
out={}
for fn in sorted(glob.glob(D+'oem_archive/ISS_OEM_*.txt')):
    lab=os.path.basename(fn)[8:21]; t,X,t0=readoem(fn)
    m=(t>=t0)&(t<=t0+3*u.day); t=t[m]; X=X[m]
    rep=CartesianRepresentation(X[:,:3].T*u.km,differentials=CartesianDifferential(X[:,3:].T*u.km/u.s))
    te=GCRS(rep,obstime=t).transform_to(TEME(obstime=t))
    r=te.cartesian.xyz.to(u.km).value.T; v=te.cartesian.differentials['s'].d_xyz.to(u.km/u.s).value.T
    best=None
    for l1,l2 in TLE:
        s=Satrec.twoline2rv(l1,l2); d=abs(s.jdsatepoch+s.jdsatepochF-t0.jd)
        if best is None or d<best[0]: best=(d,s)
    s0=best[1]; ep=s0.jdsatepoch+s0.jdsatepochF-2433281.5
    x0=np.array([s0.no_kozai,s0.ecco,s0.inclo,s0.nodeo,s0.argpo,s0.mo,s0.bstar])
    def prop(x):
        s=Satrec(); s.sgp4init(WGS72,'i',25544,ep,x[6],0,0,x[1],x[4],x[2],x[5],x[0],x[3])
        e,rr,vv=s.sgp4_array(t.jd1,t.jd2); return rr
    sol=least_squares(lambda x:(prop(x)-r).ravel(),x0,x_scale=np.array([1e-6,1e-6,1e-5,1e-5,1e-3,1e-3,1e-5]),method='lm')
    dr=r-prop(sol.x); sp=np.linalg.norm(v,axis=1); dt=-np.sum(dr*v,1)/sp**2*1e3      # ms; >0: ephemeris behind SGP4
    tu=t.unix; dtu=dt-dt.mean()
    f=np.fft.rfftfreq(len(dtu),d=np.median(np.diff(tu))/86400.); P=np.abs(np.fft.rfft(dtu))**2
    k=np.argmax(P[1:])+1
    pas=[]
    for p in Q['passes']:
        if tu[0]+3600<p['t_pass']<tu[-1]-3600:
            pas.append(dict(pass_=p['pass_'],order=p['order'],fit_resid_ms=float(np.interp(p['t_pass'],tu,dt)),dt_ms=p['dt_ms']))
    out[lab]=dict(window_start=t0.isot,n_nodes=int(len(t)),along_rms_ms=float(dtu.std()),along_p2p_ms=float(dt.max()-dt.min()),
                  pos_rms_km=float(np.sqrt(np.mean(np.sum(dr**2,1)))),peak_cycles_per_day=float(f[k]),passes=pas)
    print(lab,'rms %.1f ms p2p %.0f peak %.2f c/d, %d passes'%(dtu.std(),dt.max()-dt.min(),f[k],len(pas)))
# pooled comparison over all passes covered (each pass once: the window that starts latest before it)
allp={}
for lab in sorted(out):
    for p in out[lab]['passes']: allp[p['pass_']]=p
A=np.array([p['fit_resid_ms'] for p in allp.values()]); B=np.array([p['dt_ms'] for p in allp.values()])
order={str(o):dict(n=int(sum(p['order']==o for p in allp.values())),fit_resid_mean_ms=float(np.mean([p['fit_resid_ms'] for p in allp.values() if p['order']==o])),
                   dt_mean_ms=float(np.mean([p['dt_ms'] for p in allp.values() if p['order']==o]))) for o in sorted({p['order'] for p in allp.values()})}
sl,ic=np.polyfit(A,B,1)
summ=dict(n_passes=len(allp),corr_fit_resid_vs_dt=float(np.corrcoef(A,B)[0,1]),dt_rms_ms=float(np.sqrt(np.mean((B-B.mean())**2))),
          dt_minus_fit_resid_rms_ms=float(np.sqrt(np.mean((B-A-(B-A).mean())**2))),slope_dt_on_fit_resid=float(sl),by_order=order,
          along_rms_ms_range=[min(v['along_rms_ms'] for v in out.values()),max(v['along_rms_ms'] for v in out.values())],
          peak_cycles_per_day=sorted({round(v['peak_cycles_per_day'],2) for v in out.values()}))
json.dump(dict(summary=summ,files=out),open('../results/sgp4_oem_fit.json','w'),indent=1)
print(json.dumps(summ,indent=1))
