"""Can the BRAMS network detect a given satellite?  Feasibility from the radar equation, the object's shape and attitude.

    sigma_min = z_th sqrt(n_b) / (T_eff * C/N0(sigma = 1 m^2))               [m^2]   (Eq. sigmamin of the paper)
    n_b       = 2 T  (0.5-s blocks over the usable arc T, |f_D| in 30-900 Hz)
    C/N0(1)   = Pt Gt Gr lambda^2 / ((4 pi)^3 R^4 k T_sky),  R = sqrt(R1 R2)  [Hz per m^2]
    T_eff     = echo energy / peak power: the usable arc (at most 50 s) for a steady echo; for a glint of full width T_g
                (first null to first null) the integral of its sinc^2 lobe over the arc, which tends to T_g/2
    M         = sigma_max / sigma_min                                         margin; M < 1: impossible in any attitude

Pocket form (Oudsbergen gain, z_th = 7):  sigma_min ~ 100 m^2 (R / 400 km)^4 (50 s / T_eff) (T / 50 s)^(1/2) F_rx.

Shapes (lambda = 6.0 m): flat plate of area A seen face-on 4 pi A^2/lambda^2, glint T_g = 2 theta_0 / Omega with
theta_0 = lambda / 2L and Omega = v R_E / [R (R_E + R)] (bisector rate of a nadir-pointing facet near the zenith, L the
plate side, sqrt(A) if not given), or lambda P/(2 pi L) for an object tumbling with period P; large complex structure
(ISS): sigma_max given as an energy-equivalent cross section over the arc, no glint limit; resonant wire 0.86 lambda^2;
sphere of radius a << lambda: 9 pi a^2 (k a)^4.

Library:  passes(tle_lines, t0, t1, ...) -> list of predicted passes with sigma_min, sigma_max, M and sunlit flag
          expected(passes, group) -> expected detections with the detection probability curve P_d(M) of results/feasibility_curve.json
CLI:      python3 brams_feasibility.py TLE_FILE NAME --start 2026-11-10 --days 7 --plate-area 116 [--side 10] [--rx BEOUDS,BEHUMA]
"""
import json, os, argparse, datetime as dtm, numpy as np
from sgp4.api import Satrec
from sgp4.propagation import gstime
HERE=os.path.dirname(os.path.abspath(__file__))
TH=json.load(open(os.path.join(HERE,'../results/detect_theory.json')))
LAM=TH['lam_m']; KB=1.380649e-23; Z_TH=7.0             # threshold of the ISS search; the blind search uses each receiver's own
GDB=dict(TH['Geff_dB']); GDB.setdefault('BEHUMA',GDB['BEOUDS']); GDB.setdefault('BEUCCE',GDB['BEOTTI'])   # assumed, not calibrated
GEO=json.load(open(os.path.join(HERE,'../data/stations_geom.json')))
def ecef(lat,lon,h):
    a=6378.137;f=1/298.257223563;e2=f*(2-f);la,lo=np.radians(lat),np.radians(lon);N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
TX=ecef(GEO['TX']['lat'],GEO['TX']['lon'],GEO['TX']['alt']/1e3)
RX={k:ecef(v['lat'],v['lon'],v['alt']/1e3) for k,v in GEO.items() if k!='TX'}

# Elevation dependence of the gains (antenna maps from aircraft, results/antenna_map_tobit.json): receive gain per
# receiver and elevation bin, transmit gain per elevation bin, relative to the near-zenith geometry of the SPACEMOBILE-008
# calibration (receiver and beacon elevations above 45 deg). Receivers not calibrated on SPACEMOBILE-008 (Humain, Uccle)
# take Oudsbergen's value plus the difference of their aircraft-map offsets.
AM=json.load(open(os.path.join(HERE,'../results/antenna_map_tobit.json'))); EB=AM['el_edges']; EBT=[0,5,10,15,20,30,45,90]
G0=dict(TH['Geff_dB'])
for st in ('BEHUMA','BEUCCE'): G0[st]=TH['Geff_dB']['BEOUDS']+AM['off'][st]-AM['off']['BEOUDS']
def _b(edges,x,n): return int(min(max(np.searchsorted(edges,x,side='right')-1,0),n-1))
def gain_db(rx,el_rx=None,el_tx=None):
    if el_rx is None: return GDB[rx]                 # zenith calibration (model of the paper, Eq. cn0)
    xr=AM['el'][rx]['x']; return G0[rx]+xr[_b(EB,el_rx,len(xr))]-xr[-1]+AM['eltx'][_b(EBT,el_tx,len(AM['eltx']))]-AM['eltx'][-1]
def cn0_per_m2(R_km,rx='BEOUDS',el_rx=None,el_tx=None):
    return TH['Pt_W']*10**(gain_db(rx,el_rx,el_tx)/10)*LAM**2/((4*np.pi)**3*(R_km*1e3)**4*KB*TH['Tsky_K'])
def sigma_min(R_km,T_eff,rx='BEOUDS',el_rx=None,el_tx=None,T_arc=None,z_th=Z_TH):
    """sigma needed for z >= z_th: E/N0 >= z_th sqrt(n_b), n_b = 2 T_arc blocks of 0.5 s (T_arc defaults to T_eff)."""
    Ta=min(T_arc if T_arc is not None else T_eff,50.0)
    return z_th*np.sqrt(2*Ta)/min(T_eff,50.0)/cn0_per_m2(R_km,rx,el_rx,el_tx)
def elev(site,X):
    lat,lon=GEO[site]['lat'],GEO[site]['lon']; la,lo=np.radians(lat),np.radians(lon)
    U=np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)]); d=X-(TX if site=='TX' else RX[site])
    return float(np.degrees(np.arcsin(d@U/np.linalg.norm(d))))
def plate_sigma(A): return 4*np.pi*A**2/LAM**2
RE=6371.0
def glint_s(R_km,side_m,v_kms=7.35):
    """Full width of the glint of a nadir-pointing facet (first null to first null): 2 theta_0 / Omega."""
    Om=v_kms*RE/(R_km*(RE+R_km)); return (LAM/side_m)/Om
from scipy.special import sici
def glint_energy_s(T_arc,T_g):
    """Energy-equivalent duration of a sinc^2 glint lobe (nulls at +-T_g/2) centred in an arc of length T_arc (s)."""
    if not np.isfinite(T_g): return min(T_arc,50.0)
    a=T_g/2; X=min(T_arc,50.0)/(2*a)
    F=sici(2*np.pi*X)[0]/np.pi-np.sin(np.pi*X)**2/(np.pi**2*X)
    return float(2*a*F)

def _sun(jd):
    n=jd-2451545.0; L=np.radians((280.46+0.9856474*n)%360); g=np.radians((357.528+0.9856003*n)%360)
    lam=L+np.radians(1.915)*np.sin(g)+np.radians(0.02)*np.sin(2*g); eps=np.radians(23.439-4e-7*n)
    return np.array([np.cos(lam),np.cos(eps)*np.sin(lam),np.sin(eps)*np.sin(lam)])

def passes(tle,t0,t1,sigma_max=None,plate_area=None,side=None,rx=('BEOUDS','BEHUMA','BEOTTI','FRMAUB','BEUCCE'),step=5.0,elev_gain=False,spin_period=None):
    """Predicted passes of one object between datetimes t0 and t1 (UTC). tle = (line1, line2). One entry per receiver and
    pass: R (km), usable arc T (s, |f_D| in 30-900 Hz), glint (s), T_eff, sigma_min, sigma_max, margin M, sunlit."""
    s=Satrec.twoline2rv(*tle)
    if sigma_max is None: sigma_max=plate_sigma(plate_area)
    if plate_area is not None and side is None: side=np.sqrt(plate_area)
    ts=np.arange(t0.timestamp(),t1.timestamp(),step); jd=ts/86400+2440587.5; ji=np.floor(jd-0.5)+0.5; fr=jd-ji
    out=[]
    for c0 in range(0,len(ts),20000):
        sl=slice(c0,c0+20000); e,r,v=s.sgp4_array(ji[sl],fr[sl]); r=np.asarray(r)
        g=np.array([gstime(a+b) for a,b in zip(ji[sl],fr[sl])]); cg,sg=np.cos(g),np.sin(g)
        X=np.stack([r[:,0]*cg+r[:,1]*sg,-r[:,0]*sg+r[:,1]*cg,r[:,2]],-1); X[np.asarray(e)!=0]=np.nan
        R1=np.linalg.norm(X-TX,axis=1)
        for st in rx:
            R2=np.linalg.norm(X-RX[st],axis=1); rr=np.sqrt(R1*R2); fd=-np.gradient(R1+R2,step)*1e3/LAM
            ok=(np.abs(fd)>30)&(np.abs(fd)<900)&(rr<2000)
            idx=np.nonzero(ok)[0]
            if not len(idx): continue
            for seg in np.split(idx,np.nonzero(np.diff(idx)>int(120/step))[0]+1):   # one pass, even if |f_D| < 30 Hz near closest approach
                T=len(seg)*step
                if T<20: continue
                k=seg[np.nanargmin(rr[seg])]; R=float(rr[k])
                tg=(glint_s(R,side) if side else np.inf) if not spin_period else LAM*spin_period/(2*np.pi*side); Te=glint_energy_s(T,tg)
                er,et=elev(st,X[k]),elev('TX',X[k]); sm=sigma_min(R,Te,st,er,et,T_arc=T) if elev_gain else sigma_min(R,Te,st,T_arc=T); jdk=jd[sl][k]; su=_sun(jdk); p=r[k]; x=p@su
                out.append(dict(rx=st,t=dtm.datetime.fromtimestamp(ts[sl][k],dtm.timezone.utc).isoformat()[:19],R_km=round(R,1),
                    T_s=T,el_rx=round(er,1),el_tx=round(et,1),glint_s=None if not side else round(tg,1),T_eff=round(Te,1),sigma_min=round(sm,1),
                    sigma_max=round(sigma_max,1),M=round(sigma_max/sm,3),sunlit=bool(x>0 or np.linalg.norm(p-x*su)>6371.0)))
    return out

_C=None
def pd_of_M(M,group):
    """Detection probability for margin M from the curve fitted on Phase 1 (code/feasibility_validate.py), groups
    'structure' (ISS-like, M from its typical effective cross section), 'plate_shadow', 'plate_sunlit'."""
    global _C
    if _C is None: _C=json.load(open(os.path.join(HERE,'../results/feasibility_curve.json')))
    return curve_pd(_C,group,M)
def curve_pd(C,group,M):
    """Logistic in log10 M fitted for M >= 1; false-alarm rate (1e-3) below M = 1."""
    if M<1: return 1e-3
    a,b=C['groups'][group]['a'],C['groups'][group]['b']; return float(max(1/(1+np.exp(-(a+b*np.log10(M)))),1e-3))
def expected(ps,group,sunlit_group=None):
    """Expected detections: sum of P_d(M) over the passes; for nadir plates whose attitude depends on sunlight, pass the
    two groups ('plate_shadow', 'plate_sunlit') as group and sunlit_group."""
    return sum(pd_of_M(p['M'],(sunlit_group if (sunlit_group and p['sunlit']) else group)) for p in ps)

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('tle'); a.add_argument('name'); a.add_argument('--start',required=True)
    a.add_argument('--days',type=float,default=7); a.add_argument('--plate-area',type=float); a.add_argument('--side',type=float)
    a.add_argument('--sigma-max',type=float); a.add_argument('--rx',default='BEOUDS,BEHUMA,BEOTTI,FRMAUB,BEUCCE')
    a.add_argument('--group',default=None,help='structure | plate (default: plate if --plate-area else structure)')
    a.add_argument('--spin-period',type=float,help='tumbling object: glint full width lambda P / (2 pi L)')
    a.add_argument('--elev-gain',action='store_true',help='elevation-dependent gains from the aircraft maps (not better in the validation)')
    q=a.parse_args()
    L=[l.rstrip() for l in open(q.tle) if l.strip()]; tle=None
    for k in range(len(L)-2):
        if (L[k][2:] if L[k].startswith('0 ') else L[k]).strip()==q.name and L[k+1].startswith('1 '): tle=(L[k+1],L[k+2]); break
    if tle is None: raise SystemExit('object not found in '+q.tle)
    t0=dtm.datetime.fromisoformat(q.start).replace(tzinfo=dtm.timezone.utc); t1=t0+dtm.timedelta(days=q.days)
    ps=passes(tle,t0,t1,sigma_max=q.sigma_max,plate_area=q.plate_area,side=q.side,rx=q.rx.split(','),elev_gain=q.elev_gain,spin_period=q.spin_period)
    grp=q.group or ('plate' if q.plate_area else 'structure')
    for p in ps: p['Pd']=pd_of_M(p['M'],grp if grp=='structure' else ('plate_sunlit' if p['sunlit'] else 'plate_shadow'))
    ok=[p for p in ps if p['M']>=1]
    print(f'{q.name}: {len(ps)} passes in band, {len(ok)} with margin M >= 1')
    for p in sorted(ok,key=lambda p:(-p['Pd'],-p['M']))[:15]:
        print('  %s %s R=%6.0f km T_eff=%4.0f s sigma_min=%7.0f m2 M=%7.2f %s Pd=%.2f'%(p['t'],p['rx'],p['R_km'],p['T_eff'],p['sigma_min'],p['M'],'sunlit' if p['sunlit'] else 'shadow',p['Pd']))
    print('expected detections in %.0f days: %.1f (false alarms over the same passes: %.1f)'%(q.days,sum(max(p['Pd']-1e-3,0) for p in ps),1e-3*len(ps)))
