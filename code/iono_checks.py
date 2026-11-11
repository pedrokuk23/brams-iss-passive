# Checks of the ionospheric part of the along-track fit:
# (1) WGS84 geodetic height instead of spherical height (|r| - 6371 km) inside the Chapman layer;
# (2) identifiability of delta-t versus the ionospheric term: condition number of the column-scaled design matrix,
#     correlation of the two estimates, and the relative size of the parts of the two signatures that cannot be told apart;
# (3) injection tests: residuals simulated with ionospheres different from the fitted model (other hmF2, other H,
#     horizontal gradient, topside-heavy profile), no timing error, fitted with the published model -> bias of delta-t.
import numpy as np, json, os
exec(open('iono_fit.py').read().split("tr=json.load")[0])
tr=json.load(open('../results/iono_tracks_v2.json')); D=os.environ.get('WAVDIR','../data/raw_wav/')
PPS={r['file'][11:31]:r['off0_ms'] for r in json.load(open('../results/pps_audit.json'))}
NMIN=19
c=299792458.0; f0=49.97e6; K=40.3/(c*f0); RE=6371.0
def cleaned(a):
    a=np.array(a); med=np.median(a[:,1]); return a[abs(a[:,1]-med)<4]
def geodetic_h(P):   # ECEF km -> WGS84 ellipsoidal height km (Bowring, one iteration is < 1 m error at these heights)
    a=6378.137; f=1/298.257223563; b=a*(1-f); e2=f*(2-f); ep2=(a*a-b*b)/(b*b)
    x,y,z=P[:,0],P[:,1],P[:,2]; p=np.hypot(x,y); th=np.arctan2(z*a,p*b)
    lat=np.arctan2(z+ep2*b*np.sin(th)**3,p-e2*a*np.cos(th)**3)
    for _ in range(2):
        N=a/np.sqrt(1-e2*np.sin(lat)**2); h=p/np.cos(lat)-N; lat=np.arctan2(z,p*(1-e2*N/(N+h)))
    N=a/np.sqrt(1-e2*np.sin(lat)**2); return p/np.cos(lat)-N, np.degrees(lat)
def prof(h,lat,hm=300.,H=55.,grad=0.,kind='chapman'):
    z=(h-hm)/H
    if kind=='chapman': v=np.exp(0.5*(1-z-np.exp(-z)))
    else:   # topside-heavy: scale height doubles above the peak
        H2=np.where(h>hm,2*H,H); z=(h-hm)/H2; v=np.exp(0.5*(1-z-np.exp(-z)))
    return v*(1+grad*(lat-50.6)*111.2/100)    # grad = fractional change per 100 km northward
def stec_g(A,B,geod,**kw):
    s=np.linspace(0,1,400)[:,None]; P=A[None,:]+s*(B-A)[None,:]
    if geod: h,lat=geodetic_h(P)
    else: h=np.linalg.norm(P,axis=1)-RE; lat=np.degrees(np.arcsin(P[:,2]/np.linalg.norm(P,axis=1)))
    return np.trapezoid(prof(h,lat,**kw),dx=1/399)*np.linalg.norm(B-A)*1e3
geo={}
def g(key):
    if key in geo: return geo[key]
    md,x,pps=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    a=cleaned(tr[key]); tc=a[:,0]; t=np.arange(tc.min()-5,tc.max()+5,0.5); X=iss_ecef(t0+t)
    Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
    mt=-np.interp(tc,t,dfd); obs=a[:,1]-mt*PPS[key]*1e-3
    geo[key]=(obs,mt,tc,t,X,TX,RX); return geo[key]
def icol(key,geod,**kw):
    obs,mt,tc,t,X,TX,RX=g(key)
    S=np.array([stec_g(TX,p,geod,**kw)+stec_g(p,RX,geod,**kw) for p in X]); return K*np.interp(tc,t,np.gradient(S,t))
def lsq(A,y):
    p,*_=np.linalg.lstsq(A,y,rcond=None); r=y-A@p; rms=np.sqrt(np.mean(r**2))
    C=np.linalg.inv(A.T@A)*rms**2*len(y)/max(len(y)-2,1); return p,C,rms
passes=sorted({k[:13] for k in tr}); out={}
for p in passes:
    keys=[k for k in tr if k.startswith(p) and len(cleaned(tr[k]))>=NMIN]
    if not keys: continue
    obs=np.concatenate([g(k)[0] for k in keys]); mt=np.concatenate([g(k)[1] for k in keys])
    mis=np.concatenate([icol(k,False) for k in keys]); mig=np.concatenate([icol(k,True) for k in keys])
    ps,Cs,rs=lsq(np.column_stack([mt,mis]),obs); pg,Cg,rg=lsq(np.column_stack([mt,mig]),obs)
    # identifiability
    A=np.column_stack([mt,mis]); As=A/np.linalg.norm(A,axis=0); cond=np.linalg.cond(As)
    corr=Cs[0,1]/np.sqrt(Cs[0,0]*Cs[1,1])
    # fraction of the ionospheric signature that is indistinguishable from a timing error (projection on mt)
    proj=(mis@mt)/(mt@mt); frac=np.linalg.norm(proj*mt)/np.linalg.norm(mis)
    # injections: true delta-t = 0, ionosphere with VTEC ~ 5 TECU (Nm = 2.2e11 m^-3 for H = 55 km), noise = fitted rms
    rng=np.random.default_rng(1); inj={}
    Ninj=2.2e11
    for name,kw in {'hm250':dict(hm=250.),'hm420':dict(hm=420.),'H40':dict(H=40.),'H80':dict(H=80.),
                    'grad-10%':dict(grad=-0.10),'grad+10%':dict(grad=0.10),'topside2H':dict(kind='topside'),
                    'hm420_grad-10%':dict(hm=420.,grad=-0.10)}.items():
        sim=Ninj*np.concatenate([icol(k,True,**kw) for k in keys])
        b0=lsq(np.column_stack([mt,mis]),sim)[0][0]*1e3                     # noise-free bias, ms
        bn=[lsq(np.column_stack([mt,mis]),sim+rng.normal(0,rs,len(sim)))[0][0]*1e3 for _ in range(200)]
        inj[name]=dict(bias_ms=float(b0),sd_noise_ms=float(np.std(bn)))
    out[p]=dict(receivers=[k[14:] for k in keys],n=int(len(obs)),dt_sph_ms=float(ps[0]*1e3),dt_geod_ms=float(pg[0]*1e3),
        d_dt_ms=float((pg[0]-ps[0])*1e3),Nm_sph=float(ps[1]),Nm_geod=float(pg[1]),rms_sph=float(rs),rms_geod=float(rg),
        cond=float(cond),corr_dt_Nm=float(corr),iono_frac_like_timing=float(frac),injection=inj)
    print(p,'dt sph %.1f geod %.1f diff %.2f ms  Nm %.2e->%.2e  cond %.1f corr %.2f frac %.2f'%(ps[0]*1e3,pg[0]*1e3,(pg[0]-ps[0])*1e3,ps[1],pg[1],cond,corr,frac))
    print('   injection bias ms:',{k:round(v['bias_ms'],1) for k,v in inj.items()})
allb=[abs(v['bias_ms']) for o in out.values() for v in o['injection'].values()]
out['summary']=dict(max_abs_d_dt_geodetic_ms=max(abs(o['d_dt_ms']) for o in out.values()),
    cond_range=[min(o['cond'] for o in out.values() if 'cond' in o),max(o['cond'] for o in out.values() if 'cond' in o)],
    corr_range=[min(o['corr_dt_Nm'] for o in out.values() if 'cond' in o),max(o['corr_dt_Nm'] for o in out.values() if 'cond' in o)],
    max_abs_injection_bias_ms=max(allb),median_abs_injection_bias_ms=float(np.median(allb)))
print(out['summary'])
json.dump(out,open('../results/iono_checks.json','w'),indent=1)
