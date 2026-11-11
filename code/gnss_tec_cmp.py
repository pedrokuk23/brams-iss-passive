# Comparação do termo ionosférico ajustado com o TEC vertical de GNSS (mapas ROB/EPN, 0,5°, 5 min)
# 1) VTEC GNSS nos pontos de perfuração (300 km) dos dois trechos TX-ISS-RX ao longo do arco rastreado;
# 2) VTEC total do perfil de Chapman ajustado e fração abaixo da ISS;
# 3) ajuste de δt com N_mF2 fixado pelo VTEC GNSS (hmF2 = 300, 350, 420 km).
import numpy as np, json, datetime as dtm
exec(open('iono_sens_v2.py').read().split("out={}; maxchg=0")[0])
# ---- mapas GNSS
G={}
for l in open('../data/gnss_tec/rob_vtec_subset.txt'):
    if l.startswith('#') or not l.strip(): continue
    k,v=l.strip().split(':'); G[k]=np.array([[int(a) for a in r.split(',')] for r in v.strip(';').split(';')])/10.
LAT=np.arange(48.5,52.51,0.5); LON=np.arange(2.5,7.01,0.5)
def tkey(d):  # datetime -> chave e tempo do mapa
    doy=d.timetuple().tm_yday; return f"{doy}{chr(65+d.hour)}{d.minute:02d}"
def vtec(t_unix,la,lo):
    d=dtm.datetime.fromtimestamp(t_unix,dtm.timezone.utc); m0=d.replace(minute=d.minute//5*5,second=0,microsecond=0); m1=m0+dtm.timedelta(minutes=5)
    w=(d-m0).total_seconds()/300; out=[]
    for m in (m0,m1):
        k=tkey(m); A=G[k]
        i=np.interp(la,LAT,np.arange(len(LAT))); j=np.interp(lo,LON,np.arange(len(LON)))
        i0,j0=int(min(i,len(LAT)-2)),int(min(j,len(LON)-2)); fi,fj=i-i0,j-j0
        out.append(A[i0,j0]*(1-fi)*(1-fj)+A[i0+1,j0]*fi*(1-fj)+A[i0,j0+1]*(1-fi)*fj+A[i0+1,j0+1]*fi*fj)
    return out[0]*(1-w)+out[1]*w
def llh(X):
    r=np.linalg.norm(X); return np.degrees(np.arcsin(X[2]/r)),np.degrees(np.arctan2(X[1],X[0])),r-RE
def ipp(A,B,h=300.):   # ponto do segmento A->B na altitude h
    s=np.linspace(0,1,2001)[:,None]; P=A[None,:]+s*(B-A)[None,:]; k=np.argmin(abs(np.linalg.norm(P,axis=1)-RE-h)); return P[k]
hh=np.linspace(60,1500,3000)
res={}
for p,keys in passes.items():
    if not keys: continue
    v=[];hs=[];grad=[]
    for key in keys:
        a,tc,t,X,TX,RX,dfd=geo[key]; md,_,_=bwav.read(D+f'RAD_BEDOUR_{key}_SYS001.wav'); t0=md['start_us']/1e6
        for tt in np.linspace(tc.min(),tc.max(),7):
            Xs=X[np.argmin(abs(t-tt))]; hs.append(llh(Xs)[2])
            for P in (ipp(TX,Xs),ipp(Xs,RX)):
                la,lo,_=llh(P); v.append(vtec(t0+tt,la,lo))
                gN=(vtec(t0+tt,la+0.25,lo)-vtec(t0+tt,la-0.25,lo))/0.5/111.2; gE=(vtec(t0+tt,la,lo+0.25)-vtec(t0+tt,la,lo-0.25))/0.5/(111.2*np.cos(np.radians(la)))
                grad.append((gN,gE))
    V=float(np.mean(v)); hiss=float(np.mean(hs)); g=np.mean(grad,axis=0)
    r=dict(receivers=[k[14:] for k in keys],vtec_gnss=V,vtec_range=[float(min(v)),float(max(v))],h_iss_km=hiss,grad_tecu_per_100km=[float(g[0]*100),float(g[1]*100)],fits={})
    base=fit(p,300)[3]
    for hm in [300,350,420]:
        tot=np.trapezoid(chap(hh,hm,55),hh)*1e3; below=np.trapezoid(chap(hh[hh<hiss],hm,55),hh[hh<hiss])*1e3   # m (por unidade de Nm)
        d,Nm,rms,dev=fit(p,hm)
        Nfix=V*1e16/tot; d2,_,rms2,dev2=fit(p,hm,Nfix=Nfix)
        r['fits'][hm]=dict(free=dict(Nm=float(Nm),vtec_total=float(Nm*tot/1e16),vtec_below_iss=float(Nm*below/1e16),rms=float(rms),dt_minus_oem=[float(x) for x in dev]),
                           frac_below_iss=float(below/tot),
                           gnss_fixed=dict(Nm=float(Nfix),rms=float(rms2),dt_minus_oem=[float(x) for x in dev2],change_ms=float(np.max(np.abs(np.array(dev2)-base)))))
    res[p]=r
    f=r['fits']
    print(f"{p} GNSS VTEC {V:.1f} TECU ({min(v):.1f}-{max(v):.1f}), grad N {g[0]*100:+.2f} E {g[1]*100:+.2f} TECU/100 km, h_ISS {hiss:.0f} km")
    for hm in [300,350,420]:
        q=f[hm]; print(f"   hm={hm}: free VTEC tot {q['free']['vtec_total']:.1f} below ISS {q['free']['vtec_below_iss']:.1f} (frac {q['frac_below_iss']:.2f}) rms {q['free']['rms']:.2f} | GNSS-fixed rms {q['gnss_fixed']['rms']:.2f} dt change {q['gnss_fixed']['change_ms']:.1f} ms")
json.dump(res,open('../results/gnss_tec_cmp.json','w'),indent=1)
