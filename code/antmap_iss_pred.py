# Nível esperado do eco da ISS em cada estação a partir do mapa de antena obtido com aviões
import numpy as np, json, sys, datetime as dt
sys.argv=['x','X']
exec(open('station_run.py').read().split('print(f"== {ST} ==")')[0])
FIT=json.load(open('../results/antenna_map_tobit.json'))
C=json.load(open('../data/stations/coords.json')); TXg=C['tx']; TX=ecef(TXg[0],TXg[1],0.236)
def enu(lat,lon,X0,X):
    la,lo=np.radians(lat),np.radians(lon); d=X-X0
    E=np.array([-np.sin(lo),np.cos(lo),0]); N=np.array([-np.sin(la)*np.cos(lo),-np.sin(la)*np.sin(lo),np.cos(la)]); U=np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
    e,n,u=d@E,d@N,d@U; return np.degrees(np.arctan2(e,n))%360, np.degrees(np.arctan2(u,np.hypot(e,n)))
EB=np.array(FIT['el_edges'])
def grx(st,az,el):
    a=FIT['az'][st]['high' if el>=10 else 'low'][int(az//30)%12]
    return FIT['off'][st]+a+FIT['el'][st]['x'][min(np.searchsorted(EB,el,side='right')-1,7)]
def gtx(elt): return FIT['eltx'][min(np.searchsorted([0,5,10,15,20,30,45,90],elt,side='right')-1,6)]
def gb(b): return FIT['beta'][min(int(b//10),17)]
TBL=json.load(open('../results/iss_detection_table.json')); Z={(r['pass_'],r['st']):r['z'] for r in TBL}
out=[]
for p in sorted({r['pass_'] for r in TBL}):
    t0=dt.datetime.strptime(p,'%Y%m%d_%H%M').replace(tzinfo=dt.timezone.utc).timestamp()
    t=np.arange(t0-300,t0+600,1.0); X=iss_ecef(t)
    for st,(la0,lo0) in C['stations'].items():
        if (p,st) not in Z: continue
        RX=ecef(la0,lo0,0.1); R1=np.linalg.norm(X-TX,axis=1); R2=np.linalg.norm(X-RX,axis=1)
        fd=-np.gradient(R1+R2,t)/lam; m=(np.abs(fd)<1060)&(np.abs(fd)>3)
        rows=[]
        for i in np.where(m)[0]:
            az,el=enu(la0,lo0,RX,X[i]); _,elt=enu(TXg[0],TXg[1],TX,X[i])
            if el<0 or elt<0: continue
            u1=(X[i]-TX)/R1[i]; u2=(X[i]-RX)/R2[i]; b=np.degrees(np.arccos(np.clip(u1@u2,-1,1)))
            pl=-20*np.log10(R1[i]*R2[i])
            rows.append((pl,grx(st,az,el),gtx(elt),gb(b),az,el,b))
        if not rows: continue
        A=np.array(rows)
        lin=lambda v:10*np.log10(np.sum(10**(v/10)))
        out.append(dict(pass_=p,st=st,z=Z[(p,st)],n=len(A),
            pl_only=lin(A[:,0]),with_rx=lin(A[:,0]+A[:,1]),with_rx_tx=lin(A[:,0]+A[:,1]+A[:,2]),
            az=[float(A[:,4].min()),float(A[:,4].max())],el=[float(A[:,5].min()),float(A[:,5].max())],grx_mean=float(A[:,1].mean())))
json.dump(out,open('../results/antmap/iss_pred.json','w'))
for p in sorted({o['pass_'] for o in out}):
    ref=[o for o in out if o['pass_']==p and o['st']=='BEOTTI'][0]
    for o in [o for o in out if o['pass_']==p]:
        print(f"{p} {o['st']}: z={o['z']:5.1f}  rel.Ottignies: só distância {o['pl_only']-ref['pl_only']:+5.1f} dB | +antena RX {o['with_rx']-ref['with_rx']:+5.1f} | +TX {o['with_rx_tx']-ref['with_rx_tx']:+5.1f}   el {o['el'][0]:.0f}-{o['el'][1]:.0f} az {o['az'][0]:.0f}-{o['az'][1]:.0f} Grx {o['grx_mean']:+.1f}")
