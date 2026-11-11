# Which attitude explains the sunlight effect? (exploratory) For each Phase 2 pass of a Starlink
# V2-mini or DTC satellite, at closest approach to the transmitter-receiver pair: the unit bisector n of the directions
# from the satellite to the transmitter and to the receiver (the normal a flat reflector must have to return a specular
# echo), and the angle between n and four candidate normals of a large flat surface:
#   nadir      the Earth-facing side of the bus or of arrays held flat (normal toward the Earth's centre)
#   sun        a plate facing the Sun (either face: the angle to the line through the Sun)
#   sun_ax_v   a plate that tracks the Sun by rotating about the along-track axis (normal = Sun direction projected on
#              the plane perpendicular to the velocity)
#   sun_ax_h   the same, rotating about the cross-track (orbit-normal) axis
# Detection probability of sunlit and shadowed passes is then given in bins of each angle; a candidate that is the
# attitude of the reflector should concentrate the detections at small angles.
import json, glob, numpy as np
exec(open('stats_cluster.py').read().split("out=dict(n_boot")[0])       # U2m (Phase 2 own thresholds), sat2, geometry tools
def ang(a,b): return float(np.degrees(np.arccos(np.clip(abs(a@b)/np.linalg.norm(a)/np.linalg.norm(b),-1,1))))
def geom(s,f):
    jd0=jd_of(f); best=None; rx=RX[f[-6:]]
    for dt in np.arange(0,300,2):
        jd=jd0+dt/86400; e,p,v=s.sgp4(int(jd),jd-int(jd))
        if e: continue
        p=np.array(p); v=np.array(v); gm=np.radians((280.46061837+360.98564736629*(jd-2451545))%360)
        Rz=np.array([[np.cos(gm),np.sin(gm),0],[-np.sin(gm),np.cos(gm),0],[0,0,1]]); pe=Rz@p
        rr=np.linalg.norm(pe-TX)*np.linalg.norm(pe-rx)
        if best is None or rr<best[0]: best=(rr,p,v,pe,Rz,jd)
    _,p,v,pe,Rz,jd=best; su=sunvec(jd)                       # inertial (TEME ~ equator of date) vectors
    tx=Rz.T@TX; rxi=Rz.T@rx                                   # stations rotated to the same frame
    n=(tx-p)/np.linalg.norm(tx-p)+(rxi-p)/np.linalg.norm(rxi-p)
    vh=v/np.linalg.norm(v); h=np.cross(p,v); h/=np.linalg.norm(h)
    pv=su-(su@vh)*vh; ph=su-(su@h)*h
    return dict(nadir=ang(n,-p),sun=ang(n,su),sun_ax_v=ang(n,pv),sun_ax_h=ang(n,ph))
rows=[]
for u in U2m:
    if u['cls'] not in ('Starlink V2-mini','Starlink DTC','Starlink V1.x'): continue
    try: g=geom(sat2(u['r']['obj'],u['file']),u['file'])
    except KeyError: continue
    rows.append(dict(cls=u['cls'],grp=u['grp'],det=bool(u['det']),**g))
out={'bins_deg':[0,5,10,20,40,90]}
for c in ('Starlink V2-mini','Starlink DTC','Starlink V1.x'):
    for gp in ('plate_shadow','plate_sunlit'):
        q=[r for r in rows if r['cls']==c and r['grp']==gp]
        for key in ('nadir','sun','sun_ax_v','sun_ax_h'):
            res=[]
            for a,b in zip(out['bins_deg'][:-1],out['bins_deg'][1:]):
                s=[r for r in q if a<=r[key]<b]; k=sum(r['det'] for r in s); res.append(dict(lo=a,hi=b,n=len(s),k=int(k),pd=k/len(s) if s else None,ci=cp(k,len(s)) if s else None))
            out[f'{c}|{gp}|{key}']=res
json.dump(out,open('../results/sun_geometry.json','w'),indent=1)
for k,v in out.items():
    if k=='bins_deg': continue
    print('%-42s'%k,' '.join('%g-%g:%d/%d'%(r['lo'],r['hi'],r['k'],r['n']) for r in v))
