# Exploratory (not in PROTOCOL.md): detection probability with the satellite in sunlight or in the Earth's shadow at its
# closest approach to the receiver-beacon pair. Applied to Phase 2 (where the effect was noticed) and to Phase 1 as an
# independent check. Shadow: cylindrical Earth shadow (radius 6371 km), low-precision solar position (Astronomical
# Almanac); closest approach from the element set on a 5-s grid over the file, with the file start from its name.
import json, glob, re, datetime, collections, numpy as np
from scipy.stats import beta
from sgp4.api import Satrec
def cp(k,n,a=0.05): return (0.0 if k==0 else float(beta.ppf(a/2,k,n-k+1)), 1.0 if k==n else float(beta.ppf(1-a/2,k+1,n-k)))
G=json.load(open('../data/stations_geom.json'))
def ecef(lat,lon,h):
    a=6378.137;f=1/298.257223563;e2=f*(2-f);la,lo=np.radians(lat),np.radians(lon);N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
TX=ecef(G['TX']['lat'],G['TX']['lon'],G['TX']['alt']/1e3)
RX={k:ecef(v['lat'],v['lon'],v['alt']/1e3) for k,v in G.items() if k!='TX'}
def sunvec(jd):
    n=jd-2451545.0; L=np.radians((280.46+0.9856474*n)%360); g=np.radians((357.528+0.9856003*n)%360)
    lam=L+np.radians(1.915)*np.sin(g)+np.radians(0.02)*np.sin(2*g); eps=np.radians(23.439-4e-7*n)
    return np.array([np.cos(lam),np.cos(eps)*np.sin(lam),np.sin(eps)*np.sin(lam)])
def jd_of(f): return (datetime.datetime.strptime(f[:13],'%Y%m%d_%H%M')-datetime.datetime(1970,1,1)).total_seconds()/86400+2440587.5
def sunlit(s,f):
    jd0=jd_of(f); best=None; rx=RX[f[-6:]]
    for dt in np.arange(0,300,5):
        jd=jd0+dt/86400; e,p,v=s.sgp4(int(jd),jd-int(jd))
        if e: continue
        p=np.array(p); gm=np.radians((280.46061837+360.98564736629*(jd-2451545))%360)
        pe=np.array([[np.cos(gm),np.sin(gm),0],[-np.sin(gm),np.cos(gm),0],[0,0,1]])@p
        rr=np.linalg.norm(pe-TX)*np.linalg.norm(pe-rx)
        if best is None or rr<best[0]: best=(rr,p,jd)
    _,p,jd=best; su=sunvec(jd); x=p@su
    return bool(x>0 or np.linalg.norm(p-x*su)>6371.0)
def sets_from(files):
    S={}
    for f in files:
        L=[l.rstrip() for l in open(f) if l.strip()]
        for k in range(len(L)-2):
            if L[k+1].startswith('1 ') and L[k+2].startswith('2 ') and not L[k].startswith(('1 ','2 ')):
                nm=L[k].strip(); nm=nm[2:] if nm.startswith('0 ') else nm; S.setdefault(nm.strip(),[]).append((L[k+1],L[k+2]))
    return S
def picker(S):
    EP={}
    def sat(nm,f):
        if nm not in EP: EP[nm]=np.array([Satrec.twoline2rv(*q).jdsatepoch+Satrec.twoline2rv(*q).jdsatepochF for q in S[nm]])
        return Satrec.twoline2rv(*S[nm][int(np.argmin(abs(EP[nm]-jd_of(f))))])
    return sat
out={}
# Phase 2
g={}; exec(open('phase2_analyze.py').read().split("# echo track of an owner")[0],g); P2=g['P']
D2={(d['file'],d['obj']) for d in json.load(open('../results/phase2_summary.json'))['detections']}
sat2=picker(sets_from(sorted(glob.glob('../data/tle/phase2/*.txt'))))
# Phase 1 (catalog snapshot plus the objects with their own element sets, as in blind_enumerate.py)
g1={}; exec(open('blind_analyze.py').read().split("# echo track of an owner")[0],g1); P1=g1['P']
D1={(d['file'],d['obj']) for d in json.load(open('../results/blind_summary.json'))['detections']}
S1=sets_from(['../data/tle/leo.txt']+sorted(glob.glob('../data/tle/*.txt')))
sat1=picker(S1)
CL=('ISS complex','AST SpaceMobile','Starlink DTC','Starlink V2-mini','Starlink V1.x')
for ph,PP,DD,sat in (('phase2',P2,D2,sat2),('phase1',P1,D1,sat1)):
    res=collections.defaultdict(lambda:[0,0])
    for r in PP:
        if r['class'] not in CL: continue
        try: lit=sunlit(sat(r['obj'],r['file']),r['file'])
        except KeyError: continue
        b=next(i for i,(a,c) in enumerate(zip([0,30,100,300,1e3,3e3,1e4],[30,100,300,1e3,3e3,1e4,1e12])) if a<=r['sigma_min']<c)
        h=int(r['file'][9:11])
        for key in ((r['class'],lit,'all'),(r['class'],lit,b),(r['class'],lit,'h%02d'%(h//3*3))):
            res[key][0]+=(r['file'],r['obj']) in DD; res[key][1]+=1
    out[ph]=[dict(cls=c,sunlit=l,bin=b,k=v[0],n=v[1],ci=cp(*v)) for (c,l,b),v in sorted(res.items(),key=lambda x:str(x[0]))]
    for d in out[ph]:
        if d['bin']=='all': print(ph,'%-17s %-7s %4d/%5d = %.3f [%.3f,%.3f]'%(d['cls'],'sunlit' if d['sunlit'] else 'shadow',d['k'],d['n'],d['k']/d['n'],*d['ci']))
json.dump(out,open('../results/eclipse.json','w'),indent=0)
