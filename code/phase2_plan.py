# Phase 2 (validation campaign, 19 Oct - 2 Nov 2026): which 5-min recordings to retrieve and which objects to search in them.
# Target classes: ISS complex, AST SpaceMobile, Starlink (DTC, V2-mini, and V1.x as a control). For every object, day and
# receiver, the predicted geometry on a 10-s grid gives the 5-min slots in which the object stays >= 20 s in the usable
# Doppler band, and the cross section sigma_min it would need for z = 7 (radar equation, PROTOCOL item 5).
# A pass is physically detectable when sigma_min <= sigma_max, the specular peak 4 pi A^2 / lambda^2 of the total flat
# area A of its type (in Phase 1, such passes held 358 of the 370 detections, and the others were detected at 0.14%,
# the false-alarm rate). Recordings retrieved: every slot with a physically detectable target pass at Oudsbergen or Humain,
# and at all five receivers every slot with an ISS-complex pass. In each retrieved recording every target pass within
# 1500 km is searched, so that the impossible passes serve as controls.
import json, re, glob, numpy as np, datetime as dtm
from sgp4.api import Satrec, SatrecArray
from sgp4.propagation import gstime
LAM=299792458./49.97e6
G=json.load(open('../data/stations_geom.json'))
T=json.load(open('../results/detect_theory.json')); GE={k:10**(v/10) for k,v in T['Geff_dB'].items()}
GE['BEHUMA']=GE['BEOUDS']; GE['BEUCCE']=GE['BEOTTI']      # not calibrated on SPACEMOBILE-008: assumed for planning only
def ecef(lat,lon,h):
    a=6378.137;f=1/298.257223563;e2=f*(2-f);la,lo=np.radians(lat),np.radians(lon);N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
TX=ecef(G['TX']['lat'],G['TX']['lon'],G['TX']['alt']/1000); RX={k:ecef(v['lat'],v['lon'],v['alt']/1000) for k,v in G.items() if k!='TX'}
ISSC=('ISS (','POISK','PROGRESS-MS','SOYUZ-MS','CREW DRAGON','DRAGON CRS','CYGNUS','HTV-X')
DTC={l.split()[0] for l in open('../data/tle/leo.txt') if '[DTC]' in l}   # Space-Track names lack the CelesTrak [DTC] tag
def cls(n):
    if n in DTC: return 'Starlink DTC'
    if n.startswith('ISS') and not n.startswith('ISS OBJECT') or any(n.startswith(p) for p in ISSC): return 'ISS complex'
    if n.startswith('SPACEMOBILE') or n.startswith('BLUEWALKER'): return 'AST SpaceMobile'
    if '[DTC]' in n: return 'Starlink DTC'
    m=re.match(r'STARLINK-(\d+)',n)
    if m: return 'Starlink V2-mini' if int(m.group(1))>=30000 else 'Starlink V1.x'
    return None
sets={}
for f in sorted(glob.glob('../data/tle/phase2/*.txt')):
    L=[l.rstrip('\n') for l in open(f)]
    for k in range(len(L)-2):
        if L[k].startswith('0 ') and L[k+1].startswith('1 ') and L[k+2].startswith('2 '):
            nm=L[k][2:].strip(); c=cls(nm)
            if c: sets.setdefault(nm,{})[L[k+1][18:32]]=(L[k+1],L[k+2])
names=sorted(sets); print('objects',len(names),{c:sum(cls(n)==c for n in names) for c in {cls(n) for n in names}},flush=True)
JD0={y:dtm.datetime(2000+y,1,1,tzinfo=dtm.timezone.utc).timestamp()/86400+2440587.5 for y in range(20,40)}
days=[dtm.datetime(2026,10,19,tzinfo=dtm.timezone.utc)+dtm.timedelta(days=i) for i in range(15)]
plan=[]
for d0 in days:
    jdm=d0.timestamp()/86400+2440587.5+0.5
    sats=[];nm_ok=[]
    for nm in names:
        ep=list(sets[nm]); jde=np.array([JD0[int(x[:2])]+float(x[2:])-1 for x in ep])
        sats.append(Satrec.twoline2rv(*sets[nm][ep[int(np.argmin(abs(jde-jdm)))]])); nm_ok.append(nm)
    t=d0.timestamp()+np.arange(0,86400,10.0); jd=t/86400+2440587.5; jdi=np.floor(jd-0.5)+0.5; fr=jd-jdi
    g=np.array([gstime(a+b) for a,b in zip(jdi,fr)]); cg,sg=np.cos(g),np.sin(g)
    for i0 in range(0,len(sats),1500):
        arr=SatrecArray(sats[i0:i0+1500]); e,r,v=arr.sgp4(jdi,fr)
        X=np.stack([r[...,0]*cg+r[...,1]*sg,-r[...,0]*sg+r[...,1]*cg,r[...,2]],-1); X[e!=0]=np.nan
        R1=np.linalg.norm(X-TX,axis=-1)
        for st,P in RX.items():
            R2=np.linalg.norm(X-P,axis=-1); rr=np.sqrt(R1*R2); fd=-np.gradient(R1+R2,10.0,axis=-1)*1e3/LAM
            ok=(np.abs(fd)>30)&(np.abs(fd)<900)&(rr<=1500)
            if not ok.any(): continue
            for i,s in zip(*np.nonzero(ok.reshape(ok.shape[0],288,30).sum(-1)>=2)):   # >= 20 s in band within the slot
                w=ok[i,s*30:(s+1)*30]; rmin=float(np.nanmin(np.where(w,rr[i,s*30:(s+1)*30],np.inf))); nv=int(w.sum())*10
                cn0=T['Pt_W']*GE[st]*T['lam_m']**2/((4*np.pi)**3*(rmin*1e3)**4*1.380649e-23*T['Tsky_K'])
                smin=70/min(nv,50)/cn0; nm=nm_ok[i0+i]
                lamw=299792458/49.97e6; c=cls(nm)
                smx=1e5 if c=='ISS complex' else 4*np.pi*((64 if nm.startswith('BLUEWALKER') else 223) if c=='AST SpaceMobile' else 116 if c in ('Starlink DTC','Starlink V2-mini') else 29)**2/lamw**2
                # selection: physically possible (sigma_min <= sigma_max) AND detection probability measured in Phase 1 of at
                # least ~1% (ten times the false-alarm rate): sigma_min <= 1000 m^2 for DTC and V2-mini, <= 3000 m^2 for AST,
                # <= 300 m^2 for V1.x (control class), any for the ISS complex within 1500 km
                lim={'ISS complex':1e9,'AST SpaceMobile':3000,'Starlink DTC':1000,'Starlink V2-mini':1000,'Starlink V1.x':300}[c]
                if smin<=min(lim,smx):
                    slot=(d0+dtm.timedelta(minutes=5*int(s))).strftime('%Y%m%d_%H%M')
                    plan.append((f'{slot}_{st}',nm,c,round(rmin,1),nv,round(smin,1),round(smx,1)))
    print(d0.date(),len(plan),flush=True)
json.dump(plan,open('../results/phase2_plan_passes.json','w'))
lam=299792458/49.97e6
def smax(p):
    if p['cls']=='ISS complex': return 1e5
    if p['cls']=='AST SpaceMobile': return 4*np.pi*(64 if p['obj'].startswith('BLUEWALKER') else 223)**2/lam**2
    if p['cls'] in ('Starlink DTC','Starlink V2-mini'): return 4*np.pi*116**2/lam**2
    return 4*np.pi*29**2/lam**2
plan=[dict(file=a,obj=b,cls=c,rrmin=d,n_valid_s=e,sigma_min=f,sigma_max=g) for a,b,c,d,e,f,g in plan]
json.dump(plan,open('../results/phase2_plan_passes.json','w'))
sel=set()
for p in plan:
    st=p['file'][-6:]
    if p['cls']=='ISS complex' or st in ('BEOUDS','BEHUMA'): sel.add(p['file'])
sel=sorted(sel); json.dump(sel,open('../results/phase2_files.json','w'))
import collections
print('files',len(sel),collections.Counter(f[-6:] for f in sel),'~GB %.1f'%(len(sel)*3.66e-3))
pp=[p for p in plan if p['file'] in set(sel)]
print('passes to search',len(pp),collections.Counter(p['cls'] for p in pp))
