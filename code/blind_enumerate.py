# Blind search, step 1: every pass of every catalogued object in every recording. For each (object, recording) pair we keep
# the geometry summary needed to decide whether a search is meaningful: seconds of the file during which the predicted
# Doppler lies in the usable band (30 Hz < |f_D| < 900 Hz), and the smallest geometric mean range during those seconds.
import os, glob, json, numpy as np, datetime as dtm
from sgp4.api import Satrec, SatrecArray, jday
from sgp4.propagation import gstime
import bwav
c=299792458.; lam=c/49.97e6
def ecef(lat,lon,h):
    a=6378.137;f=1/298.257223563;e2=f*(2-f);la,lo=np.radians(lat),np.radians(lon);N=a/np.sqrt(1-e2*np.sin(la)**2)
    return np.array([(N+h)*np.cos(la)*np.cos(lo),(N+h)*np.cos(la)*np.sin(lo),(N*(1-e2)+h)*np.sin(la)])
def read_sets(path):
    L=[l.rstrip() for l in open(path) if l.strip()]; out=[]
    for k in range(len(L)-2):
        if L[k+1].startswith('1 ') and L[k+2].startswith('2 '):
            nm=L[k].strip(); nm=nm[2:] if nm.startswith('0 ') else nm; out.append((nm,L[k+1],L[k+2]))
    return out
cat=read_sets('../data/tle/leo.txt')
special={}                                     # objects with several element sets: pick the epoch closest to each file
for f in ['iss.txt','bluebird_gp_history.txt','envisat_gp_history.txt']:
    for nm,l1,l2 in read_sets('../data/tle/'+f):
        if f=='iss.txt': nm='ISS'
        special.setdefault(nm,[]).append((l1,l2))
cat=[x for x in cat if x[0] not in special and not x[0].startswith('ISS (ZARYA)')]
names=[x[0] for x in cat]; arr=SatrecArray([Satrec.twoline2rv(x[1],x[2]) for x in cat])
files=sorted(glob.glob(os.environ.get('WAVDIR','../data/raw_wav/')+'*.wav')+glob.glob(os.environ.get('EVDIR','../data/raw_wav_envisat/')+'*.wav'))
out=[]
for fp in files:
    md,_,_=bwav.read(fp); t0=md['start_us']/1e6
    TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
    ts=t0+np.arange(0,300.01,1.0); jd=[];fr=[]
    for t in ts:
        d=dtm.datetime.fromtimestamp(t,dtm.timezone.utc); a,b=jday(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond/1e6); jd.append(a); fr.append(b)
    jd=np.array(jd); fr=np.array(fr); g=np.array([gstime(a+b) for a,b in zip(jd,fr)])
    def summar(rr_all,fd_all,nm_list,extra=None):
        valid=(np.abs(fd_all)>30)&(np.abs(fd_all)<900)
        nv=valid.sum(1); rrm=np.where(valid,rr_all,np.inf).min(1)
        for i in np.flatnonzero((nv>=5)&(rrm<1500)):
            out.append(dict(file=os.path.basename(fp)[11:31],obj=nm_list[i],n_valid_s=int(nv[i]),rrmin_valid=float(rrm[i]),rrmin=float(rr_all[i].min())))
    def geom(r):
        X=np.stack([r[...,0]*np.cos(g)+r[...,1]*np.sin(g),-r[...,0]*np.sin(g)+r[...,1]*np.cos(g),r[...,2]],-1)
        R1=np.linalg.norm(X-TX,axis=-1); R2=np.linalg.norm(X-RX,axis=-1); fd=-np.gradient(R1+R2,1.0,axis=-1)*1e3/lam      # ranges in km, lam in m
        return np.sqrt(R1*R2),fd
    e,r,v=arr.sgp4(jd,fr); ok=(e==0).all(1); rr,fd=geom(r); rr[~ok]=np.inf
    summar(rr,fd,names)
    for nm,sets in special.items():
        ep=np.array([Satrec.twoline2rv(*s).jdsatepoch+Satrec.twoline2rv(*s).jdsatepochF for s in sets]); s=sets[int(np.argmin(abs(ep-(jd[0]+fr[0]))))]
        sat=Satrec.twoline2rv(*s); rs=np.array([sat.sgp4(a,b)[1] for a,b in zip(jd,fr)])[None]; rr1,fd1=geom(rs); summar(rr1,fd1,[nm])
    print(os.path.basename(fp)[11:31],len(out),flush=True)
json.dump(out,open('../results/blind_passes.json','w'))
