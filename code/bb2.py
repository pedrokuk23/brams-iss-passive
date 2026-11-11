# Busca de ecos dos BlueBird (AST SpaceMobile) com a mesma estatística da ISS (normalização na grade), grade ampliada
import numpy as np, json, sys, glob, os, datetime as dtm
from sgp4.api import Satrec, jday
from sgp4.propagation import gstime
exec(open('synth.py').read().split("if __name__")[0])
L=[l.rstrip() for l in open('../data/tle/leo.txt') if l.strip()]
TLE={L[k].strip():(L[k+1],L[k+2]) for k in range(len(L)-2) if L[k+1].startswith('1 ') and L[k+2].startswith('2 ')}
def sat_ecef(sat,ts):
    out=[]
    for t in ts:
        d=dtm.datetime.fromtimestamp(t,dtm.timezone.utc); jd,fr=jday(d.year,d.month,d.day,d.hour,d.minute,d.second+d.microsecond/1e6)
        e,r,v=sat.sgp4(jd,fr); g=gstime(jd+fr); out.append([r[0]*np.cos(g)+r[1]*np.sin(g),-r[0]*np.sin(g)+r[1]*np.cos(g),r[2]])
    return np.array(out)
class GeoSat:
    def __init__(s,path,name):
        md,_,_=bwav.read(path); s.t0=md['start_us']/1e6; s.fs=md['samplerate']
        s.TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); s.RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
        s.t=np.arange(-60,360,0.25); s.X=sat_ecef(Satrec.twoline2rv(*TLE[name]),s.t0+s.t)
        s.R1=np.linalg.norm(s.X-s.TX,axis=1); s.R2=np.linalg.norm(s.X-s.RX,axis=1); s.fd=-np.gradient(s.R1+s.R2,s.t)/lam
        s.rr=np.sqrt(s.R1*s.R2); s.name=name
TAUS=np.arange(-15,15.01,0.5); DFS=np.arange(-24,24.01,0.5)
