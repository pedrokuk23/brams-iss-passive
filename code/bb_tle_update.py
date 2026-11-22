# BlueBird SPACEMOBILE-008: deslocamento ao longo da órbita, no instante do eco (5/11 03:11 UTC), entre o elemento de
# 3/11 usado na busca e os elementos emitidos depois (gp_history), expresso como δt (δt>0: satélite atrás do previsto).
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
from sgp4.api import Satrec, jday
import numpy as np, json
L=[l.rstrip() for l in open('../data/tle/bluebird_gp_history.txt') if l.strip()]
S={}
for i in range(0,len(L),3):
    s=Satrec.twoline2rv(L[i+1],L[i+2]); S.setdefault(L[i][2:].strip(),{})[s.jdsatepoch+s.jdsatepochF]=s
sat=S['SPACEMOBILE-008']; eps=sorted(sat); old=sat[eps[0]]
jd=sum(jday(2026,11,5,3,11,0))
e,r0,v0=old.sgp4(int(jd),jd-int(jd)); r0=np.array(r0); v0=np.array(v0); vh=v0/np.linalg.norm(v0)
out=[]
for ep in eps[1:]:
    e,r1,_=sat[ep].sgp4(int(jd),jd-int(jd)); dr=np.array(r1)-r0; al=float(dr@vh)
    out.append(dict(epoch_jd=ep,hours_from_echo=(ep-jd)*24,along_km=al,dt_ms=-al/np.linalg.norm(v0)*1e3,other_km=float(np.linalg.norm(dr-al*vh))))
    print('%+7.1f h  along %+.3f km  dt %+6.1f ms'%(out[-1]['hours_from_echo'],al,out[-1]['dt_ms']))
json.dump(dict(reference_epoch_jd=eps[0],echo_utc='2026-11-05T03:11:00',newer=out),open('../results/bb_tle_update.json','w'),indent=1)
