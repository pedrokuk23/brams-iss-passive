# Estatísticas da calibração com aviões (Tabela II, Fig. 3) a partir de uni_*.json (saída de air_uniform.py)
import json,glob,numpy as np,collections
Z=6; rows=[]
for f in sorted(glob.glob('../results/uni_*.json')): rows+=json.load(open(f))
order=[('20261103','BEOTTI'),('20261103','BEOUDS'),('20261103','FRMAUB'),('20261103','BEHUMA'),('20261104','BEOTTI'),('20261104','BEOUDS'),('20261104','FRMAUB'),('20261104','BEHUMA'),('20261104','BEUCCE')]
tab={}; det_all=[]; by=collections.defaultdict(set)
for d,st in order:
    R=[r for r in rows if r['day']==d and r['st']==st]
    W=[r for r in R if abs(r['tau'])<=1]; D=[r['tau_ref'] for r in W if r['z']>=Z]
    for r in W:
        if r['z']>=Z: by[(d,r['ac'])].add(st)
    med=np.median(D); mad=1.4826*np.median(np.abs(np.array(D)-med))
    tab[f'{d}|{st}']=dict(N=len(R),win=len(W),det=len(D),med_ms=med*1e3,mad_ms=mad*1e3,se_ms=1.253*mad/np.sqrt(len(D))*1e3); det_all+=D
det_all=np.array(det_all); m=np.median(det_all)
out=dict(N=len(rows),win=sum(v['win'] for v in tab.values()),chance=2/18,Z=Z,tab=tab,n_det=len(det_all),med_all_ms=m*1e3,
         mad_all_ms=1.4826*np.median(np.abs(det_all-m))*1e3,unique=len(by),by_count={str(k):v for k,v in collections.Counter(len(s) for s in by.values()).items()})
json.dump(out,open('../results/aircraft_uniform.json','w'),indent=1); print(json.dumps({k:v for k,v in out.items() if k!='tab'}))
