# Phase 2 of the blind search with Maubeuge reported separately (Sec. 10.2). The nulls of Maubeuge are not
# exchangeable between its recordings (its out-of-sample false-alarm rate is 1.6%), so its detections are kept out of
# the class statistics and the pooled false-alarm rate is computed over the other four receivers with the same
# split-half procedure (far_checks.py, same seed). Reads results/phase2_summary.json and phase2_byfile.json; writes
# results/maubeuge_separate.json.
import json, collections as C
from far_checks import far_checks, cp
S=json.load(open('../results/phase2_summary.json')); B=json.load(open('../results/phase2_byfile.json'))
M='FRMAUB'; STS=[s for s in ['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE'] if s!=M]
N=[dict(file=f,z=z,_in=i,_ex=e) for st,f,z,i,e in S['nulls']]
expl={id(r):r['_ex'] for r in N}; th={k:v['z_th'] for k,v in S['thresholds'].items()}
fc=far_checks(N,expl,th,STS,lambda r:r['_in'])
det=[d for d in S['detections'] if not d['file'].endswith(M)]; dm=[d for d in S['detections'] if d['file'].endswith(M)]
units=C.Counter(); k=C.Counter()
issrec=set()
for f,L in B.items():
    if f.endswith(M): continue
    for p in L:
        if not p.get('possible',True): continue
        if p['cls']=='ISS complex': issrec.add(f)
        else: units[p['cls']]+=1
units['ISS complex']=len(issrec)
for d in det:
    k[d['class']]+=0 if d['class']=='ISS complex' else 1
k['ISS complex']=len({d['file'] for d in det if d['class']=='ISS complex'})
classes={c:dict(n=units[c],k=k[c],pd=k[c]/units[c] if units[c] else None,ci=cp(k[c],units[c]) if units[c] else None) for c in units}
out=dict(far_pooled_4rx=fc['split_half']['all'],far_by_rx={s:fc['split_half'][s] for s in STS},classes_without_maubeuge=classes,
         maubeuge=dict(n_detections=len(dm),by_class=dict(C.Counter(d['class'] for d in dm)),
                       n_recordings=len({f for f in B if f.endswith(M)})))
json.dump(out,open('../results/maubeuge_separate.json','w'),indent=1); print(json.dumps(out,indent=1))
