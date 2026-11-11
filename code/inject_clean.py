# Injection test without the four recordings whose track was not clean: two tracks of ISS modules
# (ZVEZDA, POISK), which carry the station's own echo, and two Starlink tracks at Oudsbergen with a real echo already
# detected on them. Reads results/inject_test.json, writes results/inject_test_clean.json.
import json, numpy as np
d=json.load(open('../results/inject_test.json'))
BAD={'20261020_0345_BEHUMA','20261022_0835_BEOUDS','20261021_1925_BEOUDS','20261026_0730_BEOUDS'}
out={'excluded':sorted(BAD),'summary':{}}
for env in ('steady','glint'):
    for lv in sorted({x['level'] for x in d['rows']}):
        r=[x for x in d['rows'] if x['env']==env and x['level']==lv and x['file'] not in BAD]
        out['summary'].setdefault(env,{})[str(lv)]=dict(n=len(r),pd=float(np.mean([x['det'] for x in r])),
            median_z_over_pred=float(np.median([x['z']/x['z_pred'] for x in r])))
json.dump(out,open('../results/inject_test_clean.json','w'),indent=1)
for env,v in out['summary'].items(): print(env,{k:round(w['pd'],2) for k,w in v.items()})
