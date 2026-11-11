# Headline numbers of the paper, published versus fresh run.
import json,sys,os
A,B=sys.argv[1],sys.argv[2]
L=lambda d,f: json.load(open(os.path.join(d,f)))
def row(name,a,b): print(f'{name[:24]:24s} {a} | {b}')
for lab,f in (('P2',lambda d: L(d,'phase2_summary.json')),):
    a,b=f(A),f(B)
    for c in a['classes']:
        ka,na=a['classes'][c]['all']['k'],a['classes'][c]['all']['n']        # 'all' holds the class total; the
        kb,nb=b['classes'][c]['all']['k'],b['classes'][c]['all']['n']        # other keys are per receiver
        row('P2 '+c,f'{ka}/{na}',f'{kb}/{nb}')
    for s in a['thresholds']: row('z_th '+s,round(a['thresholds'][s]['z_th'],3),round(b['thresholds'][s]['z_th'],3))
    row('P2 detections',len(a['detections']),len(b['detections']))
a,b=L(A,'cluster_stats.json')['feasibility'],L(B,'cluster_stats.json')['feasibility']
for k in ('observed','predicted','pred_95','obs_95'): row('feas '+k,a[k],b[k])
a,b=L(A,'feasibility_validation.json'),L(B,'feasibility_validation.json')
row('M<1 P1',a['impossible']['phase1']['k'],b['impossible']['phase1']['k']); row('M<1 P2',a['impossible']['phase2']['k'],b['impossible']['phase2']['k'])
a,b=L(A,'shadow_strata.json'),L(B,'shadow_strata.json')
for k in a: row('MH '+k,round(a[k]['mh_odds_ratio'],2),round(b[k]['mh_odds_ratio'],2))
a,b=L(A,'null_holdout.json'),L(B,'null_holdout.json')
row('holdout det/n',f"{sum(r['detect'] for r in a)}/{len(a)}",f"{sum(r['detect'] for r in b)}/{len(b)}")
a,b=L(A,'antenna_map_tobit.json'),L(B,'antenna_map_tobit.json')
row('tobit s',round(a['s'],3),round(b['s'],3)); row('tobit beta max|d|','',round(max(abs(x-y) for x,y in zip(a['beta'],b['beta'])),3))
a,b=L(A,'blind_summary.json'),L(B,'blind_summary.json'); row('blind_summary same',a==b,'')
