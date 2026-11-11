import json,glob,numpy as np
R=[]
for f in sorted(glob.glob('../results/catscan/cat3_*.json')): R+=json.load(open(f))
files=sorted(set(r['file'] for r in R))
for fk in files:
    S=[r for r in R if r['file']==fk]; iss=[r for r in S if r['obj']=='ISS (ZARYA)']
    for r in S: r['issc']=bool(iss) and abs(r['rr']-iss[0]['rr'])<2 and abs(r['peak']-iss[0]['peak'])<1e-6
C=[r for r in R if not r['issc']]; I=[r for r in R if r['obj']=='ISS (ZARYA)']
p=np.array([r['peak'] for r in C]); c=np.array([r['ctrl'] for r in C])
print('pairs',len(R),'files',len(files),'non-ISS pairs',len(C),'objects',len(set(r['obj'] for r in C)))
print('ctrl quantiles 50/99/99.9/max',np.percentile(c,[50,99,99.9]),c.max())
thr=c.max()
for th in [3,4,5,thr]: print(f'th {th:.2f}: peak>th {np.sum(p>th)} ctrl>th {np.sum(c>th)}')
print('ISS:'); 
for r in I: print(' ',r['file'],round(r['peak'],1),round(r['ctrl'],1),r['df'],round(r['rr']))
for r in sorted(C,key=lambda r:-r['peak'])[:20]: print(f"{r['file']} {r['obj'][:24]:24s} peak={r['peak']:5.1f} ctrl={r['ctrl']:4.1f} df={r['df']:+5.1f} rr={r['rr']:6.0f} nb={r['nb']}")
json.dump(R,open('../results/catalog_scan.json','w'))   # list of (file, object) records, read by close_scan.py and glint_close.py
