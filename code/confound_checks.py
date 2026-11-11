# Checks of alternative explanations of the sunlight effect:
# (1) orbital state: V2-mini passes binned by the mean altitude of the satellite (from the element set), since
#     satellites raising or lowering their orbits fly other attitudes than those in the operational shells;
# (2) local time: exceedance rate of the sign-reversed nulls of Phase 2 by 3-h UTC window (meteor activity and
#     interference vary with local time; a real change of the false-alarm rate would show here);
# (3) the same shadow/sunlit contrast within each altitude bin, (a) restricted to the near-zenith passes (receiver
#     elevation above 60 deg) and (b) as a Mantel-Haenszel odds ratio stratified by receiver, margin M and elevation.
import json, numpy as np
exec(open('stats_cluster.py').read().split("out=dict(n_boot")[0])       # U2m (Phase 2 own thresholds), sat2, cp
RE=6378.137; MU=398600.4418
def alt_km(s):
    n=s.no_kozai/60.0                                                  # rad/s
    return (MU/n**2)**(1/3)-RE
out={'alt_bins_km':[0,450,500,525,545,700]}
rows=[]
for u in U2m:
    if u['cls']!='Starlink V2-mini': continue
    try: a=alt_km(sat2(u['r']['obj'],u['file']))
    except KeyError: continue
    rows.append(dict(grp=u['grp'],det=bool(u['det']),alt=a,el=u['r'].get('el_rx'),st=u['file'][-6:],M=u.get('M')))
res=[]
for lo,hi in zip(out['alt_bins_km'][:-1],out['alt_bins_km'][1:]):
    d={}
    for g in ('plate_shadow','plate_sunlit'):
        s=[r for r in rows if r['grp']==g and lo<=r['alt']<hi]; k=sum(r['det'] for r in s)
        d[g]=dict(n=len(s),k=int(k),pd=k/len(s) if s else None,ci=cp(k,len(s)) if s else None)
    res.append(dict(lo=lo,hi=hi,**d))
out['v2mini_by_altitude']=res
def mh(rs):
    # Mantel-Haenszel common odds ratio (shadow vs sunlit) over strata receiver x margin bin x elevation bin
    def key(r):
        m=r['M'] if r['M'] is not None else 0; e=r['el'] if r['el'] is not None else 0
        return (r['st'],int(np.digitize(m,[1,10,30])),int(np.digitize(e,[30,60])))
    G={}
    for r in rs: G.setdefault(key(r),[]).append(r)
    num=den=0.0; ns=0; npass=0
    for g in G.values():
        a=sum(r['det'] for r in g if r['grp']=='plate_shadow'); b=sum(not r['det'] for r in g if r['grp']=='plate_shadow')
        c=sum(r['det'] for r in g if r['grp']=='plate_sunlit'); d=sum(not r['det'] for r in g if r['grp']=='plate_sunlit')
        n=a+b+c+d
        if (a+b)==0 or (c+d)==0: continue
        num+=a*d/n; den+=b*c/n; ns+=1; npass+=n
    return dict(or_mh=num/den if den>0 else None,n_strata=ns,n_pass=npass)
res3=[]
for lo,hi in zip(out['alt_bins_km'][:-1],out['alt_bins_km'][1:]):
    s_=[r for r in rows if lo<=r['alt']<hi]
    if not s_: continue
    d={}
    for g in ('plate_shadow','plate_sunlit'):
        z=[r for r in s_ if r['grp']==g and r['el'] is not None and r['el']>60]; k=sum(r['det'] for r in z)
        d[g+'_el_gt60']=dict(n=len(z),k=int(k),pd=k/len(z) if z else None,ci=cp(k,len(z)) if z else None)
    res3.append(dict(lo=lo,hi=hi,**d,mh=mh(s_)))
out['v2mini_by_altitude_adjusted']=res3
S=json.load(open('../results/phase2_summary.json')); th={k:v['z_th'] for k,v in S['thresholds'].items()}
hr={}
for st,f,z,ins,ex in S['nulls']:
    if not ins or ex: continue
    h=int(f[9:11])//3*3; e=hr.setdefault(h,[0,0]); e[0]+=1; e[1]+=int(z>=th[st])
out['null_exceed_by_utc3h']={f'{h:02d}-{h+3:02d}':dict(n=v[0],k=v[1],rate=v[1]/v[0],ci=cp(v[1],v[0])) for h,v in sorted(hr.items())}
json.dump(out,open('../results/confound_checks.json','w'),indent=1)
for r in res: print(r['lo'],r['hi'],{g:(r[g]['k'],r[g]['n']) for g in ('plate_shadow','plate_sunlit')})
for r in res3: print(r)
for k,v in out['null_exceed_by_utc3h'].items(): print(k,v['k'],v['n'])
