# Earth shadow versus sunlight, stratified (exploratory, PROTOCOL item 8): for each class of flat-array satellite, the
# detection probability in shadow and in sunlight within strata of receiver, margin M, elevation of the pass at the
# receiver and local hour, so that the contrast cannot come from a different mix of these factors. Intervals: bootstrap
# over recordings. Mantel-Haenszel common odds ratio over the strata of receiver x margin x elevation (as before) and
# x hour (3-h UTC bins), with its Robins-Breslow-Greenland interval, a bootstrap over satellites and the Breslow-Day
# test of homogeneity; odds ratios by elevation and by hour.
import json, collections, numpy as np
exec(open('stats_cluster.py').read().split("out=dict(n_boot")[0])     # U1, U2m with grp (shadow/sunlit), M, r (el_rx)
def strat(U,key):
    g=collections.defaultdict(lambda:{'plate_shadow':[],'plate_sunlit':[]})
    for u in U: g[key(u)][u['grp']].append(u)
    return g
def pr(L):
    if not L: return None
    k=sum(u['det'] for u in L); return dict(n=len(L),k=int(k),pd=k/len(L),ci=boot_rate(L) if len(L)>=20 else None)
def mbin(M): return '<10' if M<10 else '10-30' if M<30 else '>=30'
def ebin(e): return '<30' if e<30 else '30-60' if e<60 else '>=60'
from statsmodels.stats.contingency_tables import StratifiedTable
def tables(V,key):
    T=[]
    for k,v in strat(V,key).items():
        a=sum(u['det'] for u in v['plate_shadow']); b=len(v['plate_shadow'])-a; c=sum(u['det'] for u in v['plate_sunlit']); dd=len(v['plate_sunlit'])-c
        if (a+b) and (c+dd) and (a+c) and (b+dd): T.append(np.array([[a,b],[c,dd]],float))
    return T
def mh(V,key):
    T=tables(V,key)
    if not T: return None
    st=StratifiedTable(T); n=int(sum(t.sum() for t in T))
    try: bd=st.test_equal_odds(adjust=True); bdp=float(bd.pvalue)
    except Exception: bdp=None
    try: lo,hi=st.oddsratio_pooled_confint()
    except Exception: lo=hi=None
    return dict(or_mh=float(st.oddsratio_pooled),ci=[float(lo),float(hi)] if lo is not None else None,breslow_day_p=bdp,n_strata=len(T),n=n,
                k_shadow=int(sum(t[0,0] for t in T)),n_shadow=int(sum(t[0].sum() for t in T)),k_sunlit=int(sum(t[1,0] for t in T)),n_sunlit=int(sum(t[1].sum() for t in T)))
def boot_mh(V,key,nb=300):
    g=collections.defaultdict(list)
    for u in V: g[u['r']['obj']].append(u)
    G=list(g.values()); r=[]
    for _ in range(nb):
        W=[u for i in rng.integers(0,len(G),len(G)) for u in G[i]]
        try: r.append(StratifiedTable(tables(W,key)).oddsratio_pooled)
        except Exception: pass
    r=np.array([x for x in r if np.isfinite(x)])
    return [float(x) for x in np.percentile(r,[2.5,97.5])] if len(r) else None
out={}
for ph,U in (('phase1',U1),('phase2',U2m)):
    for c in ('Starlink V2-mini','Starlink DTC','AST SpaceMobile','Starlink V1.x'):
        V=[u for u in U if u['cls']==c]; d={}
        for name,key in (('receiver',lambda u:u['file'][-6:]),('margin',lambda u:mbin(u['M'])),('elevation',lambda u:ebin(u['r']['el_rx'])),
                         ('hour_utc',lambda u:'%02d-%02d'%(int(u['file'][9:11])//3*3,int(u['file'][9:11])//3*3+3))):
            d[name]={k:{g:pr(v[g]) for g in v} for k,v in sorted(strat(V,key).items())}
        # Mantel-Haenszel odds ratio shadow vs sunlit over receiver x margin x elevation
        num=den=0.0
        for k,v in strat(V,lambda u:(u['file'][-6:],mbin(u['M']),ebin(u['r']['el_rx']))).items():
            a=sum(u['det'] for u in v['plate_shadow']); b=len(v['plate_shadow'])-a; cc=sum(u['det'] for u in v['plate_sunlit']); dd=len(v['plate_sunlit'])-cc; n=a+b+cc+dd
            if n and v['plate_shadow'] and v['plate_sunlit']: num+=a*dd/n; den+=b*cc/n
        d['mh_odds_ratio']=num/den if den else None
        # with the hour of the day (3-h UTC bins, nearly the local sidereal time bins over these two weeks, which also
        # sets the sky noise) added to the strata; Breslow-Day test of a common odds ratio; odds ratio by elevation and by
        # hour; bootstrap over satellites (all passes of one satellite resampled together) for the confidence interval
        hk=lambda u:int(u['file'][9:11])//3
        key4=lambda u:(u['file'][-6:],mbin(u['M']),ebin(u['r']['el_rx']),hk(u))
        d['mh_with_hour']=mh(V,key4)
        d['mh_03_06_utc']=mh([u for u in V if hk(u)==1],lambda u:(u['file'][-6:],mbin(u['M']),ebin(u['r']['el_rx'])))
        d['or_by_elevation']={e:mh([u for u in V if ebin(u['r']['el_rx'])==e],lambda u:(u['file'][-6:],mbin(u['M']),hk(u))) for e in ('<30','30-60','>=60')}
        d['or_by_hour']={'%02d-%02d'%(3*h,3*h+3):mh([u for u in V if hk(u)==h],lambda u:(u['file'][-6:],mbin(u['M']),ebin(u['r']['el_rx']))) for h in range(8)}
        if d['mh_with_hour']: d['mh_with_hour']['ci_boot_satellite']=boot_mh(V,key4)
        out[f'{ph}|{c}']=d
json.dump(out,open('../results/shadow_strata.json','w'),indent=1)
for k,d in out.items():
    print('==',k,'MH OR %.1f'%d['mh_odds_ratio'] if d.get('mh_odds_ratio') else '','| with hour',d['mh_with_hour'],'| 03-06',d['mh_03_06_utc'])
    print('   by elevation',{e:(round(v['or_mh'],2),[round(x,2) for x in v['ci']] if v and v['ci'] else None,v['breslow_day_p']) if v else None for e,v in d['or_by_elevation'].items()})
    if 'V2-mini' in k or 'DTC' in k:
        for name in ('receiver','margin','elevation'):
            for s,v in d[name].items():
                a,b=v.get('plate_shadow'),v.get('plate_sunlit')
                f=lambda x:'%d/%d=%.3f'%(x['k'],x['n'],x['pd']) if x else '-'
                print('   %-9s %-7s shadow %-18s sunlit %s'%(name,s,f(a),f(b)))
