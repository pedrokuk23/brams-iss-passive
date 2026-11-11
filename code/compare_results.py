# Fresh run versus published results, tolerant of record order and of last-digit floating-point differences.
# Lists of records (dicts) are matched by their string-valued fields (file, object, station...); numbers are compared by
# relative difference. For each file: matched values, share within 1e-6 / 1e-3 / 1e-2, the largest difference, and
# records present in only one of the two runs.
import json, os, sys, glob, math
A,B=sys.argv[1],sys.argv[2]
def key(d): return tuple(sorted((k,v) for k,v in d.items() if isinstance(v,(str,bool)) or v is None))
def flat(x,p=''):
    if isinstance(x,dict):
        for k,v in x.items(): yield from flat(v,f'{p}/{k}')
    elif isinstance(x,list):
        if x and all(isinstance(e,dict) for e in x):
            seen={}
            for e in x:
                k=key(e); seen[k]=seen.get(k,-1)+1
                yield from flat(e,f'{p}[{"|".join(str(v) for _,v in k)}#{seen[k]}]')
        else:
            for i,v in enumerate(x): yield from flat(v,f'{p}[{i}]')
    else: yield p,x
rows=[]; missing=[]
for fa in sorted(glob.glob(A+'/**/*.json',recursive=True)):
    rel=os.path.relpath(fa,A); fb=os.path.join(B,rel)
    if not os.path.exists(fb): missing.append(rel); continue
    if os.path.getsize(fa)>float(os.environ.get('MAXB',3e8)): rows.append((rel,'skipped (size %.0f MB)'%(os.path.getsize(fa)/1e6))); continue
    try: a=dict(flat(json.load(open(fa)))); b=dict(flat(json.load(open(fb))))
    except Exception as e: rows.append((rel,'unreadable: '+str(e)[:50])); continue
    common=set(a)&set(b); r=[]; other=0
    for k in common:
        x,y=a[k],b[k]
        if isinstance(x,bool) or isinstance(y,bool) or not isinstance(x,(int,float)) or not isinstance(y,(int,float)):
            other+= x!=y; continue
        if isinstance(x,float) and math.isnan(x) and isinstance(y,float) and math.isnan(y): r.append((0,k,x,y)); continue
        r.append((abs(x-y)/max(abs(x),abs(y),1e-9),k,x,y))
    n=len(r); f=lambda t: sum(1 for d in r if d[0]<=t)
    w=max(r) if r else (0,'',0,0)
    rows.append((rel,f'{n} num: {f(1e-6)} <1e-6, {f(1e-3)} <1e-3, {f(1e-2)} <1e-2; max {w[0]:.1e} at {w[1][:60]} ({w[2]} -> {w[3]}); text diffs {other}; only published {len(set(a)-set(b))}, only new {len(set(b)-set(a))}'))
o=open(os.path.join(B,'..','compare_report2.md'),'w')
o.write('# Fresh run versus published results (matched by record)\n\n')
for r in rows: o.write(f'- **{r[0]}**: {r[1]}\n')
o.write(f'\nPublished files not produced by the fresh run ({len(missing)}): '+', '.join(missing)+'\n'); o.close()
print(open(os.path.join(B,'..','compare_report2.md')).read())
