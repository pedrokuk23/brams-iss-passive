# Check of the manuscript against results/paper_numbers.json: every registered text form must appear in the indicated
# part of the manuscript (main text, abstract, supplementary), and no superseded value may remain. The manuscript
# sources are read from PAPERDIR (default ../paper/); without them the check is skipped.
import json, os, re, sys
P=os.environ.get('PAPERDIR','../paper/')
main=os.path.join(P,'brams_asr.tex'); supp=os.path.join(P,'brams_asr_supplementary.tex')
if not os.path.exists(main): print('check_numbers: manuscript sources not found in',P,'- skipped'); sys.exit(0)
T={'main':open(main).read()}
T['supp']=open(supp).read() if os.path.exists(supp) else ''
m=re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}',T['main'],re.S); T['abstract']=m.group(1) if m else ''
def norm(s): return re.sub(r'\\amarelo\{|\\begin\{novo\}|\\end\{novo\}','',s)
N={k:norm(v) for k,v in T.items()}
R=json.load(open('../results/paper_numbers.json')); bad=0
for e in R['entries']:
    for w in e['where']:
        if e['text'] not in N[w]:
            print('MISSING in %-8s %-20s %s'%(w,e['key'],e['text'])); bad+=1
for s in R['stale']:
    for w in ('main','supp'):
        if s in T[w]: print('STALE   in %-8s %s'%(w,s)); bad+=1
print('check_numbers: %d entries, %d problems'%(R['n'],bad)); sys.exit(1 if bad else 0)
