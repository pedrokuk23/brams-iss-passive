# Junta os blocos de 1 s de todas as estações, remove blocos com outro avião a <2 Hz no mesmo instante
import json,numpy as np,sys,gzip,os
D=sys.argv[1] if len(sys.argv)>1 else '../results/antmap'
openj=lambda p: json.load(gzip.open(p+'.gz','rt')) if os.path.exists(p+'.gz') else json.load(open(p)); rows=[]
for st in ['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']:
    d=openj(f'{D}/am_{st}.json')
    byf={}
    for tr in d['tracks']: byf.setdefault(tr['file'],[]).append(tr)
    for f,trs in byf.items():
        idx={}
        for j,tr in enumerate(trs):
            for t,fd,r1,r2 in zip(tr['t'],tr['fd'],tr['r1'],tr['r2']): idx.setdefault(t,[]).append((j,fd,20*np.log10(r1*r2)))
        for j,tr in enumerate(trs):
            for i,t in enumerate(tr['t']):
                own=20*np.log10(tr['r1'][i]*tr['r2'][i])
                cont=any(k!=j and abs(fd-tr['fd'][i])<2.5 and pl<own+15 for k,fd,pl in idx[t])
                rows.append((st,f,tr['ac'],tr['type'],tr['tau'],tr['score'],t,tr['snr'][i],tr['fd'][i],tr['r1'][i],tr['r2'][i],
                    tr['az_rx'][i],tr['el_rx'][i],tr['az_tx'][i],tr['el_tx'][i],tr['hdg'][i],tr['beta'][i],int(cont)))
cols='st file ac type tau score t snr fd r1 r2 az_rx el_rx az_tx el_tx hdg beta cont'.split()
json.dump(dict(cols=cols,rows=rows),open(f'{D}/blocks.json','w'))
import collections
c=collections.Counter((r[0],r[-1]) for r in rows); print(c)
