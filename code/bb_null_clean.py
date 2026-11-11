# Distribuição nula da busca BlueBird (grade ampliada) nos receptores SDR, com remoção das linhas CW
exec(open('bb2.py').read())
from linefilt import clean
G=[('20261105_0310','SPACEMOBILE-008'),('20261104_0015','SPACEMOBILE-009'),('20261104_0040','SPACEMOBILE-010')]
passkeys=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
import os
DL=os.environ.get('BBDIR',D); out=[]
for ST in ['BEOUDS','FRMAUB']:
    geos=[(GeoSat(DL+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav',n),k,n) for k,n in G]
    files=[D+f for f in json.load(open('bluebird_files.json'))['null'][ST]]   # fixed list: recordings without BlueBird or ISS passes
    for f in files:
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=clean(bb,ts,fs)
        for g,k,n in geos:
            z,tau,df=search(bb,ts,fs,fb,g,TAUS,DFS); r=dict(st=ST,file=f.split('/')[-1],geo=k,sat=n,z=z,tau=tau,df=df,clean=True); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/bb_null_clean.json','w'))
