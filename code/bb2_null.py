import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
exec(open('bb2.py').read())
# geometria: passagem do SPACEMOBILE-008 de 05/11 03:10 e do -009 de 04/11 00:15 e do -010 de 04/11 00:40 em cada estação
G=[('20261105_0310','SPACEMOBILE-008'),('20261104_0015','SPACEMOBILE-009'),('20261104_0040','SPACEMOBILE-010')]
passkeys=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
out=[]
for ST in ['BEOTTI','BEOUDS','FRMAUB']:
    geos=[(GeoSat(f'{WAVDIR_}RAD_BEDOUR_{k}_{ST}_SYS001.wav',n),k,n) for k,n in G]
    files=[D+f for f in json.load(open('bluebird_files.json'))['null'][ST]]   # fixed list: recordings without BlueBird or ISS passes
    for f in files:
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs)
        for g,k,n in geos:
            z,tau,df=search(bb,ts,fs,fb,g,TAUS,DFS); r=dict(st=ST,file=f.split('/')[-1],geo=k,sat=n,z=z,tau=tau,df=df); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/bb2_null.json','w'))
