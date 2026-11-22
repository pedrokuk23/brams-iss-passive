# Distribuição nula: estatística máxima da busca (τ, Δf) com a geometria das passagens da ISS em gravações sem ISS
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
exec(open('synth.py').read().split("if __name__")[0])
passkeys=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
out=[]
for ST in ['BEOTTI','BEOUDS','FRMAUB']:
    geos=[Geo(f'{k}_{ST}') for k in passkeys if os.path.exists(D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav')]
    files=[D+f for f in json.load(open('null_files.json'))['null_iss'][ST]]   # fixed list of recordings without the ISS
    for f in files:
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs)
        for g,k in zip(geos,passkeys):
            z,tau,df=search(bb,ts,fs,fb,g); r=dict(st=ST,file=f.split('/')[-1],geo=k,z=z,tau=tau,df=df); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/null_iss.json','w'))
