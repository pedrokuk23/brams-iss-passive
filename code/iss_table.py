import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
exec(open('synth.py').read().split("if __name__")[0])
out=[]
for k in ['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']:
    for ST in ['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']:
        f=D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav'
        if not os.path.exists(f): continue
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); g=Geo(f'{k}_{ST}')
        z,tau,df=search(bb,ts,fs,fb,g); out.append(dict(pass_=k,st=ST,z=z,tau=tau,df=df)); print(k,ST,round(z,1),tau,df,flush=True)
json.dump(out,open('../results/iss_detection_table.json','w'),indent=1)
