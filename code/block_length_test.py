# Sensibilidade da estatística de detecção da ISS à duração do bloco coerente T_b
exec(open('synth.py').read().split("if __name__")[0])
from linefilt import maybe_clean
out=[]
for k in ['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']:
    for ST in ['BEOTTI','BEOUDS','BEHUMA']:
        f=D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav'
        if not os.path.exists(f): continue
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs); g=Geo(f'{k}_{ST}')
        r=dict(pass_=k,st=ST)
        for Tb in [0.25,0.5,1.0,2.0]:
            z,tau,nu=search(bb,ts,fs,fb,g,Tb=Tb); r[str(Tb)]=[z,tau,nu]
        out.append(r); print(k,ST,{t:round(r[t][0],1) for t in ['0.25','0.5','1.0','2.0']},flush=True)
json.dump(out,open('../results/block_length_test.json','w'),indent=1)
