# Validação independente do critério de detecção (congelado antes desta execução):
# z >= 7, |tau| <= 0.5 s e pico fora da borda da grade de frequência (|Δf| < 12 Hz).
# Arquivos: amostra de 30 em 30 min de 04/11 (05:00-21:00 UTC) que NÃO entrou em null_iss.json nem em null_v2.json.
import os as _os; _os.environ['BRAMS_TX_HEADER']='1'   # historical configuration: transmitter position written in each file, as in the original run (README, Configurations)
exec(open('synth.py').read().split("if __name__")[0])
from linefilt import maybe_clean
passkeys=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
DAYDIR=os.environ.get('DAYDIR',D)
used=set(r['file'] for r in json.load(open('../results/null_v2.json')))|set(r['file'] for r in json.load(open('../results/null_iss.json')))
out=[]
for ST in ['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']:
    geos=[(k,Geo(f'{k}_{ST}')) for k in passkeys if os.path.exists(D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav')]
    nm=lambda f: f.split('/')[-1]
    files=[f for f in sorted(glob.glob(DAYDIR+f'RAD_BEDOUR_20261104_*_{ST}_SYS001.wav'))
           if nm(f)[22:24] in ('00','30') and 5<=int(nm(f)[20:22])<=21 and nm(f) not in used]
    print('#',ST,len(files),'arquivos',flush=True)
    for f in files:
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
        for k,g in geos:
            z,tau,df=search(bb,ts,fs,fb,g)
            r=dict(st=ST,file=f.split('/')[-1],geo=k,z=z,tau=tau,df=df,detect=bool(z>=7 and abs(tau)<=0.5 and abs(df)<12)); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/null_holdout.json','w'))
