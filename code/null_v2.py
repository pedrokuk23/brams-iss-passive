# Distribuição nula com remoção de linhas CW nos receptores SDR: busca (τ, Δf) com a geometria das passagens
# da ISS aplicada a gravações sem ISS. Uso: STS="BEOUDS FRMAUB" OUT=... python3 null_v2.py
exec(open('synth.py').read().split("if __name__")[0])
from linefilt import maybe_clean
passkeys=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
DIRS=[D,os.environ.get('DAYDIR',D)]   # gravações noturnas e amostra diurna de 04/11 (por padrão, ambas em ../data/raw_wav/)
DAY=['0100','0300','0600','0900','1200','1500','1800','2100']      # amostra do dia 04/11
out=[]
for ST in os.environ['STS'].split():
    geos=[(k,Geo(f'{k}_{ST}')) for k in passkeys if os.path.exists(D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav')]
    files=[(DIRS[1] if '_20261104_' in f and f[20:24] in DAY else DIRS[0])+f for f in json.load(open('null_files.json'))['null_v2'][ST]]   # fixed list: night recordings without the ISS and the day sample of 4 Nov
    for f in files:
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=maybe_clean(bb,ts,fs)
        for k,g in geos:
            z,tau,df=search(bb,ts,fs,fb,g); r=dict(st=ST,file=f.split('/')[-1],geo=k,z=z,tau=tau,df=df,clean=fs>6000); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open(os.environ['OUT'],'w'))
