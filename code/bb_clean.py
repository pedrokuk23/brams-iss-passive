# BlueBird: repete a busca de grade ampliada nos receptores SDR depois da remoção das linhas CW
exec(open('bb2.py').read())
import os
from linefilt import clean
prev={(r['file'],r['sat']):r for r in json.load(open('../results/bb2_results.json'))}
out=[]
for (f,n),r0 in sorted(prev.items()):
    if r0['rrmin']>600 or 'BEOTTI' in f: continue
    p=os.environ.get('BBDIR',D)+f
    md,x,_=bwav.read(p); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); bb=clean(bb,ts,fs)
    g=GeoSat(p,n); z,tau,df=search(bb,ts,fs,fb,g,TAUS,DFS)
    r=dict(file=f,sat=n,rrmin=r0['rrmin'],z_raw=r0['z'],tau_raw=r0['tau'],df_raw=r0['df'],z=z,tau=tau,df=df); out.append(r); print(json.dumps(r),flush=True)
json.dump(out,open('../results/bb_clean.json','w'),indent=1)
