# Estatística de detecção da ISS em cada passagem e receptor, sem e com remoção das linhas CW (SDR)
exec(open('synth.py').read().split("if __name__")[0])
from linefilt import lines, clean
out=[]
for k in ['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']:
    for ST in ['BEOTTI','BEOUDS','FRMAUB','BEHUMA','BEUCCE']:
        f=D+f'RAD_BEDOUR_{k}_{ST}_SYS001.wav'
        if not os.path.exists(f): continue
        md,x,_=bwav.read(f); fs=md['samplerate']; fb,bb,ts=baseband(x,fs); g=Geo(f'{k}_{ST}')
        z0,t0_,d0=search(bb,ts,fs,fb,g)
        fl,lv=lines(bb,fs); top=sorted(zip(lv,fl))[::-1][:3]
        s0=np.sqrt(np.median(np.abs(bb)**2)/np.log(2))
        if fs>6000:
            bb2=clean(bb,ts,fs,fl); z1,t1,d1=search(bb2,ts,fs,fb,g); s1=np.sqrt(np.median(np.abs(bb2)**2)/np.log(2))
        else: z1,t1,d1,s1=z0,t0_,d0,s0
        r=dict(pass_=k,st=ST,fs=fs,fb=fb,z_raw=z0,tau_raw=t0_,df_raw=d0,z=z1,tau=t1,df=d1,noise_gain_db=20*np.log10(s0/s1),
               strongest_lines=[(round(b,1),round(a,1)) for a,b in top],nlines=len(fl))
        out.append(r); print(k,ST,round(z0,1),t0_,'->',round(z1,1),t1,d1,round(r['noise_gain_db'],1),r['strongest_lines'][:2],flush=True)
json.dump(out,open('../results/iss_detection_table_v2.json','w'),indent=1)
