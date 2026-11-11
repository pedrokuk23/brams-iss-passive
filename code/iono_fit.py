import os; WAVDIR_=os.environ.get('WAVDIR','../data/raw_wav/')   # folder of the Phase 1 recordings
import numpy as np, json, sys
from sgp4.propagation import gstime
sys.argv=['x','X']
base=open('station_run.py').read().split('print(f"== {ST} ==")')[0]
exec(base)
import bwav
D=WAVDIR_
tr=json.load(open('../results/iono_tracks.json'))
c=299792458.0; f0=49.97e6; K=40.3/(c*f0)   # Hz por (el/m^2/s)
RE=6371.0
def chapman(h,Nm,hm=300.,H=55.):
    z=(h-hm)/H; return Nm*np.exp(0.5*(1-z-np.exp(-z)))
def stec(A,B,Nm):
    s=np.linspace(0,1,400)[:,None]; P=A[None,:]+s*(B-A)[None,:]
    h=np.linalg.norm(P,axis=1)-RE; L=np.linalg.norm(B-A)*1e3
    return np.trapezoid(chapman(h,Nm),dx=1/399)*L
def cleaned(a):
    a=np.array(a); med=np.median(a[:,1]); return a[abs(a[:,1]-med)<4]
passes={'0220 03/11':['20261103_0220_BEOTTI','20261103_0220_BEOUDS'],'2355 03/11':['20261103_2355_BEOTTI','20261103_2355_BEOUDS'],
        '0040 03/11':['20261103_0040_BEOTTI'],'0130 04/11':['20261104_0130_BEOUDS'],'2305 04/11':['20261104_2305_BEOUDS']}
res={}
for pname,keys in passes.items():
    obs=[];mod_t=[];mod_i=[]
    for key in keys:
        if key not in tr or len(tr[key])<5: continue
        st=key.split('_')[2]; fn=D+f'RAD_BEDOUR_{key}_SYS001.wav'
        md,x,pps=bwav.read(fn); t0=md['start_us']/1e6
        TX=ecef(md['beacon_lat'],md['beacon_long'],md['beacon_alt']/1000); RX=ecef(md['ant_lat'],md['ant_long'],md['ant_alt']/1000)
        a=cleaned(tr[key]); tc=a[:,0]
        t=np.arange(tc.min()-5,tc.max()+5,0.5); X=iss_ecef(t0+t)
        Rb=np.linalg.norm(X-TX,axis=1)+np.linalg.norm(X-RX,axis=1); fd=-np.gradient(Rb,t)/lam; dfd=np.gradient(fd,t)
        S=np.array([stec(TX,p,1.0)+stec(p,RX,1.0) for p in X])   # por unidade de Nm (m^-3)
        dS=np.gradient(S,t)
        obs+=list(a[:,1]); mod_t+=list(-np.interp(tc,t,dfd)); mod_i+=list(K*np.interp(tc,t,dS))
    obs=np.array(obs); mt=np.array(mod_t); mi=np.array(mod_i)
    def lsq(cols):
        A=np.column_stack(cols); p,*_=np.linalg.lstsq(A,obs,rcond=None); r=obs-A@p; return p,np.sqrt(np.mean(r**2))
    pc,rc=lsq([np.ones_like(obs)])
    pt,rt=lsq([mt])
    pi,ri=lsq([mi])
    pti,rti=lsq([mt,mi])
    res[pname]=dict(n=len(obs),rms_const=rc,const=pc[0],rms_tempo=rt,dt_s=pt[0],dt_km=pt[0]*7.66,rms_iono=ri,NmF2=pi[0],foF2=np.sqrt(abs(pi[0])/1.24e10)/1e6*np.sign(pi[0]),rms_ambos=rti)
    r=res[pname]
    print(f"{pname}: n={r['n']:3d} | const {r['const']:+.2f} Hz rms {r['rms_const']:.2f} | tempo δt={r['dt_s']*1000:+.0f} ms ({r['dt_km']:+.2f} km) rms {r['rms_tempo']:.2f} | iono NmF2={r['NmF2']:.2e} (foF2={r['foF2']:+.1f} MHz) rms {r['rms_iono']:.2f} | tempo+iono rms {r['rms_ambos']:.2f}")
json.dump(res,open('../results/iono_fit.json','w'),indent=1)
