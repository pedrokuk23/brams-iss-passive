# Detectability model: bistatic radar equation limited by galactic noise, scattering of simple shapes, glint duration,
# and the link between echo energy and the detection statistic; calibrated and tested with the measured echoes.
import json, numpy as np
c=299792458.; f0=49.97e6; lam=c/f0; kB=1.380649e-23
Pt=120.                                   # W, beacon power at the amplifier output (BRAMS beacon page, 2020 upgrade)
Fa=52-23*np.log10(f0/1e6); Tsky=290*10**(Fa/10)     # ITU-R P.372 median galactic noise: Fa = 52 - 23 log10(f_MHz) dB
K4=(4*np.pi)**3
def cn0(sigma,R,G): return Pt*G*lam**2*sigma/(K4*R**4*kB*Tsky)          # C/N0 in Hz for R = sqrt(R1 R2) in m
def plate_sigma(A): return 4*np.pi*A**2/lam**2
cn=json.load(open('../results/detectability_cn0.json'))
bb={r['label'].split()[-1]:r for r in cn if r['label'].startswith('SPACE')}
sig_bb=plate_sigma(223.)
G={st:10**(r['cn0_peak_dBHz']/10)/cn0(sig_bb,r['rr_at_peak']*1e3,1.) for st,r in bb.items()}
out=dict(Tsky_K=float(Tsky),Fa_dB=float(Fa),Pt_W=Pt,lam_m=lam,sigma_bluebird_m2=float(sig_bb),
         Geff_dB={k:float(10*np.log10(v)) for k,v in G.items()})
# z versus E/N0 for all measured echoes
zz=[(r['z_obs'],r['z_pred'],r['label']) for r in cn]
rat=np.array([a/b for a,b,_ in zz if b>0]); out['z_obs_over_pred']=dict(median=float(np.median(rat)),p10=float(np.percentile(rat,10)),p90=float(np.percentile(rat,90)),n=len(rat))
# effective cross section of the ISS at Oudsbergen and the other receivers, using that receiver's G from the BlueBird
iss=[]
for r in cn:
    if not r['label'].startswith('ISS'): continue
    st=r['label'].split()[-1]
    if st not in G: continue
    s=10**(r['cn0_peak_dBHz']/10)/cn0(1.,r['rr_at_peak']*1e3,G[st]); iss.append(dict(label=r['label'],sigma_m2=float(s),A_eq_m2=float(np.sqrt(s*lam**2/(4*np.pi)))))
out['iss_sigma']=iss
# sensitivity: whole-window statistic needs E/N0 >= z_th*sqrt(nb); a glint of duration T within the window has E = C*T
nb=100; zth=7.; EN0min=zth*np.sqrt(nb); out['EN0_min_window']=float(EN0min)
def sigma_min(R,T,Gst): return EN0min/min(T,nb*0.5)/cn0(1.,R,Gst)
# smallest flat facet, nadir-pointing, near-zenith pass at range h: glint of full width T_g = 2 theta_0 / Omega with
# theta_0 = lam/2L and Omega = v R_E / [h (R_E + h)] (v speed relative to the Earth); its sinc^2 lobe gives the energy
# C_peak * T_eff, T_eff = integral of the lobe over the 50-s window (T_g/2 for a glint much shorter than the window)
from scipy.special import sici
v=7.35e3; RE=6371e3
def Tg_of(h,L): return (lam/L)/(v*RE/(h*(RE+h)))
def Teff_of(Tg,T=nb*0.5):
    a=Tg/2; X=T/(2*a); return 2*a*(sici(2*np.pi*X)[0]/np.pi-np.sin(np.pi*X)**2/(np.pi**2*X))
def Lmin(h,Gst):
    for L in np.arange(0.5,60,0.05):
        if cn0(plate_sigma(L*L),h,Gst)*Teff_of(Tg_of(h,L))>=EN0min: return float(L)
out['Lmin_m']={st:{str(int(h/1e3)):Lmin(h,G[st]) for h in (350e3,420e3,530e3,780e3)} for st in G}
# check of the glint energy on SPACEMOBILE-008: energy-equivalent duration E/C_peak of the echo against T_g/2
out['bluebird_glint']={st:dict(T_eq_s=float(10**((r['E_N0_dB']-r['cn0_peak_dBHz'])/10)),Tg_s=float(Tg_of(r['rr_at_peak']*1e3,np.sqrt(223.))),
    Teff_model_s=float(Teff_of(Tg_of(r['rr_at_peak']*1e3,np.sqrt(223.)),r['nb']*0.5))) for st,r in bb.items()}
# energy-equivalent cross section of the ISS: the echo fluctuates, so its energy over the arc, E/N0, is converted into
# the constant cross section that would give the same energy over a 50-s arc (used for the margin of the feasibility test)
for x in iss:
    r=next(q for q in cn if q['label']==x['label']); x['T_eq_s']=float(10**((r['E_N0_dB']-r['cn0_peak_dBHz'])/10)); x['sigma_energy_50s_m2']=x['sigma_m2']*x['T_eq_s']/50.
out['iss_T_eq_s']=dict(median=float(np.median([10**((r['E_N0_dB']-r['cn0_peak_dBHz'])/10) for r in cn if r['label'].startswith('ISS')])),
    n=sum(1 for r in cn if r['label'].startswith('ISS')))
out['iss_sigma_energy_50s_m2']=dict(median=float(np.median([x['sigma_energy_50s_m2'] for x in iss])),
    p10=float(np.percentile([x['sigma_energy_50s_m2'] for x in iss],10)),p90=float(np.percentile([x['sigma_energy_50s_m2'] for x in iss],90)))
# reference cross sections
a=0.5; ka=2*np.pi*a/lam
out['ref_sigma_m2']=dict(rayleigh_sphere_1m=float(9*np.pi*a**2*ka**4),resonant_wire=float(0.86*lam**2),
     plate_5m=float(plate_sigma(25)),plate_15m=float(plate_sigma(225)),cyl_2m_by_8m=float(2*np.pi*2.1*8**2/lam))
# Envisat upper limits from the glint search (null maximum per receiver as threshold) and the window statistic
gl=json.load(open('../results/glint_envisat.json')); ev=json.load(open('../results/envisat.json'))
ul={}
for st in ['BEOUDS','FRMAUB','BEOTTI']:
    nz=max(r['zg'] for r in gl if r['case']=='null' and st in r['file'])
    passes=[r for r in ev if r['case']=='pass' and r['st']==st]
    Rmin=min(r['rrmin'] for r in passes)*1e3
    ul[st]=dict(z_null_max=float(nz),R_km=float(Rmin/1e3),
        sigma_ul_2s=float(nz*np.sqrt(4)/2.0/cn0(1.,Rmin,G[st])),sigma_ul_16s=float(nz*np.sqrt(32)/16.0/cn0(1.,Rmin,G[st])),
        sigma_ul_window=float(sigma_min(Rmin,50,G[st])))
out['envisat_upper_limits']=ul
json.dump(out,open('../results/detect_theory.json','w'),indent=1)
print(json.dumps({k:v for k,v in out.items() if k!='iss_sigma'},indent=1))
print('ISS sigma range %.0f-%.0f m2'%(min(r['sigma_m2'] for r in iss),max(r['sigma_m2'] for r in iss)))
for r in iss: print(r['label'],'%.0f'%r['sigma_m2'])
