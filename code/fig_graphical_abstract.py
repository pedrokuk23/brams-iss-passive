# Graphical abstract (Elsevier, about 13:5): network, blind search per class against false alarms, V2-mini sunlight
# effect, ISS timing against NASA predictions (Phase 1 pairs and Phase 2 passes with predictions < 1 day old)
import numpy as np, sys, json, glob, os
sys.argv=['x','X']
exec(open('station_run.py').read().split('print(f"== {ST} ==")')[0])
import bwav
exec(open('figstyle.py').read())
D=os.environ.get('WAVDIR','../data/raw_wav/')
C=json.load(open('../data/stations/coords.json')); TXll=C['tx']
fig,ax=plt.subplots(1,4,figsize=(11.6,3.4),gridspec_kw=dict(width_ratios=[1,1.1,0.6,1.0],wspace=0.62))
a=ax[0]
for key in ['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']:
    import datetime as _dt; t0=_dt.datetime.strptime(key,'%Y%m%d_%H%M').replace(tzinfo=_dt.timezone.utc).timestamp()   # file start (ground track only)
    X=iss_ecef(t0+np.arange(-200,500,10.)); r=np.linalg.norm(X,axis=1)
    a.plot(np.degrees(np.arctan2(X[:,1],X[:,0])),np.degrees(np.arcsin(X[:,2]/r)),color='#9a9a95',lw=0.7,zorder=1)
a.plot(TXll[1],TXll[0],'*',ms=10,color='#222',zorder=3)
a.annotate('TX',(TXll[1],TXll[0]),xytext=(4,-9),textcoords='offset points',fontsize=7)
for st,(la,lo) in C['stations'].items():
    a.plot(lo,la,MK[st],ms=5,color=COL[st],mec='white',mew=0.6,zorder=3)
    off={'BEUCCE':(-30,4),'FRMAUB':(-42,-2),'BEOTTI':(4,-3),'BEOUDS':(-20,5),'BEHUMA':(4,-6)}[st]
    a.annotate(NAME[st],(lo,la),xytext=off,textcoords='offset points',fontsize=7)
a.set_xlim(2.5,7); a.set_ylim(49.3,51.7); a.set_aspect(1/np.cos(np.radians(50.3)))
a.set_title('Open meteor radar: 49.97-MHz\nbeacon, 5 receivers, ISS tracks',fontsize=7.5); a.tick_params(labelsize=7)
a.set_xlabel('Longitude (°E)',fontsize=7); a.set_ylabel('Latitude (°N)',fontsize=7)
for s in ['top','right']: a.spines[s].set_visible(False)
b=ax[1]; P2=json.load(open('../results/phase2_summary.json'))
cls=[('ISS complex','ISS'),('AST SpaceMobile','AST BlueBird'),('Starlink DTC','Starlink DTC'),('Starlink V2-mini','Starlink V2-mini'),
     ('Starlink V1.x','Starlink V1.x'),('Other','Other objects')]
for i,(c,lab) in enumerate(cls):
    v=P2['classes'][c]['all']; lo,hi=v['ci']
    b.errorbar(100*v['pd'],i,xerr=[[100*(v['pd']-lo)],[100*(hi-v['pd'])]],fmt='o',color='#2b5d8a',ms=4,lw=1,capsize=2)
    b.annotate(f"{v['k']:,}/{v['n']:,}",(100*hi,i),xytext=(4,-3),textcoords='offset points',fontsize=6,color='#555')
b.axvspan(0.12,0.16,color='#c0392b',alpha=0.35,lw=0); b.annotate('false alarms\n0.12–0.16%',(0.14,len(cls)-0.4),ha='center',fontsize=6,color='#c0392b')
b.set_xscale('log'); b.set_xlim(0.08,200); b.set_ylim(len(cls)-0.3,-0.7)
b.set_yticks(range(len(cls))); b.set_yticklabels([l for _,l in cls],fontsize=7); b.tick_params(labelsize=7)
b.set_xlabel('Passes detected (%)',fontsize=7)
b.set_title('Blind search, 24,284 passes of a\nvalidation campaign: detections per class',fontsize=7.5)
for s_ in ['top','right']: b.spines[s_].set_visible(False)
e=ax[2]; SS=json.load(open('../results/shadow_strata.json'))['phase2|Starlink V2-mini']['elevation']['>=60']
for i,(g,lab,col) in enumerate([('plate_shadow','Earth\nshadow','#26324a'),('plate_sunlit','Sunlit','#e3b23c')]):
    v=SS[g]; e.bar(i,100*v['pd'],color=col,width=0.6,yerr=[[100*(v['pd']-v['ci'][0])],[100*(v['ci'][1]-v['pd'])]],capsize=2,error_kw=dict(lw=0.8))
    e.annotate(f"{100*v['pd']:.1f}%\n({v['k']}/{v['n']:,})",(i,100*v['ci'][1]),xytext=(0,3),textcoords='offset points',ha='center',fontsize=6)
e.set_xticks([0,1]); e.set_xticklabels(['Earth\nshadow','Sunlit'],fontsize=7); e.set_ylim(0,32); e.tick_params(labelsize=7)
e.set_ylabel('Passes detected (%)',fontsize=7); e.set_title('Starlink V2-mini,\nnear-zenith passes',fontsize=7.5)
for s_ in ['top','right']: e.spines[s_].set_visible(False)
c=ax[3]; F=json.load(open('../results/iono_fit_v2.json'))
for k,v in F.items():
    if len(k)>13: c.plot(v['oem_ms'],v['dt_ms'],MK[k[14:]],color=COL[k[14:]],ms=4,ls='')
Q=json.load(open('../results/phase2_iss_summary.json'))
fr=[x for x in Q['prediction'] if x['prop_age_h']<24]
c.plot([x['dt_oem_ms'] for x in fr],[x['dt_ms'] for x in fr],'o',mfc='none',mec='#444',ms=5,mew=0.8)
c.plot([-200,200],[-200,200],':',color='#888',lw=0.6); c.set_xlim(-200,200); c.set_ylim(-200,200); c.tick_params(labelsize=7)
c.set_xlabel('Offset predicted by NASA (ms)',fontsize=7); c.set_ylabel('Measured TLE offset $\\delta t$ (ms)',fontsize=7)
c.set_title('ISS timing vs NASA predictions:\n18 ms (Phase 1), 39 ms rms (Phase 2, ○)',fontsize=7.5)
for s_ in ['top','right']: c.spines[s_].set_visible(False)
h=[plt.Line2D([],[],marker=MK[s],color=COL[s],ls='',ms=4) for s in COL]
fig.legend(h,[NAME[s] for s in COL],title='Receivers (map and Phase 1 ISS timing)',title_fontsize=7,loc='lower center',ncol=5,fontsize=7,bbox_to_anchor=(0.5,-0.14),frameon=False)
plt.savefig('../results/figures/graphical_abstract.png',dpi=300,bbox_inches='tight',pad_inches=0.05)
