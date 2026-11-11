# Table S2 of the Supplementary material (Section S10): every ISS pass of both campaigns against the latest NASA OEM file
# created before it, with the two ages of that prediction and the residual of the SGP4 fit to the OEM (sgp4_oem_fit.py).
# Reads results/iss_oem_combined.json and results/sgp4_oem_fit.json; writes results/table_s10.tex.
import json
C=json.load(open('../results/iss_oem_combined.json')); G=json.load(open('../results/sgp4_oem_fit.json'))
fit={}
for lab in sorted(G['files']):
    for p in G['files'][lab]['passes']: fit[p['pass_']]=p['fit_resid_ms']
ab={'BEOTTI':'O','BEOUDS':'D','BEHUMA':'H','BEUCCE':'U','FRMAUB':'M'}
L=[r'\begin{tabular}{llcrlrrrrr}',r'\toprule',
   r'Pass (UTC) & Ph. & Receivers & $\delta t$ & OEM file & Age$_\mathrm{start}$ & Age$_\mathrm{created}$ & $\delta t_\mathrm{OEM}$ & Difference & SGP4 fit \\',
   r' & & & (ms) & & (h) & (h) & (ms) & (ms) & (ms) \\',r'\midrule']
for r in sorted(C['rows'],key=lambda r:r['pass_']):
    p=r['pass_']; t='%s %s %s:%s'%(p[6:8],{'10':'Oct','11':'Nov'}[p[4:6]],p[9:11],p[11:13])
    f=r['oem'].replace('T',' ')[:13]; f='%s %s'%(f[6:8]+' '+{'10':'Oct','11':'Nov'}[f[4:6]],f[9:13])
    sg=('%+.0f'%fit[p]) if p in fit and r['campaign']==2 else '--'
    L.append(r'%s & %d & %s & $%+.0f$ & %s & %.0f & %.0f & $%+.0f$ & $%+.0f$ & %s \\'%(t,r['campaign'],''.join(ab[x] for x in sorted(r['receivers'],key=list(ab).index)),
             r['dt_ms'],f,r['age_start_h'],r['age_created_h'],r['dt_oem_ms'],r['diff_ms'],('$'+sg+'$') if sg!='--' else sg))
L+=[r'\bottomrule',r'\end{tabular}']
open('../results/table_s10.tex','w').write('\n'.join(L)+'\n'); print(len(C['rows']),'rows')
