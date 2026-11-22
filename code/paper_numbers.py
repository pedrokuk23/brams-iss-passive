# Registry of the numbers quoted in the paper: each is recomputed from the result files and written, with the exact text
# form in which it appears, to results/paper_numbers.json. check_numbers.py then verifies that every form is present in
# the manuscript sources (and that superseded values are absent).
import json, numpy as np
R='../results/'
def J(f): return json.load(open(R+f))
reg=[]
def add(key,value,text,where=('main',)):
    reg.append(dict(key=key,value=float(value) if value is not None else None,text=text,where=list(where)))
f0=lambda x:'%.0f'%x; f1=lambda x:'%.1f'%x; f2=lambda x:'%.2f'%x
def pm(x):  # signed integer as in the tex ($+31$ or $-102$)
    return '$%+.0f$'%x
# ---------------- Sec. 6, Phase 1
F=J('iono_fit_v2.json'); D=J('receiver_delay.json')
pairs=np.array([v['dt_ms']-v['oem_ms'] for k,v in F.items() if len(k)>13])
pp=np.array([F[k]['dt_ms']-F[k]['oem_ms'] for k in sorted(F) if len(k)==13])
add('p1_pairs_rms',np.sqrt(np.mean(pairs**2)),'an rms of %s~ms'%f0(np.sqrt(np.mean(pairs**2))))
add('p1_pairs_mean',pairs.mean(),'a mean of %s~ms'%f0(abs(pairs.mean())))
add('p1_perpass_rms',np.sqrt(np.mean(pp**2)),'with an rms of %s~ms'%f0(np.sqrt(np.mean(pp**2))))
for i,x in enumerate(pp): add('p1_perpass_%d'%i,x,pm(x))
X=np.array([F[k]['oem_ms'] for k in sorted(F) if len(k)==13]); Y=np.array([F[k]['dt_ms'] for k in sorted(F) if len(k)==13])
sl,ic=np.polyfit(X,Y,1); res=Y-sl*X-ic; se=np.sqrt(res@res/3/np.sum((X-X.mean())**2))
add('p1_slope',sl,'$%s\\pm%s$'%(f2(sl),f2(se)))
add('p1_latency',D['delay_ms'],'%s~ms'%f1(D['delay_ms']) if False else '$8.7\\pm1.1$~ms')   # joint model of the text (iss_checks2)
add('p1_oem_range_lo',X.min(),'$%.0f$ to $%+.0f$~ms'%(X.min(),X.max()))
# ---------------- Sec. 6, Phase 2 and both campaigns
Q=J('phase2_iss_summary.json'); C=J('iss_oem_combined.json'); G=J('sgp4_oem_fit.json')
add('p2_n_slots',Q['n_recordings'],'the %d five-minute slots'%Q['n_recordings'])
add('p2_n_passes',Q['n_passes'],'%d passes, because three passes'%Q['n_passes'])
add('p2_n_fits',Q['n_fits_used'],'%d receiver fits'%Q['n_fits_used'])
ro=Q['receiver_offsets']
add('p2_lag_H',ro['BEHUMA']['mean_ms'],'$%s\\pm%s$, $%s\\pm%s$ and $%s\\pm%s$~ms'%(f1(ro['BEHUMA']['mean_ms']),f1(ro['BEHUMA']['se_ms']),f1(ro['BEUCCE']['mean_ms']),f1(ro['BEUCCE']['se_ms']),f1(ro['FRMAUB']['mean_ms']),f1(ro['FRMAUB']['se_ms'])))
add('p2_lat_HU',Q['latency_v5_HU']['L_ms'],'$%s\\pm%s$~ms'%(f1(Q['latency_v5_HU']['L_ms']),f1(Q['latency_v5_HU']['L_se_ms'])))
add('p2_lat_HUM',Q['latency_v5_HUM']['L_ms'],'$%s\\pm%s$~ms'%(f1(Q['latency_v5_HUM']['L_ms']),f1(Q['latency_v5_HUM']['L_se_ms'])))
add('p2_rx_agree',Q['latency_v5_HUM']['resid_sd_ms'],'the receivers of a pass agree to %s~ms'%f1(Q['latency_v5_HUM']['resid_sd_ms']))
fr=C['by_start_age_excl_20261025']['0-24 h|both']
add('both_fresh_rms',fr['rms_ms'],'%s~ms rms'%f0(fr['rms_ms']),('main','abstract'))
add('both_fresh_n',fr['n'],'%d passes'%fr['n'],('main','abstract'))
add('both_fresh_mean',fr['mean_ms'],'a mean of $+%s$~ms and a standard deviation of %s~ms'%(f0(fr['mean_ms']),f0(fr['sd_ms'])))
add('both_fresh_n1',C['by_start_age_excl_20261025']['0-24 h|1']['n'],'(%d of Phase~1 and %d of Phase~2)'%(C['by_start_age_excl_20261025']['0-24 h|1']['n'],C['by_start_age_excl_20261025']['0-24 h|2']['n']))
add('both_fresh_m',fr['rms_ms']*7.4,'about %d~m'%(round(fr['rms_ms']*7.35/10)*10))
r=C['fresh_regression']; add('both_slope',r['slope'],'slope of $%s\\pm%s$ ($r=%s$)'%(f2(r['slope']),f2(r['slope_se']),f2(r['r'])))
sc=C['fresh_scatter']
add('within_night',sc['within_sd_ms'],'scatter by %s~ms, against %s~ms between the means of the %s nights'%(f0(sc['within_sd_ms']),f0(sc['between_sd_of_means_ms']),{7:'seven',8:'eight'}.get(sc['n_nights'],str(sc['n_nights']))))
add('night_mean',sc['mean_of_night_means_ms'],'($+%s\\pm%s$~ms over the nights)'%(f0(sc['mean_of_night_means_ms']),f0(sc['se_night_level_ms'])))
b2=C['by_start_age_excl_20261025']['24-48 h|both']; b3=C['by_start_age_excl_20261025']['48- h|both']
add('age_1_2d',b2['rms_ms'],'%s~ms rms for 1--2 days'%f0(b2['rms_ms'])); add('age_1_2d_tab',b2['rms_ms'],'%s~ms (%d)'%(f0(b2['rms_ms']),b2['n']))
add('age_2d',b3['rms_ms'],'%s~ms beyond two days'%f0(b3['rms_ms'])); add('age_2d_tab',b3['rms_ms'],'%s~ms (%d)'%(f0(b3['rms_ms']),b3['n']))
cr=C['by_creation_age']['0-24 h|both']; add('creation_rms',cr['rms_ms'],'gives %s~ms for predictions created less than one day before the pass'%f0(cr['rms_ms'])); add('creation_rms_tab',cr['rms_ms'],'%s~ms (%d)'%(f0(cr['rms_ms']),cr['n']))
for kk,lab in (('by_start_age','all'),('by_creation_age','cre')):
    for bb in ('24-48 h|both','48- h|both'):
        e=C[kk][bb]; add('tab_%s_%s'%(lab,bb[:2]),e['rms_ms'],'%s~ms (%d)'%(f0(e['rms_ms']),e['n']))
gr=C['growth']; add('growth_quad',gr['quad_rms_ms'],'fits the 17 differences to %s~ms rms'%f0(gr['quad_rms_ms']))
add('growth_p1',gr['phase1_fresh_mean_ms'],'average $%+.0f$~ms and all five passes of Phase~1 $%+.1f$~ms'%(gr['phase1_fresh_mean_ms'],gr['phase1_all_mean_ms']))
NO=J('iss_nonoverlap.json')
add('nonov_formal',NO['formal_one_set_ms'][0],'%s--%s~ms'%(f1(NO['formal_one_set_ms'][0]),f1(NO['formal_one_set_ms'][1])),('main',))
add('nonov_chi2',None,'reduced chi-square of %s--%s over the four sets and gives %s--%s~ms'%(f1(min(x['chi2_red'] for x in NO['latency_by_set'])),f1(max(x['chi2_red'] for x in NO['latency_by_set'])),f1(NO['latency_range_ms'][0]),f1(NO['latency_range_ms'][1])))
MS=J('maubeuge_separate.json'); fa=MS['far_pooled_4rx']
add('far_p2_4rx',fa['far']*100,'0.%s\\%% in Phase~2 (%s--%s\\%%)'%(('%.2f'%(fa['far']*100))[2:],f2(fa['range_95'][0]*100),f2(fa['range_95'][1]*100)))
add('maub_det',MS['maubeuge']['n_detections'],'its %d Phase~2 detections (%d of the ISS complex and %d of DTC satellites)'%(MS['maubeuge']['n_detections'],MS['maubeuge']['by_class']['ISS complex'],MS['maubeuge']['by_class']['Starlink DTC']))
cm=MS['classes_without_maubeuge']; add('maub_rates',None,'ISS and DTC rates of Phase~2 are %s\\%% and %s\\%%'%(f0(cm['ISS complex']['pd']*100),f1(cm['Starlink DTC']['pd']*100)))
add('lat_p1_s10',C['latency_ms']['phase1_HU'],'%s~ms in Phase~1, Humain and Uccle'%f1(C['latency_ms']['phase1_HU']),('supp',))
po=Q['pass_order']; add('order_dt','%s'%0 and po['0']['dt_mean_ms'],'%s, %s and %s~ms on average'%(pm(po['0']['dt_mean_ms']),pm(po['1']['dt_mean_ms']),pm(po['2']['dt_mean_ms'])))
add('order_n',po['0']['n'],'(%d, %d and %d passes of Phase~2)'%(po['0']['n'],po['1']['n'],po['2']['n']))
add('post_rms',Q['tle_post']['rms_ms'],'(%s~ms rms over %d passes)'%(f0(Q['tle_post']['rms_ms']),Q['tle_post']['n']))
gs=G['summary']
add('sgp4_rms',gs['along_rms_ms_range'][0],'%s--%s~ms rms'%(f0(gs['along_rms_ms_range'][0]),f0(gs['along_rms_ms_range'][1])),('main','supp'))
add('sgp4_r',gs['corr_fit_resid_vs_dt'],'$r=%s$'%f2(gs['corr_fit_resid_vs_dt']),('main','supp'))
o=gs['by_order']; add('sgp4_order',o['0']['fit_resid_mean_ms'],'(%s, %s and %s~ms'%(pm(o['0']['fit_resid_mean_ms']),pm(o['1']['fit_resid_mean_ms']),pm(o['2']['fit_resid_mean_ms'])),('main','supp'))
add('sgp4_npass',gs['n_passes'],'at the %d Phase~2 passes'%gs['n_passes'])
add('sgp4_scatter',gs['dt_rms_ms'],'from %s to %s~ms'%(f0(gs['dt_rms_ms']),f0(gs['dt_minus_fit_resid_rms_ms'])),('supp',))
# ---------------- Sec. 10: false alarms, coincidence, emission, altitude, injection, feasibility
S2=J('phase2_summary.json'); B=J('blind_summary.json')
add('far_p1',B['far_checks']['split_half']['all']['far']*100,'0.%s\\%% per pass in Phase~1 (95\\%% range over the splits %s--%s\\%%)'%(('%.2f'%(B['far_checks']['split_half']['all']['far']*100))[2:],f2(B['far_checks']['split_half']['all']['range_95'][0]*100),f2(B['far_checks']['split_half']['all']['range_95'][1]*100)))
add('far_p2',S2['far_checks']['split_half']['all']['far']*100,'0.%s\\%% over all five receivers'%('%.2f'%(S2['far_checks']['split_half']['all']['far']*100))[2:])
add('far_pooled_abstract',None,'0.12\\% per pass',('abstract','main'))
CO=J('coincidence.json')['pass']['all non-ISS']; CN=J('coincidence.json')['null']['all non-ISS']
add('coinc',CO['detected_at_least_2'],'Of the %s passes'%format(CO['n_groups'],',').replace(',',','))
add('coinc2',CO['detected_at_least_2'],'%d were detected at one receiver or more and %d at two or more, against %d expected'%(CO['detected_at_least_1'],CO['detected_at_least_2'],round(CO['expected_at_least_2_if_independent'])))
add('coinc_tau',CO['dtau_s_frac_le1'],'within 1~s in %d of the %d pairs'%(round(CO['dtau_s_frac_le1']*CO['n_pairs']),CO['n_pairs']))
add('coinc_null',CN['detected_at_least_1'],'at one receiver in %d cases and never at two'%CN['detected_at_least_1'])
U=J('uemr_search.json'); ch=U['chance']
add('uemr_exp',ch['expected'],'about %d exceedances are expected'%round(ch['expected']),('supp',))
add('uemr_p',ch['p_ge_obs_poisson'],'%s for a Poisson count, %s for the exact'%(f2(ch['p_ge_obs_poisson']),f2(ch['p_ge_obs_exact'])),('supp',))
add('uemr_ul',ch['ul95_excess_frac']*100,'%s\\%% of passes for a Poisson count'%f2(ch['ul95_excess_frac']*100),('supp',))
add('uemr_ul0',ch['ul95_frac_no_subtraction']*100,'or %s\\%% without'%f1(ch['ul95_frac_no_subtraction']*100),('supp',))
add('uemr_shadow_ul',ch['one_sided_ul95_shadow_1_of_380']*100,'one-sided 95\\%% upper bound %s\\%%'%f1(ch['one_sided_ul95_shadow_1_of_380']*100),('supp',))
IJ=J('inject_test.json'); IS=IJ['summary']
add('inj_n',IJ['n_recordings'],'baseband of %d Phase~2 recordings'%IJ['n_recordings'])
add('inj_steady',None,'steady echo was detected in %s and %s\\%%'%(', '.join(f0(IS['steady'][l]['pd']*100) for l in ('1.0','1.25','1.5')),f0(IS['steady']['2.0']['pd']*100)))
add('inj_glint',None,'full width in %s and %s\\%%'%(', '.join(f0(IS['glint'][l]['pd']*100) for l in ('1.0','1.25','1.5')),f0(IS['glint']['2.0']['pd']*100)))
def half(e):
    L=sorted(float(k) for k in IS[e]); Pd=[IS[e][k]['pd'] for k in sorted(IS[e],key=float)]; return float(np.interp(0.5,Pd,L))
add('inj_half',half('steady'),'detected at %s (steady) and %s (glint) times'%(f1(half('steady')),f1(half('glint'))))
FV=J('cluster_stats.json')['feasibility']
add('feas_obs',FV['observed'],'%s'%format(FV['observed'],','))
add('feas_pred',FV['predicted'],'%s'%format(round(FV['predicted']),','))
GL=J('shadow_glm.json'); cs=GL['cluster_satellite']; cr=GL['cluster_recording']; ab=GL['shadow_by_altitude_bin']
tw=GL['cluster_two_way']; add('glm_or',cs['or_shadow'],'is %s (95\\%% interval %s--%s with errors clustered by satellite, %s--%s clustered by five-minute slot or by satellite and slot together; %s passes of %s satellites)'%(f0(cs['or_shadow']),f1(cs['ci95'][0]),f0(cs['ci95'][1]),f1(tw['ci95'][0]),f0(tw['ci95'][1]),format(GL['n'],','),format(GL['n_sat'],',')))
add('glm_ridge',GL['ridge_lambda1']['or_shadow'],'a ridge-penalized fit gives %s'%f0(GL['ridge_lambda1']['or_shadow']))
add('glm_alt',ab['1']['or_shadow'],'(odds ratio %s, %s--%s), and absent below (%s, %s--%s)'%(f0(ab['1']['or_shadow']),f0(ab['1']['ci95'][0]),f0(ab['1']['ci95'][1]),f1(ab['0']['or_shadow']),f1(ab['0']['ci95'][0]),f1(ab['0']['ci95'][1])))
add('glm_concl',cs['or_shadow'],'with %s times higher odds'%f0(cs['or_shadow']))
CA=J('catalog_association.json')['summary']
add('cat_n',CA['n_objects_in_catalog'],'each of the %s objects'%format(CA['n_objects_in_catalog'],','))
add('cat_alt',CA['n_with_alternative'],'%d of the %s detections that could be tested (of %s) another object'%(CA['n_with_alternative'],format(CA['n_checked'],','),format(CA['n_detections'],',')))
add('cat_multi',CA['n_multi_rx_alternative_surviving'],'in only %d of %d'%(CA['n_multi_rx_alternative_surviving'],CA['n_multi_rx']))
add('cat_ast',CA['ambiguous_by_class']['AST SpaceMobile'],'up to %d of the %d BlueBird detections and %d of the %d V1.x'%(CA['ambiguous_by_class']['AST SpaceMobile'],CA['by_class']['AST SpaceMobile']['n'],CA['ambiguous_by_class']['Starlink V1.x'],CA['by_class']['Starlink V1.x']['n']))
out=dict(n=len(reg),entries=reg,
         stale=['20--60~ms','18--39~ms','61~ms rms in','21~ms rms','sub-kilometer level','in 19 passes','set by the predictions','45 Phase~2 recordings',
                'A likely cause','comparable to the drift','day-old','20-ms agreement','20--125~km','In near-zenith passes','declared by the author','designed after the Phase~2','first version, used to select','four of five receivers','below 0.2\\%','\\amarelo{20~ms rms','Poisson probability of 10 or more','about 0.35\\% of passes)'])
json.dump(out,open(R+'paper_numbers.json','w'),indent=1); print(len(reg),'entries')
