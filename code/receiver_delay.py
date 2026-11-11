# Atraso fixo dos receptores que gravam a versão 5 do formato (Humain, Uccle): diferença de δt em relação à média
# de Ottignies e Oudsbergen (versão 4) em cada passagem, correção comum e teste deixando uma passagem de fora.
import json, numpy as np
F=json.load(open('../results/iono_fit_v2.json'))
P=['20261103_0040','20261103_0220','20261103_2355','20261104_0130','20261104_2305']
rec={k:v for k,v in F.items() if len(k)>13}
ref={p:np.mean([F[p+'_'+s]['dt_ms'] for s in ['BEOTTI','BEOUDS']]) for p in P}
diff={k:rec[k]['dt_ms']-ref[k[:13]] for k in rec if k[14:] in ('BEHUMA','BEUCCE')}
v=np.array(list(diff.values())); D=v.mean()
loo={k:diff[k]-np.mean([diff[j] for j in diff if j!=k]) for k in diff}
dt={k:rec[k]['dt_ms']-(D if k in diff else 0) for k in rec}
dev=[dt[k]-np.mean([dt[j] for j in dt if j[:13]==k[:13]]) for k in dt]
oem=[dt[k]-rec[k]['oem_ms'] for k in dt]
out=dict(differences_ms=diff,delay_ms=float(D),sem_ms=float(v.std(ddof=1)/np.sqrt(len(v))),loo_residual_ms=loo,
         rms_dev_after_ms=float(np.sqrt(np.mean(np.square(dev)))),rms_oem_after_ms=float(np.sqrt(np.mean(np.square(oem)))))
json.dump(out,open('../results/receiver_delay.json','w'),indent=1); print(json.dumps(out,indent=1))
