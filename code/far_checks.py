# Checks of the false-alarm rule of the blind search, run at the end of blind_analyze.py and phase2_analyze.py.
# (1) The rate measured on the nulls that set the threshold equals 1e-3 by construction, so it is also measured
#     out of sample: the recordings of each receiver are split at random into two halves, the threshold is set on the
#     nulls of one half with the same rule and the exceedances are counted on the nulls of the other half (500 splits).
# (2) Effective rate with the nulls attributed to real echoes removed from the denominator as well.
# (3) Sensitivity: thresholds set on all nulls, without removing those attributed to real echoes.
import numpy as np
from scipy.stats import beta

def cp(k, n, a=0.05):
    return (0.0 if k == 0 else float(beta.ppf(a/2, k, n-k+1)), 1.0 if k == n else float(beta.ppf(1-a/2, k+1, n-k)))

def rule(zs, n, far):
    zs = np.sort(zs)[::-1]; kmax = int(np.floor(far*n))
    return float(zs[kmax]+1e-9) if len(zs) > kmax else 0.0

def far_checks(N, expl, th, STS, inside, far=1e-3, nsplit=500, seed=7):
    rng = np.random.default_rng(seed); out = {'split_half': {}, 'effective': {}}
    tot_k = tot_n = 0; per_split = np.zeros(nsplit); per_split_n = np.zeros(nsplit)
    for st in STS:
        nn = [r for r in N if r['file'].endswith(st)]
        if not nn: continue
        files = np.array(sorted({r['file'] for r in nn})); ks = ns = 0
        for s in range(nsplit):
            A = set(rng.choice(files, len(files)//2, replace=False))
            a = [r for r in nn if r['file'] in A]; b = [r for r in nn if r['file'] not in A]
            t = rule([r['z'] for r in a if inside(r) and not expl.get(id(r))], len(a), far)
            k = sum(1 for r in b if inside(r) and not expl.get(id(r)) and r['z'] >= t)
            ks += k; ns += len(b); per_split[s] += k; per_split_n[s] += len(b)
        out['split_half'][st] = dict(far=ks/ns, n_per_split=ns/nsplit)
        tot_k += ks; tot_n += ns
        ne = sum(1 for r in nn if expl.get(id(r))); n_eff = len(nn)-ne
        k = sum(1 for r in nn if inside(r) and not expl.get(id(r)) and r['z'] >= th[st])
        out['effective'][st] = dict(n_null=len(nn), n_explained=ne, n_eff=n_eff, exceed=k, far=k/n_eff, ci=cp(k, n_eff),
                                    ci_upper_zero_exceed=float(beta.ppf(0.95, 1, n_eff)) if k == 0 else None)
    r = per_split/per_split_n
    out['split_half']['all'] = dict(far=tot_k/tot_n, range_95=[float(x) for x in np.percentile(r, [2.5, 97.5])], nsplit=nsplit)
    return out
