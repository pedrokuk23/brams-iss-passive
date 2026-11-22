"""Checks of the threshold rule of the blind search, the association
invariant of the published detections, the class totals used by code/key_numbers.py, and the glint-energy and
sensitivity formulas of the detectability model. Run from the repository root:  python -m pytest -q tests"""
import json, os, sys
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(ROOT, 'code'); RES = os.path.join(ROOT, 'results')
sys.path.insert(0, CODE)
from far_checks import rule, far_checks       # noqa: E402

def load(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p): pytest.skip(name + ' not present')
    return json.load(open(p))

def test_threshold_rule_allows_floor_far_n_exceedances():
    z = np.arange(1, 2001, dtype=float)            # 2000 nulls, rule 1e-3 -> 2 exceedances allowed
    t = rule(z, len(z), 1e-3)
    assert np.sum(z >= t) == 2
    assert rule(z[:999], 999, 1e-3) > z[:999].max()  # fewer than 1000 nulls: no exceedance allowed

def test_split_half_far_close_to_nominal_on_iid_nulls():
    rng = np.random.default_rng(0)
    N = [dict(file=f'{i:05d}_X', z=float(v)) for i in range(400) for v in rng.exponential(size=50)]
    out = far_checks(N, {}, {'X': rule([r['z'] for r in N], len(N), 1e-3)}, ['X'], lambda r: True, nsplit=50)
    assert 0.3e-3 < out['split_half']['all']['far'] < 3e-3

@pytest.mark.parametrize('summ', ['phase2_summary.json'])
def test_published_thresholds_follow_the_rule(summ):
    S = load(summ)
    for st, v in S['thresholds'].items():
        zs = [z for s, f, z, ins, ex in S['nulls'] if s == st and ins and not ex]
        assert abs(rule(zs, v['n_null'], 1e-3) - v['z_th']) < 1e-3     # the compact null list keeps 4 decimals

def test_published_detections_are_not_double_counted():
    S = load('phase2_summary.json'); by = {}
    for d in S['detections']: by.setdefault(d['file'], []).append(d)
    bad = 0
    for f, dd in by.items():
        for i in range(len(dd)):
            for j in range(i+1, len(dd)):
                a, b = dd[i], dd[j]
                bad += abs(a['t_pk']-b['t_pk']) <= 0.25 and abs(a['f_raw']-b['f_raw']) < 0.5
    assert bad == 0                                 # two detections never share the same peak

def test_class_totals_are_the_all_entry():
    S = load('phase2_summary.json')
    for c, v in S['classes'].items():
        per_rx = [x for k, x in v.items() if k != 'all']
        assert sum(x['n'] for x in per_rx) == v['all']['n'] and sum(x['k'] for x in per_rx) == v['all']['k']
    assert S['classes']['Starlink V2-mini']['all']['k'] == 1028

def test_glint_energy_tends_to_half_width():
    import brams_feasibility as BF
    assert BF.glint_energy_s(50.0, 2.0) == pytest.approx(1.0, rel=0.02)   # full width 2 s: energy-equivalent T_g/2
    assert 9.0 < BF.glint_energy_s(50.0, 20.0) < 10.0             # 20 s: slightly less, the arc truncates the side lobes
    assert BF.glint_energy_s(50.0, 1e5) == pytest.approx(50.0, rel=1e-3)   # very long glint: the arc

def test_sigma_min_scaling():
    import brams_feasibility as BF
    s50 = BF.sigma_min(400, 50.0); s25 = BF.sigma_min(400, 25.0, T_arc=50.0)
    assert 70 < s50 < 85                                          # pocket value of Eq. (11), about 77 m^2
    assert s25 / s50 == pytest.approx(2.0, rel=1e-6)              # half the energy-equivalent duration
    assert BF.sigma_min(800, 50.0) / s50 == pytest.approx(16.0, rel=1e-6)

@pytest.mark.parametrize('script,out', [('phase2_analyze.py', 'phase2_summary.json'),
                                        ('feasibility_validate.py', 'feasibility_validation.json')])
def test_analysis_scripts_run_and_reproduce(script, out, tmp_path):
    """Smoke test: the analysis scripts run without the recordings and rewrite their published result unchanged."""
    import subprocess, hashlib
    import glob, shutil
    if not os.path.exists(os.path.join(RES, 'phase2_meta.json')): pytest.skip('phase2_meta.json not present')
    if not glob.glob(os.path.join(CODE, '..', 'data', 'tle', 'phase2', '*.txt')):
        pytest.skip('Phase 2 element sets (data/tle/phase2, from auxiliary_data_phase2.zip) not present')
    pub = os.path.join(RES, out); keep = str(tmp_path / out); shutil.copy2(pub, keep)   # restored whatever happens
    before = hashlib.sha256(open(pub, 'rb').read()).hexdigest()
    try:
        r = subprocess.run([sys.executable, script], cwd=CODE, capture_output=True, text=True, timeout=1800)
        assert r.returncode == 0, r.stderr[-2000:]
        assert hashlib.sha256(open(pub, 'rb').read()).hexdigest() == before
    finally:
        shutil.copy2(keep, pub)
