"""Basic automated checks of the core formulas and of the published numbers.
Run from the repository root:  python -m pytest -q tests
Tests that need the raw BRAMS WAV files are skipped when data/raw_wav/ is empty."""
import json, os, sys, glob
import numpy as np
import pytest
from scipy import signal

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(ROOT, 'code'); RES = os.path.join(ROOT, 'results')
sys.path.insert(0, CODE)
from linefilt import lines, clean            # noqa: E402

def ecef(lat, lon, h):                        # same expression as code/station_run.py
    a = 6378.137; f = 1/298.257223563; e2 = f*(2-f); la, lo = np.radians(lat), np.radians(lon)
    N = a/np.sqrt(1-e2*np.sin(la)**2)
    return np.stack([(N+h)*np.cos(la)*np.cos(lo), (N+h)*np.cos(la)*np.sin(lo), (N*(1-e2)+h)*np.sin(la)], -1)

def test_station_code_matches():
    src = open(os.path.join(CODE, 'station_run.py')).read()
    assert 'f=1/298.257223563' in src and 'lam=299792.458/49.97e6' in src

def test_ecef_reference_points():
    assert np.allclose(ecef(0, 0, 0), [6378.137, 0, 0], atol=1e-6)
    assert np.allclose(ecef(90, 0, 0), [0, 0, 6356.752314], atol=1e-3)

def test_bistatic_doppler_sign_and_value():
    lam = 299792.458/49.97e6                  # km
    TX = np.array([0., 0., 0.]); RX = np.array([10., 0., 0.])
    t = np.linspace(0, 10, 1001); v = 0.2     # km/s, receding vertically above the midpoint
    X = np.stack([5+0*t, 0*t, 100+v*t], 1)
    R = np.linalg.norm(X-TX, axis=1)+np.linalg.norm(X-RX, axis=1)
    fd = -np.gradient(R, t)/lam
    expected = -2*v*(100/np.hypot(100, 5))/lam
    assert fd[500] < 0 and abs(fd[500]-expected)/abs(expected) < 1e-3

def test_direct_signal_suppression_response():
    fs = 6048.; t = np.arange(int(20*fs))/fs; w = int(fs)
    for f, maxloss_db in [(3.0, 1.0), (30.0, 0.1)]:
        x = np.exp(2j*np.pi*f*t) + 5.0         # echo plus a strong constant (direct signal at 0 Hz)
        y = x-signal.oaconvolve(x, np.ones(w)/w, mode='same')
        k = slice(2*w, -2*w)
        assert abs(np.mean(y[k])) < 1e-2       # direct signal removed
        loss = -20*np.log10(np.abs(np.vdot(np.exp(2j*np.pi*f*t[k]), y[k]))/len(t[k]))
        assert loss < maxloss_db

def test_line_removal_keeps_swept_echo():
    rng = np.random.default_rng(1); fs = 6048.; n = int(300*fs); t = np.arange(n)/fs
    noise = (rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
    line = 30*np.exp(2j*np.pi*500.0*t)
    f0, rate = 800., -10.                     # echo sweeping 800 -> -2200 Hz
    chirp = 0.5*np.exp(2j*np.pi*(f0*t+0.5*rate*t**2))
    fl, lv = lines(noise+line+chirp, fs)
    assert any(abs(f-500) < 0.5 for f in fl)
    out = clean(noise+line+chirp, t, fs, fl)
    resid = np.abs(np.vdot(np.exp(2j*np.pi*500*t), out)/n)
    assert resid < 30*10**(-20/20)            # line attenuated by more than 20 dB
    kept = np.abs(np.vdot(chirp, out))/np.vdot(chirp, chirp).real
    assert kept > 0.95                        # swept echo preserved

def test_published_numbers():
    F = json.load(open(os.path.join(RES, 'iono_fit_v2.json')))
    rec = {k: v for k, v in F.items() if len(k) > 13}
    d = np.array([v['dt_ms']-v['oem_ms'] for v in rec.values()])
    assert len(d) == 16 and abs(np.sqrt(np.mean(d**2))-20.2) < 0.2   # 10-point Lagrange OEM interpolation
    passes = [k for k in F if len(k) == 13]
    pp = np.array([F[p]['dt_ms']-F[p]['oem_ms'] for p in passes])
    assert abs(np.sqrt(np.mean(pp**2))-18.1) < 0.2
    D = json.load(open(os.path.join(RES, 'receiver_delay.json')))
    assert abs(D['delay_ms']-8.8) < 0.1 and abs(D['rms_dev_after_ms']-1.6) < 0.05

@pytest.mark.skipif(not glob.glob(os.path.join(ROOT, 'data', 'raw_wav', '*.wav')), reason='raw BRAMS WAV files not present')
def test_raw_checksums():
    import hashlib
    bad = []
    for l in open(os.path.join(ROOT, 'data', 'raw_wav', 'SHA256SUMS')):
        h, n = l.split()
        p = os.path.join(ROOT, 'data', 'raw_wav', n)
        if os.path.exists(p) and hashlib.sha256(open(p, 'rb').read()).hexdigest() != h:
            bad.append(n)
    assert not bad


# ---------------------------------------------------------------------------------------------------------------
# Independent numerical checks

def _pipeline_namespace():
    """Executes the geometry part of code/station_run.py (ECEF, ISS propagation) and returns its namespace."""
    cwd = os.getcwd(); os.chdir(CODE); argv = sys.argv
    try:
        sys.argv = ['x', 'X']; ns = {}
        exec(open('station_run.py').read().split('def load(')[0], ns)
    finally:
        os.chdir(cwd); sys.argv = argv
    return ns

def test_teme_to_itrs_against_astropy():
    """TEME->Earth-fixed rotation of the pipeline (GMST IAU-82, no polar motion) against astropy's TEME->ITRS
    (IAU-82 GMST plus polar motion). The difference must be at the polar-motion level (tens of metres at most)."""
    astropy = pytest.importorskip('astropy')
    from astropy.utils import iers; iers.conf.auto_download = False; iers.conf.auto_max_age = None   # offline: bundled tables only
    from astropy.coordinates import TEME, ITRS, CartesianRepresentation
    from astropy.time import Time; import astropy.units as u
    import datetime as dt
    from sgp4.api import jday
    ns = _pipeline_namespace()
    for hh, mm in [(0, 42), (1, 32), (23, 7)]:
        d = dt.datetime(2026, 11, 4, hh, mm, tzinfo=dt.timezone.utc)
        X = ns['iss_ecef']([d.timestamp()])[0]
        jd, fr = jday(d.year, d.month, d.day, d.hour, d.minute, d.second)
        e, r, v = ns['sat'].sgp4(jd, fr); t = Time(jd, fr, format='jd', scale='utc')
        Y = TEME(CartesianRepresentation(r*u.km), obstime=t).transform_to(ITRS(obstime=t)).cartesian.xyz.to(u.km).value
        assert np.linalg.norm(X-Y)*1e3 < 30          # m; 30 m along track is about 4 ms of timing

def test_ecef_against_astropy():
    astropy = pytest.importorskip('astropy')
    from astropy.coordinates import EarthLocation; import astropy.units as u
    ns = _pipeline_namespace()
    for lat, lon, h in [(50.097, 4.588, 0.2), (50.795, 4.581, 0.12), (51.0, 5.5, 0.05), (-23.2, -45.9, 0.6)]:
        L = EarthLocation.from_geodetic(lon, lat, h*1e3)
        assert np.allclose(ns['ecef'](lat, lon, h), [L.x.to(u.km).value, L.y.to(u.km).value, L.z.to(u.km).value], atol=1e-6)

def test_bistatic_doppler_numeric_vs_analytic_velocity():
    """The pipeline differentiates the bistatic range numerically. Compare with the analytic expression
    f_D = -(1/lambda) (u_T + u_R) . V, using the SGP4 velocity rotated to the Earth-fixed frame (with the w x r term)."""
    import datetime as dt
    from sgp4.api import jday
    from sgp4.propagation import gstime
    ns = _pipeline_namespace(); lam = ns['lam']
    TX = ns['ecef'](50.097, 4.588, 0.2); RX = ns['ecef'](50.795, 4.581, 0.12)
    t0 = dt.datetime(2026, 11, 4, 1, 32, tzinfo=dt.timezone.utc).timestamp()
    t = t0+np.arange(-2, 2.01, 0.5); X = ns['iss_ecef'](t)
    Rb = np.linalg.norm(X-TX, axis=1)+np.linalg.norm(X-RX, axis=1); fd_num = -np.gradient(Rb, t)/lam
    d = dt.datetime.fromtimestamp(t0, dt.timezone.utc); jd, fr = jday(d.year, d.month, d.day, d.hour, d.minute, d.second)
    e, r, v = ns['sat'].sgp4(jd, fr); g = gstime(jd+fr); w = 7.292115e-5
    Rz = np.array([[np.cos(g), np.sin(g), 0], [-np.sin(g), np.cos(g), 0], [0, 0, 1]])
    x = Rz@np.array(r); V = Rz@np.array(v)-np.cross([0, 0, w], x)
    fd_an = -((x-TX)/np.linalg.norm(x-TX)+(x-RX)/np.linalg.norm(x-RX))@V/lam
    assert abs(fd_num[4]-fd_an) < 0.05               # Hz, out of a Doppler of several hundred Hz

def test_timing_error_recovered_from_residual():
    """Residual model: an object that is late by dt produces Delta f = f_D(t-dt) - f_D(t) = -fdot*dt.
    A 40 ms delay injected into a real ISS geometry must be recovered by the one-parameter least-squares fit."""
    ns = _pipeline_namespace(); lam = ns['lam']
    import datetime as dt
    TX = ns['ecef'](50.097, 4.588, 0.2); RX = ns['ecef'](50.795, 4.581, 0.12)
    t0 = dt.datetime(2026, 11, 4, 1, 31, tzinfo=dt.timezone.utc).timestamp()
    def fD(t):
        tt = np.arange(t.min()-2, t.max()+2, 0.1); X = ns['iss_ecef'](tt)
        Rb = np.linalg.norm(X-TX, axis=1)+np.linalg.norm(X-RX, axis=1); return np.interp(t, tt, -np.gradient(Rb, tt)/lam), tt, Rb
    tc = t0+np.arange(0, 120, 2.0); f, tt, Rb = fD(tc); fl, _, _ = fD(tc-0.040)
    fdot = np.interp(tc, tt, np.gradient(-np.gradient(Rb, tt)/lam, tt))
    dt_hat = np.linalg.lstsq((-fdot)[:, None], fl-f, rcond=None)[0][0]
    assert abs(dt_hat-0.040) < 0.002

def test_chapman_vertical_content():
    """Vertical integral of the alpha-Chapman layer used in the ionospheric term: Nm*H*sqrt(2*pi*e)."""
    def chap(h, hm=300., H=55.):
        z = (h-hm)/H; return np.exp(0.5*(1-z-np.exp(-z)))
    h = np.linspace(-200, 5000, 200001)
    assert abs(np.trapezoid(chap(h), h)/(55*np.sqrt(2*np.pi*np.e))-1) < 1e-3

def test_file_manifests_consistent():
    names = [l.strip() for l in open(os.path.join(ROOT, 'data', 'raw_wav', 'FILES.txt')) if l.strip()]
    sums = [l.split()[1] for l in open(os.path.join(ROOT, 'data', 'raw_wav', 'SHA256SUMS')) if l.strip()]
    assert len(names) == 293 and len(set(names)) == 293 and sorted(names) == sorted(sums)

@pytest.mark.skipif(not glob.glob(os.path.join(ROOT, 'data', 'raw_wav', '*.wav')), reason='raw BRAMS WAV files not present')
def test_raw_files_complete():
    """When raw files are present, all 293 expected files must be there (a partial set fails)."""
    names = [l.strip() for l in open(os.path.join(ROOT, 'data', 'raw_wav', 'FILES.txt')) if l.strip()]
    missing = [n for n in names if not os.path.exists(os.path.join(ROOT, 'data', 'raw_wav', n))]
    assert not missing, f'{len(missing)} expected WAV files missing, e.g. {missing[:3]}'

def test_envisat_manifest():
    names = [l.strip() for l in open(os.path.join(ROOT, 'data', 'raw_wav_envisat', 'FILES.txt')) if l.strip()]
    sums = [l.split()[1] for l in open(os.path.join(ROOT, 'data', 'raw_wav_envisat', 'SHA256SUMS')) if l.strip()]
    assert len(names) == 40 and sorted(names) == sorted(sums)

def test_statistic_energy_relation():
    """Section 9.3: z ~ (E/N0)/sqrt(nb) holds for the measured echoes within a factor of about 1.5 (median)."""
    import json
    d = json.load(open(os.path.join(ROOT, 'results', 'detect_theory.json')))
    assert 0.8 < d['z_obs_over_pred']['median'] < 1.25
