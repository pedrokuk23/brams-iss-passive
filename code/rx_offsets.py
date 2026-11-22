# Receiver timing offsets from the ISS passes (Sec. 6). With one transmitter position for every receiver (bwav.py),
# delta t of a pass differs between receivers by fixed amounts: the version-5 receivers (Humain, Uccle, Maubeuge) lag
# and Ottignies leads Oudsbergen. Model: delta t[pass, receiver] = a[pass] + b[group of the receiver], with Oudsbergen as
# the reference (b = 0), fitted by least squares; the formal errors of the individual fits give the reduced chi-square.
import numpy as np
V5=('BEHUMA','BEUCCE','FRMAUB')
def group(rx, v5=V5):
    if rx=='BEOTTI': return 'ottignies'
    if rx in v5: return 'v5'
    return None
def fit_offsets(rows, v5=V5):
    """rows: (pass, receiver, dt_ms, sdt_ms or None). Returns offsets (ms) with standard errors, residual scatter."""
    rows=[r for r in rows]; P=sorted({r[0] for r in rows}); G=sorted({group(r[1],v5) for r in rows}-{None})
    A=np.zeros((len(rows),len(P)+len(G))); y=np.array([r[2] for r in rows],float)
    for i,r in enumerate(rows):
        A[i,P.index(r[0])]=1; g=group(r[1],v5)
        if g: A[i,len(P)+G.index(g)]=1
    x,*_=np.linalg.lstsq(A,y,rcond=None); e=y-A@x; dof=len(y)-np.linalg.matrix_rank(A); s2=e@e/dof
    C=np.linalg.pinv(A.T@A)*s2
    s=np.array([r[3] if r[3] else np.nan for r in rows],float)
    out=dict(n=len(y),n_pass=len(P),dof=int(dof),resid_sd_ms=float(np.sqrt(s2)),resid_rms_ms=float(np.sqrt(np.mean(e**2))),
             chi2_red=float(np.nansum((e/s)**2)/dof) if np.isfinite(s).any() else None,offsets_ms={},offsets_se_ms={})
    for j,g in enumerate(G):
        out['offsets_ms'][g]=float(x[len(P)+j]); out['offsets_se_ms'][g]=float(np.sqrt(C[len(P)+j,len(P)+j]))
    out['residuals']=[dict(pass_=r[0],rx=r[1],res_ms=float(v)) for r,v in zip(rows,e)]
    return out
def correction(rx, off, v5=V5):
    g=group(rx,v5); return off['offsets_ms'].get(g,0.0) if g else 0.0
