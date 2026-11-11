# Remoção de linhas espectrais estreitas (interferências CW) do sinal em banda base dos receptores SDR
import numpy as np
from scipy import signal as sg
def lines(bb,fs,thr=10,fmin=5):
    """Linhas estreitas: picos do espectro (resolução de 0,37 Hz) mais de `thr` dB acima do patamar local
    (mediana móvel de ~37 Hz), fora de |f|<fmin (resíduo do sinal direto)."""
    f,P=sg.welch(bb,fs,nperseg=16384,return_onesided=False); o=np.argsort(f); f=f[o]; L=10*np.log10(P[o]+1e-30)
    floor=sg.medfilt(L,101); E=L-floor
    pk,_=sg.find_peaks(E,height=thr,distance=8); pk=[i for i in pk if abs(f[i])>fmin]
    return [float(f[i]) for i in pk],[float(E[i]) for i in pk]
def clean(bb,ts,fs,fl=None):
    if fl is None: fl,_=lines(bb,fs)
    # equivalente, no domínio da frequência, a subtrair a média móvel de 1 s em torno de cada linha:
    # resposta 1 - sinc((f - f_l) * 1 s) para cada linha f_l
    n=len(bb); f=np.fft.fftfreq(n,1/fs); H=np.ones(n)
    for f0 in fl:
        k=np.abs(f-f0)<20                     # o lóbulo relevante; fora disso o efeito é < 2%
        H[k]*=1-np.sinc((f[k]-f0)*1.0)
    return np.fft.ifft(np.fft.fft(bb)*H)
def maybe_clean(bb,ts,fs):
    """Aplica só aos receptores SDR (fs=6048); o receptor analógico de Ottignies não é alterado."""
    return clean(bb,ts,fs) if fs>6000 else bb
