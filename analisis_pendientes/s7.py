"""s7.py -- (c) Analisis IMF por IMF (CEEMDAN, ajustes del articulo) para v1 sola y para la senal de D; conteo de extremos."""
import sys, os, json, numpy as np
R = os.environ.get('M3_CODE', '.')  # raiz del paquete de codigo del repositorio
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
import mc_common as M
from PyEMD import CEEMDAN, EMD
FS = M.FS; a1, f1, p1, a2, f2, p2, a3, f3, p3 = M.CONFIGS['v1']
print('v1:', M.CONFIGS['v1'])
def comps_v1(t):
    return {'fundamental': a1*np.sin(2*np.pi*f1*t+p1), f'{f2:.2f} Hz tone': a2*np.sin(2*np.pi*f2*t+p2), f'{f3:.0f} Hz tone': a3*np.sin(2*np.pi*f3*t+p3)}
def share(imf, c):  # fraccion de la energia de la componente c contenida en la IMF (proyeccion)
    return float(np.dot(imf, c)/np.dot(c, c))
def fmed(imf):
    from scipy.signal import hilbert
    ph = np.unwrap(np.angle(hilbert(imf))); f = np.diff(ph)*FS/2/np.pi; n = len(f); return float(np.median(f[n//10:-n//10]))
OUT = {}
for name, T in [('v1', 0.5), ('D', 0.6)]:
    t = M.tvec(T); C = comps_v1(t)
    if name == 'D':
        t_, y, env = M.transient_signal(); C['120 Hz transient'] = env*np.sin(2*np.pi*120*(t-0.2))
    else:
        y = M.base(t, 'v1')
    c = CEEMDAN(trials=30, epsilon=0.05, parallel=False); c.noise_seed(0); imfs = c.ceemdan(y, t)
    e = EMD(spline_kind='cubic'); e.MAX_ITERATION = 15; e.emd(y, t, max_imf=10); im2, res2 = e.get_imfs_and_residue()
    rows = []
    for i, imf in enumerate(imfs):
        rows.append(dict(imf=i+1, fmed=fmed(imf), rms=float(np.sqrt(np.mean(imf**2))), **{k: share(imf, v) for k, v in C.items()}))
    OUT[name] = dict(n_imfs=len(imfs), rows=rows, emd_n=len(im2), emd_fmed=[fmed(q) for q in im2])
    print(f"\n{name}: CEEMDAN {len(imfs)} IMFs (incl. residue as last?)  | EMD {len(im2)} IMFs + residue, median IF {[round(fmed(q),1) for q in im2]}")
    for r in rows: print('  IMF%d fmed=%7.1f rms=%7.3f ' % (r['imf'], r['fmed'], r['rms']) + ' '.join(f"{k}:{100*r[k]:6.1f}%" for k in C))
# extremos (sin ruido, rejilla fina) para D y para v1
fs2 = 1e6; t = np.arange(int(0.6*fs2))/fs2
yb = a1*np.sin(2*np.pi*f1*t+p1) + a2*np.sin(2*np.pi*f2*t+p2) + a3*np.sin(2*np.pi*f3*t+p3)
env = np.where(t >= 0.2, 60*np.exp(-(t-0.2)/0.015), 0); yd = yb + env*np.sin(2*np.pi*120*(t-0.2))
yf = a1*np.sin(2*np.pi*f1*t+p1)
def ext(y, lo, hi):
    d = np.diff(y); s = np.sign(d); i = np.where(s[:-1] != s[1:])[0]+1; tt = t[i]; return tt[(tt >= lo) & (tt < hi)]
for lab, y in [('D (fund.+tonos+transitorio)', yd), ('v1 (sin transitorio)', yb), ('fundamental sola', yf)]:
    e_ = ext(y, 0.2, 0.23); print(lab, 'extremos en 200-230 ms:', len(e_), [round((x-0.2)*1e3, 2) for x in e_])
e_ = ext(yb, 0.0, 0.5); print('v1 extremos en 0.5 s:', len(e_))
print('fase de la fundamental en t0=200 ms [deg]:', round(np.degrees((2*np.pi*f1*0.2+p1) % (2*np.pi)), 1), ' fase inicial del transitorio: 0 (seno desde t0)')
json.dump(OUT, open('s7.json', 'w'), indent=1)
