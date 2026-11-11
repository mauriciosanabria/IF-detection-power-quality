"""d_semilla_sin_emd.py -- Escenario D: ablacion del ajuste LM guiado por EMD (11-nov-2026).

El estimador 'EMD-guided LM' de mc4_transient.py usa la EMD (CEEMDAN) solo para dar la semilla de
frecuencia del ajuste Levenberg-Marquardt sobre la senal original; el instante de inicio t0 se toma
como conocido. Aqui se repite el mismo ajuste (mismo modelo, mismas cotas, mismas 6 fases iniciales)
con la semilla obtenida SIN EMD: el pico espectral, en 60-500 Hz, de la senal original en la ventana
t0..t0+50 ms despues de restarle el ajuste de 50 Hz hecho en t > t0+150 ms.

Mismas realizaciones que mc10r_resumable.py (ruido 1032000+r) con fase cero, y que
d_fase_aleatoria.py con fase sorteada (1035000+r). No usa CEEMDAN.

Uso: python d_semilla_sin_emd.py     Salida: d_semilla_sin_emd.json junto a este archivo."""
import sys, os, json
sys.argv = sys.argv[:1]
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
import numpy as np, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import least_squares
from mc_common import FS, tvec, base, add_noise, summarize
import mc4_transient as m4
T0, TAU, FE, A0, N = m4.T0, m4.TAU, m4.FE, 60.0, 100


def lm_sin_emd(t, y):
    """Copia de la segunda mitad de mc4_transient.emd_guided, con x = y (sin IMF)."""
    late, early = t > T0+0.15, (t >= T0) & (t < T0+0.05)
    M = np.column_stack([np.sin(2*np.pi*50*t), np.cos(2*np.pi*50*t)])
    ca = np.linalg.lstsq(M[late], y[late], rcond=None)[0]; res = y - M@ca
    seg = res[early]*np.hanning(early.sum()); F = np.abs(np.fft.rfft(seg, 8192)); fr = np.fft.rfftfreq(8192, 1/FS)
    b = (fr > 60) & (fr < 500); fseed = fr[b][np.argmax(F[b])]
    w = (t >= T0) & (t <= T0+0.10); tw = t[w]-T0; yw = y[w]

    def model(p, tt):
        a1, f1, p1, A, tau, fe, pe = p
        return a1*np.sin(2*np.pi*f1*(tt+T0)+p1) + A*np.exp(-tt/tau)*np.sin(2*np.pi*fe*tt+pe)
    best = None
    for pe in np.linspace(-np.pi, np.pi, 6, endpoint=False):
        p0 = [np.hypot(*ca), 50.0, np.arctan2(ca[1], ca[0]), 30.0, 0.02, fseed, pe]
        rr = least_squares(lambda p: model(p, tw)-yw, p0, bounds=([0, 45, -np.inf, 0, 1e-3, 60, -np.inf], [400, 55, np.inf, 400, 0.2, 500, np.inf]))
        if best is None or rr.cost < best.cost: best = rr
    return float(best.x[5]), float(best.x[4]), float(fseed)


t = tvec(0.6); env = np.where(t >= T0, A0*np.exp(-(t-T0)/TAU), 0.0); out = {}
ref = {'fase_cero': os.path.join(R, 'experiments', 'journal', 'D_ce_raw.jsonl'), 'fase_aleatoria': os.path.join(HERE, 'd_fase_aleatoria_raw.jsonl')}
for caso in ('fase_cero', 'fase_aleatoria'):
    f, tau, seed = [], [], []
    for r in range(N):
        phi = 0.0 if caso == 'fase_cero' else float(np.random.default_rng(1035000+r).uniform(-np.pi, np.pi))
        y = add_noise(base(t) + env*np.sin(2*np.pi*FE*(t-T0) + phi), 30, np.random.default_rng(1032000+r))
        a, b, c = lm_sin_emd(t, y); f.append(a); tau.append(b); seed.append(c)
    f, tau, seed = np.array(f), np.array(tau), np.array(seed); ok = np.abs(f-FE) < 5
    emd = {json.loads(l)['r']: json.loads(l)['EMD-guided LM'] for l in open(ref[caso])}
    fe = np.array([emd[r][0] for r in range(N)]); te = np.array([emd[r][1] for r in range(N)])
    out[caso] = dict(n=N, sin_emd=dict(right=float(ok.mean()), f=summarize(f), tau_err=summarize(100*np.abs(tau-TAU)/TAU), semilla_Hz=summarize(seed)),
                     con_emd=dict(right=float((np.abs(fe-FE) < 5).mean()), f=summarize(fe), tau_err=summarize(100*np.abs(te-TAU)/TAU)),
                     dif_max_f_Hz=float(np.max(np.abs(f-fe))), dif_max_tau_ms=float(1e3*np.max(np.abs(tau-te))))
    v = out[caso]
    print(f"{caso}: sin EMD  correcto {100*v['sin_emd']['right']:.0f} %  f med {v['sin_emd']['f']['med']:.2f}  error de tau med {v['sin_emd']['tau_err']['med']:.2f} %  semilla med {v['sin_emd']['semilla_Hz']['med']:.1f} Hz [{v['sin_emd']['semilla_Hz']['p25']:.1f}-{v['sin_emd']['semilla_Hz']['p75']:.1f}]")
    print(f"{' '*len(caso)}  con EMD  correcto {100*v['con_emd']['right']:.0f} %  f med {v['con_emd']['f']['med']:.2f}  error de tau med {v['con_emd']['tau_err']['med']:.2f} %  | diferencia maxima entre ambos: {v['dif_max_f_Hz']:.4f} Hz, {v['dif_max_tau_ms']:.4f} ms", flush=True)
json.dump(out, open(os.path.join(HERE, 'd_semilla_sin_emd.json'), 'w'), indent=1)
