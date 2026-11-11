"""d_fase_aleatoria.py -- Escenario D con la fase inicial del transitorio sorteada, phi ~ U(-pi, pi).

Prueba exploratoria (11-nov-2026): NO es un resultado del articulo. Responde a la objecion de que
el Escenario D usa una sola fase del transitorio (phi = 0). Mismo ruido y mismas semillas de
CEEMDAN que mc10r_resumable.py (ruido 1032000+r, ensamble 22000+r); la fase usa una semilla aparte
(1035000+r), de modo que cada realizacion es la de D_ce_raw.jsonl con otra fase.

Metodos: los de la Tabla de D que no dependen del tiempo de computo (PF+Hilbert, ESPRIT, Matrix
Pencil) y los dos basados en CEEMDAN (regla de ocupacion de banda del HHT completo y ajuste LM
guiado por EMD), con el codigo de mc4_transient.py. Reanudable: una linea JSONL por realizacion.

Uso: python d_fase_aleatoria.py [control]
  sin argumentos: 100 realizaciones con fase aleatoria -> d_fase_aleatoria_raw.jsonl y .json
  control: las 5 primeras con phi = 0, para comparar con D_ce_raw.jsonl (comprobacion del entorno)."""
import sys, os, json
MODO = sys.argv[1] if len(sys.argv) > 1 else 'aleatoria'; sys.argv = sys.argv[:1]
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R); sys.path.insert(0, HERE)
import numpy as np, warnings; warnings.filterwarnings('ignore')
from mc_common import FS, tvec, base, add_noise, summarize
from methods.extra_estimators import bandpass, analytic_if, esprit
from methods.extra2 import matrix_pencil
import mc4_transient as m4
from d_fases import extra_extremos
from PyEMD import CEEMDAN
N, T0, TAU, FE, A0 = 100, m4.T0, m4.TAU, m4.FE, 60.0
RAW = os.path.join(HERE, 'd_fase_aleatoria_raw.jsonl')


def cerca(f, s):
    k = (f > 90) & (f < 150)
    if not k.any(): return [None, None]
    i = np.argmin(abs(f[k]-120)); return [float(f[k][i]), float(-1/s[k][i]) if s[k][i] < 0 else None]


def una(r, phi):
    t = tvec(0.6); env = np.where(t >= T0, A0*np.exp(-(t-T0)/TAU), 0.0)
    y = add_noise(base(t) + env*np.sin(2*np.pi*FE*(t-T0) + phi), 30, np.random.default_rng(1032000+r))
    row = dict(r=r, phi_deg=float(np.degrees(phi)), extra_tras_inicio=extra_extremos(phi), extra_con_salto=extra_extremos(phi, True))
    a, fi = analytic_if(bandpass(y, FS, 100, 145), FS); row['PF+Hilbert'] = [None if not np.isfinite(v) else float(v) for v in m4.est_env(t, a, fi)]
    seg = (t >= T0+0.002) & (t <= T0+0.10)
    fe, s, _ = esprit(y[seg], FS, 8); row['ESPRIT'] = cerca(fe, s)
    fe, s, _ = matrix_pencil(y[seg], FS, 8); row['Matrix Pencil'] = cerca(fe, s)
    ce = CEEMDAN(trials=30, epsilon=0.05, parallel=False); ce.noise_seed(22000+r); imfs = ce.ceemdan(y, t)
    late = t > T0+0.15; sc = []
    for imf in imfs:
        a, fi = analytic_if(imf, FS); e = (t >= T0) & (t < T0+0.06)
        sc.append(np.mean((fi[e] > 90) & (fi[e] < 150))*min(np.mean(a[(t >= T0) & (t < T0+0.05)])/(np.mean(a[late])+1e-9), 20))
    a, fi = analytic_if(imfs[int(np.argmax(sc))], FS)
    row['HHT, naive IMF rule'] = [None if not np.isfinite(v) else float(v) for v in m4.est_env(t, a, fi)]
    row['EMD-guided LM'] = list(m4.emd_guided(t, y, imfs)); return row


def agregar():
    rows = [json.loads(l) for l in open(RAW)]; out = dict(n=len(rows), nota='Prueba exploratoria; no es un resultado del articulo.')
    ex = np.array([x['extra_tras_inicio'] > 0 for x in rows])
    out['realizaciones_con_extremos_extra_tras_inicio'] = int(ex.sum())
    for k in ['Matrix Pencil', 'ESPRIT', 'EMD-guided LM', 'PF+Hilbert', 'HHT, naive IMF rule']:
        f = np.array([np.nan if x[k][0] is None else x[k][0] for x in rows], float); tau = np.array([np.nan if x[k][1] is None else x[k][1] for x in rows], float)
        ok = np.isfinite(f) & (abs(f-FE) < 5)
        out[k] = dict(right=float(ok.mean()), n_right=int(ok.sum()), right_con_extra=float(ok[ex].mean()) if ex.any() else None,
                      right_sin_extra=float(ok[~ex].mean()) if (~ex).any() else None, f=summarize(f), tau_err=summarize(100*abs(tau-TAU)/TAU))
        print(f"{k:22s} fenomeno correcto {100*ok.mean():5.1f} % ({ok.sum()}/{len(rows)}) | con extremos extra {out[k]['right_con_extra']} | sin {out[k]['right_sin_extra']} | f med {out[k]['f']['med']:.2f} | tau err med {out[k]['tau_err']['med']}")
    json.dump(out, open(os.path.join(HERE, 'd_fase_aleatoria.json'), 'w'), indent=1)


if __name__ == '__main__':
    if MODO == 'control':
        orig = {json.loads(l)['r']: json.loads(l) for l in open(os.path.join(R, 'experiments', 'journal', 'D_ce_raw.jsonl'))}
        for r in range(5):
            x = una(r, 0.0); print(r, 'nuevo', [round(v, 3) for v in x['HHT, naive IMF rule'] + x['EMD-guided LM'] if v is not None],
                                   '| original', [round(v, 3) for v in orig[r]['HHT, naive IMF rule'] + orig[r]['EMD-guided LM'] if v is not None], flush=True)
        sys.exit()
    hechas = {json.loads(l)['r'] for l in open(RAW)} if os.path.exists(RAW) else set()
    for r in range(N):
        if r in hechas: continue
        row = una(r, float(np.random.default_rng(1035000+r).uniform(-np.pi, np.pi)))
        with open(RAW, 'a') as fh: fh.write(json.dumps(row) + '\n')
        print(r, round(row['phi_deg'], 1), row['extra_tras_inicio'], row['HHT, naive IMF rule'][0], row['EMD-guided LM'][0], flush=True)
    agregar()
