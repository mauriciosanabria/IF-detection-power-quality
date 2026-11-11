"""ceemdan_ensayos_piloto.py -- PILOTO (11-nov-2026): efecto del numero de ensayos de CEEMDAN y de
recombinar los IMF de la banda, en el Escenario S (v1 y v3, 60 dB). NO es un resultado del articulo.

Responde a la objecion del "hombre de paja": con epsilon = 0.05 y 30 ensayos el residuo del ensamble
queda por encima de los tonos pequenos, y el pipeline toma un solo IMF. Aqui se repiten las primeras
realizaciones de mc1_stationary.py (mismo ruido: semilla 100*SNR + indice del nodo; misma semilla de
ensamble 1000 + r) con 30 y con mas ensayos, y con dos reglas:
  'un IMF'  : el IMF mas energetico con frecuencia mediana en la banda (regla del articulo);
  'suma'    : la suma de todos los IMF con frecuencia mediana en la banda.
Exito: error relativo de f2 < 1 %.

Uso: python ceemdan_ensayos_piloto.py [N_REAL] [ENSAYOS ...]   (por defecto 5 realizaciones; 30 y 200)
Salida: ceemdan_ensayos_piloto.json junto a este archivo."""
import sys, os, json, time
ARGS = [int(a) for a in sys.argv[1:]]; sys.argv = sys.argv[:1]
NR = ARGS[0] if ARGS else 5; ENS = ARGS[1:] or [30, 200]
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
import numpy as np, warnings; warnings.filterwarnings('ignore')
from mc_common import FS, add_noise
from methods.ground_truth import CONFIGS, y_model, time_vector
from methods.extra_estimators import analytic_if
from PyEMD import CEEMDAN
SNR = 60; t = time_vector(); m = (t > 0.05) & (t < 0.45); out = {}
for node in ['v1', 'v3']:
    th = CONFIGS[node]; f2 = th[4]; b2 = (max(60, f2*0.5), f2*1.5); y0 = y_model(th, t)
    rng = np.random.default_rng(100*SNR + list(CONFIGS).index(node)); ys = [add_noise(y0, SNR, rng) for _ in range(NR)]
    for trials in ENS:
        res = {'un IMF': [], 'suma': [], 'n_imf_en_banda': [], 'rms_imf1_kV': []}; t0 = time.time()
        for r, y in enumerate(ys):
            c = CEEMDAN(trials=trials, epsilon=0.05, parallel=False); c.noise_seed(1000+r); imfs = c.ceemdan(y, t)
            med = [float(np.median(analytic_if(x, FS)[1][m])) for x in imfs]; inb = [i for i, f in enumerate(med) if b2[0] <= f <= b2[1]]
            e = lambda f: None if f is None else float(100*abs(f-f2)/f2)
            res['un IMF'].append(e(med[max(inb, key=lambda i: np.std(imfs[i]))] if inb else None))
            res['suma'].append(e(float(np.median(analytic_if(sum(imfs[i] for i in inb), FS)[1][m])) if inb else None))
            res['n_imf_en_banda'].append(len(inb)); res['rms_imf1_kV'].append(float(np.sqrt(np.mean(imfs[0]**2))))
        ok = lambda v: sum(x is not None and x < 1 for x in v)
        out[f'{node}_{trials}'] = dict(nodo=node, ensayos=trials, n=NR, exito_un_imf=ok(res['un IMF']), exito_suma=ok(res['suma']), s_por_descomposicion=(time.time()-t0)/NR, **res)
        print(f"{node} {SNR} dB, {trials:4d} ensayos: exito un IMF {ok(res['un IMF'])}/{NR}, suma {ok(res['suma'])}/{NR}; IMF en banda {res['n_imf_en_banda']}; "
              f"RMS del primer IMF {np.mean(res['rms_imf1_kV']):.3f} kV; {(time.time()-t0)/NR:.1f} s por descomposicion", flush=True)
        json.dump(out, open(os.path.join(HERE, 'ceemdan_ensayos_piloto.json'), 'w'), indent=1)
