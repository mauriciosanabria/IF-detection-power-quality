"""ceemdan_ensayos.py -- CEEMDAN con mas ensayos y con recombinacion de los IMF de la banda,
en los Escenarios S y R (11-nov-2026). Responde a la objecion del "hombre de paja": con
epsilon = 0.05 y 30 ensayos el residuo del ensamble queda por encima de los tonos pequenos, y el
pipeline del articulo toma un solo IMF.

Se repiten exactamente las realizaciones del articulo:
  S  mc1_stationary.py: ruido con semilla 100*SNR + indice del nodo, ensamble 1000 + r, 15 realizaciones
     por condicion (v1, v3, v2 a 30 y 60 dB); banda [max(60, 0.5 f2), 1.5 f2].
  R  mc10_ceemdan100.py: ruido con semilla 1030, ensamble 20000 + r, 100 realizaciones a 30 dB.
con 30, 200 y 500 ensayos, y con dos reglas de seleccion:
  'un_imf'  la del articulo. S: el IMF mas energetico con frecuencia mediana en la banda.
            R: el IMF cuya frecuencia mediana (200-800 ms) esta mas cerca de 50 Hz.
  'suma'    la suma de todos los IMF con frecuencia mediana en la banda (S: la del tono; R: 40-60 Hz).
Con 30 ensayos, 'un_imf' debe reproducir el articulo (S: 0 %; R: 1.568 / 1.523 Hz).

Reanudable: una linea JSONL por descomposicion. Usa varios procesos.
Uso: python ceemdan_ensayos.py [N_PROCESOS] [agregar]
Salida: ceemdan_ensayos_raw.jsonl y ceemdan_ensayos.json junto a este archivo."""
import sys, os, json, time
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('MKL_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
import numpy as np, warnings; warnings.filterwarnings('ignore')
from mc_common import FS, rocof_signal, add_noise, summarize
from methods.ground_truth import CONFIGS, y_model, time_vector
from methods.extra_estimators import analytic_if
from PyEMD import CEEMDAN
RAW = os.path.join(HERE, 'ceemdan_ensayos_raw.jsonl'); ENSAYOS = [30, 200, 500]; NS, NR = 15, 100


def senal_S(snr, node, r):
    t = time_vector(); rng = np.random.default_rng(100*snr + list(CONFIGS).index(node)); y0 = y_model(CONFIGS[node], t)
    for _ in range(r + 1): y = add_noise(y0, snr, rng)
    return t, y


def senal_R(r):
    t, y0, ft = rocof_signal(); rng = np.random.default_rng(1030)
    for _ in range(r + 1): y = add_noise(y0, 30, rng)
    return t, y, ft


def rmse(f, ft, t, w):
    m = (t >= w[0]) & (t <= w[1]) & np.isfinite(f); return float(np.sqrt(np.mean((f[m]-ft[m])**2)))


def tarea(a):
    t0 = time.time()
    if a[0] == 'S':
        _, snr, node, r, trials = a; t, y = senal_S(snr, node, r); f2 = CONFIGS[node][4]; b2 = (max(60, f2*0.5), f2*1.5); m = (t > 0.05) & (t < 0.45)
        c = CEEMDAN(trials=trials, epsilon=0.05, parallel=False); c.noise_seed(1000 + r); imfs = c.ceemdan(y, t)
        med = [float(np.median(analytic_if(x, FS)[1][m])) for x in imfs]; inb = [i for i, f in enumerate(med) if b2[0] <= f <= b2[1]]
        e = lambda f: float(100*abs(f - f2)/f2)
        row = dict(un_imf=e(med[max(inb, key=lambda i: np.std(imfs[i]))]) if inb else None,
                   suma=e(float(np.median(analytic_if(sum(imfs[i] for i in inb), FS)[1][m]))) if inb else None)
    else:
        _, r, trials = a; t, y, ft = senal_R(r); w = (t > .2) & (t < .8)
        c = CEEMDAN(trials=trials, epsilon=0.05, parallel=False); c.noise_seed(20000 + r); imfs = c.ceemdan(y, t)
        fis = [analytic_if(x, FS)[1] for x in imfs]; med = [float(np.median(f[w])) for f in fis]; inb = [i for i, f in enumerate(med) if 40 <= f <= 60]
        best = fis[int(np.argmin([abs(f - 50) for f in med]))]; row = dict(un_imf=[rmse(best, ft, t, (.4, .55)), rmse(best, ft, t, (.6, .8))])
        if inb:
            fs_ = analytic_if(sum(imfs[i] for i in inb), FS)[1]; row['suma'] = [rmse(fs_, ft, t, (.4, .55)), rmse(fs_, ft, t, (.6, .8))]
        else: row['suma'] = None
    row.update(id='_'.join(str(x) for x in a), n_imfs=len(imfs), n_en_banda=len(inb), rms_imf1=float(np.sqrt(np.mean(imfs[0]**2))), s=round(time.time() - t0, 1)); return row


def agregar():
    rows = {json.loads(l)['id']: json.loads(l) for l in open(RAW)}; out = {'S': {}, 'R': {}}
    for trials in ENSAYOS:
        for snr in (30, 60):
            for node in ('v1', 'v3', 'v2'):
                v = [rows[f'S_{snr}_{node}_{r}_{trials}'] for r in range(NS) if f'S_{snr}_{node}_{r}_{trials}' in rows]
                if not v: continue
                ok = lambda k: int(sum(x[k] is not None and x[k] < 1 for x in v))
                out['S'][f'{snr}_{node}_{trials}'] = dict(n=len(v), exito_un_imf=ok('un_imf'), exito_suma=ok('suma'), realizaciones_con_imf_en_banda=int(sum(x['n_en_banda'] > 0 for x in v)),
                                                         rms_imf1_kV=float(np.mean([x['rms_imf1'] for x in v])))
        v = [rows[f'R_{r}_{trials}'] for r in range(NR) if f'R_{r}_{trials}' in rows]
        if v:
            d = dict(n=len(v), rms_imf1_kV=float(np.mean([x['rms_imf1'] for x in v])), s_por_descomposicion=float(np.mean([x['s'] for x in v])))
            for k in ('un_imf', 'suma'):
                w = [x[k] for x in v if x[k] is not None]
                d[k] = dict(ramp=summarize([a for a, b in w]), settled=summarize([b for a, b in w]), sin_imf_en_banda=len(v) - len(w))
            out['R'][str(trials)] = d
    json.dump(out, open(os.path.join(HERE, 'ceemdan_ensayos.json'), 'w'), indent=1)
    for k, v in out['S'].items(): print('S', k, f"un IMF {v['exito_un_imf']}/{v['n']}  suma {v['exito_suma']}/{v['n']}  con IMF en banda {v['realizaciones_con_imf_en_banda']}  RMS IMF1 {v['rms_imf1_kV']:.2f} kV")
    for k, v in out['R'].items(): print('R', k, f"n {v['n']}  un IMF rampa {v['un_imf']['ramp']['med']:.3f} asentado {v['un_imf']['settled']['med']:.3f} | suma rampa {v['suma']['ramp']['med']:.3f} asentado {v['suma']['settled']['med']:.3f} (sin IMF en banda: {v['suma']['sin_imf_en_banda']})  RMS IMF1 {v['rms_imf1_kV']:.2f} kV")


if __name__ == '__main__':
    args = sys.argv[1:]
    if 'agregar' in args: agregar(); sys.exit()
    nproc = int(args[0]) if args else 8
    T = [('S', snr, node, r, tr) for tr in ENSAYOS for snr in (30, 60) for node in ('v1', 'v3', 'v2') for r in range(NS)] + [('R', r, tr) for tr in ENSAYOS for r in range(NR)]
    hechas = {json.loads(l)['id'] for l in open(RAW)} if os.path.exists(RAW) else set()
    T = [a for a in T if '_'.join(str(x) for x in a) not in hechas]; print(f'{len(T)} descomposiciones pendientes, {nproc} procesos', flush=True)
    from multiprocessing import Pool
    with Pool(nproc) as p:
        for i, row in enumerate(p.imap_unordered(tarea, T)):
            with open(RAW, 'a') as fh: fh.write(json.dumps(row) + '\n')
            if (i + 1) % 25 == 0: print(f'  {i+1}/{len(T)}', flush=True)
    agregar()
