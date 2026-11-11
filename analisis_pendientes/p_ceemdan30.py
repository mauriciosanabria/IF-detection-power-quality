"""p_ceemdan30.py -- Escenario P: CEEMDAN con 30 ensayos y 20 realizaciones por celda (11-nov-2026).

En mc6_sweep.py, CEEMDAN usa 20 ensayos y solo las 5 primeras realizaciones de cada celda. Aqui se
repiten las MISMAS 20 realizaciones por celda (ruido y fase con semilla 66, en el mismo orden) con
la configuracion del resto del articulo (30 ensayos, epsilon = 0.05), semilla de ensamble 600 + k.
Misma regla: el IMF mas energetico cuya frecuencia mediana cae en [0.7 f2, 1.3 f2]; exito si el
error de f2 es menor que 1 %.

Control: con 'control' se repiten las 5 primeras realizaciones por celda con 20 ensayos, que deben
reproducir la columna CEEMDAN de res_mc6_sweep.json.

Reanudable: una linea JSONL por descomposicion. Usa varios procesos.
Uso: python p_ceemdan30.py [N_PROCESOS] [control] [agregar]
Salida: p_ceemdan30_raw.jsonl, p_ceemdan30.json (y p_control20_raw.jsonl) junto a este archivo."""
import sys, os, json, math
os.environ.setdefault('OMP_NUM_THREADS', '1'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '1'); os.environ.setdefault('MKL_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
import numpy as np, warnings; warnings.filterwarnings('ignore')
from mc_common import FS, A1, tvec, add_noise
from methods.extra_estimators import analytic_if
from PyEMD import CEEMDAN
RR = [0.0005, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50]; F2 = [70, 100, 150, 250, 400, 650, 1000]; NF = 20
t = tvec(0.5); m = (t > 0.05) & (t < 0.45)
CELDAS = [(r, f2) for r in RR for f2 in F2]


def senales():
    """Las 63 x 20 senales de mc6_sweep.py, en el mismo orden de consumo del generador."""
    rng = np.random.default_rng(66); Y = {}
    for i, (r, f2) in enumerate(CELDAS):
        for k in range(NF):
            Y[(i, k)] = add_noise(A1*np.sin(2*np.pi*50*t) + r*A1*np.sin(2*np.pi*f2*t + rng.uniform(-np.pi, np.pi)), 40, rng)
    return Y


def tarea(a):
    i, k, trials, y = a; r, f2 = CELDAS[i]; band = (0.7*f2, 1.3*f2)
    ce = CEEMDAN(trials=trials, epsilon=0.05, parallel=False); ce.noise_seed(600 + k); imfs = ce.ceemdan(y, t)
    c = [(float(np.std(x)), float(np.median(analytic_if(x, FS)[1][m]))) for x in imfs]; c = [(s, f) for s, f in c if band[0] <= f <= band[1]]
    f = max(c)[1] if c else None
    return dict(id=f'{i}_{k}', i=i, k=k, r=r, f2=f2, err=None if f is None else float(100*abs(f - f2)/f2))


def regimen(r, f2):
    q = f2/50; rho = r*q
    return 'iii' if rho > 1 + 1e-12 else ('i' if r*q*q < 1 - 1e-12 else 'ii')


def qk(r, f2):
    q = f2/50; rho = min(r*q, 1 - 1e-9); c = (1 + q*rho**2)/((1 + q)*rho)
    if c >= 1: return 0.0
    al = math.acos(c); return ((2/math.pi)*((q - 1)*math.atan2(rho*math.sin(al), 1 - rho*math.cos(al)) - al))/(q - 1)


def tasas(path, n):
    rows = [json.loads(l) for l in open(path)]; ok = np.zeros((len(CELDAS), n)); hay = np.zeros((len(CELDAS), n), bool)
    for x in rows:
        if x['k'] < n: ok[x['i'], x['k']] = x['err'] is not None and x['err'] < 1; hay[x['i'], x['k']] = True
    return ok, hay


def agregar():
    orig = {(c['r'], c['f2']): c for c in json.load(open(os.path.join(R, 'experiments', 'journal', 'res_mc6_sweep.json')))}
    ok, hay = tasas(os.path.join(HERE, 'p_ceemdan30_raw.jsonl'), NF); celdas = []
    for i, (r, f2) in enumerate(CELDAS):
        g = regimen(r, f2); n = int(hay[i].sum())
        celdas.append(dict(r=r, f2=f2, rho=r*f2/50, regimen=g, q=1.0 if g == 'iii' else (0.0 if g == 'i' else qk(r, f2)), n=n, exitos=int(ok[i].sum()),
                           tasa=float(ok[i].sum()/n) if n else None, tasa_original_20ensayos_5real=orig[(r, f2)]['CEEMDAN'], tasa_EMD_original=orig[(r, f2)]['EMD']))
    out = dict(ensayos=30, realizaciones_por_celda=NF, celdas=celdas)
    nb = sum(c['n'] for c in celdas if c['regimen'] == 'i'); eb = sum(c['exitos'] for c in celdas if c['regimen'] == 'i'); z = 1.959964
    out['bajo_curvatura'] = dict(realizaciones=nb, exitos=eb, cota_wilson_95_pct=100*z*z/(nb + z*z) if eb == 0 else None)
    for nombre, clave in (('CEEMDAN_30ensayos_20real', 'tasa'), ('CEEMDAN_original', 'tasa_original_20ensayos_5real'), ('EMD_original', 'tasa_EMD_original')):
        s = lambda c: c[clave] is not None and c[clave] >= 0.5; p = lambda c: c['regimen'] == 'iii'
        tp = sum(p(c) and s(c) for c in celdas); fn = sum((not p(c)) and s(c) for c in celdas); tn = sum((not p(c)) and not s(c) for c in celdas); fp = sum(p(c) and not s(c) for c in celdas)
        sub = [c for c in celdas if c['regimen'] != 'i']; trans = [c for c in celdas if c['regimen'] == 'ii']
        out[nombre] = dict(exito_por_regimen={g: sum(s(c) for c in celdas if c['regimen'] == g) for g in ('i', 'ii', 'iii')},
                           realizaciones_con_exito_bajo_curvatura=None if clave != 'tasa' else int(sum(c['exitos'] for c in celdas if c['regimen'] == 'i')),
                           regla_pendiente_acierta=tp + tn, siempre_falla_acierta=sum(not s(c) for c in celdas), acierta_en_25=sum(p(c) == s(c) for c in sub),
                           exactitud_balanceada_63=0.5*(tp/(tp + fn) + tn/(tn + fp)) if tp + fn and tn + fp else None,
                           q_min_con_exito_en_transicion=min([c['q'] for c in trans if s(c)], default=None), q_max_con_fallo_en_transicion=max([c['q'] for c in trans if not s(c)], default=None))
        print(nombre, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in out[nombre].items()})
    cambian = [c for c in celdas if (c['tasa'] >= 0.5) != (c['tasa_original_20ensayos_5real'] >= 0.5)]
    out['celdas_que_cambian_de_clase'] = [dict(r=c['r'], f2=c['f2'], regimen=c['regimen'], q=c['q'], nueva=c['tasa'], original=c['tasa_original_20ensayos_5real']) for c in cambian]
    print('celdas que cambian de clase (exito >= 50 %):', out['celdas_que_cambian_de_clase'])
    P = os.path.join(HERE, 'p_control20_raw.jsonl')
    if os.path.exists(P):
        okc, hayc = tasas(P, 5); dif = [(r, f2) for i, (r, f2) in enumerate(CELDAS) if hayc[i].all() and abs(okc[i].mean() - orig[(r, f2)]['CEEMDAN']) > 1e-9]
        out['control_20ensayos'] = dict(celdas_completas=int(hayc.all(axis=1).sum()), celdas_distintas_del_original=dif); print('control (20 ensayos, 5 realizaciones):', out['control_20ensayos'])
    json.dump(out, open(os.path.join(HERE, 'p_ceemdan30.json'), 'w'), indent=1)


if __name__ == '__main__':
    args = sys.argv[1:]
    if 'agregar' in args: agregar(); sys.exit()
    nproc = int(args[0]) if args and args[0].isdigit() else 8; control = 'control' in args
    raw = os.path.join(HERE, 'p_control20_raw.jsonl' if control else 'p_ceemdan30_raw.jsonl'); trials, n = (20, 5) if control else (30, NF)
    hechas = {json.loads(l)['id'] for l in open(raw)} if os.path.exists(raw) else set(); Y = senales()
    T = [(i, k, trials, Y[(i, k)]) for i in range(len(CELDAS)) for k in range(n) if f'{i}_{k}' not in hechas]
    print(f'{len(T)} descomposiciones pendientes ({trials} ensayos), {nproc} procesos', flush=True)
    from multiprocessing import Pool
    with Pool(nproc) as p:
        for j, row in enumerate(p.imap_unordered(tarea, T, chunksize=4)):
            with open(raw, 'a') as fh: fh.write(json.dumps(row) + '\n')
            if (j + 1) % 100 == 0: print(f'  {j+1}/{len(T)}', flush=True)
    if not control: agregar()
