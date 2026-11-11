"""d_fases.py -- Escenario D: extremos de la senal sin ruido en los primeros 30 ms del transitorio,
para 72 fases iniciales del transitorio (pasos de 5 grados). Suplementario S7.

Senal: v1 + 60 kV * exp(-(t-t0)/15 ms) * sin(2 pi 120 (t-t0) + phi), t0 = 200 ms.
Se cuentan los extremos (cambios de signo de la primera diferencia) en la ventana de 30 ms que
sigue a t0, sobre una rejilla de 1 MHz, y se comparan con los de v1 sin transitorio.

Con phi != 0 el transitorio empieza con un salto de 60 sin(phi) kV, que por si solo puede crear un
extremo en t0. Se reportan dos conteos:
  'tras_inicio'  ventana (t0+0.1 ms, t0+30 ms): no cuenta el extremo del salto (es el conteo de S7);
  'con_salto'    ventana [t0, t0+30 ms): lo cuenta.

Uso: python d_fases.py   (M3_CODE = raiz del repositorio; por defecto, la carpeta superior).
Salida: d_fases.json junto a este archivo."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
from methods.ground_truth import CONFIGS

FS2, T0, A0, TAU, FE, WIN, SALTO = 1e6, 0.20, 60.0, 0.015, 120.0, 0.030, 1e-4
a1, f1, p1, a2, f2, p2, a3, f3, p3 = CONFIGS['v1']
t = np.arange(int(round(0.19*FS2)), int(round(0.24*FS2)))/FS2
YB = a1*np.sin(2*np.pi*f1*t+p1) + a2*np.sin(2*np.pi*f2*t+p2) + a3*np.sin(2*np.pi*f3*t+p3)
ENV = np.where(t >= T0, A0*np.exp(-(t-T0)/TAU), 0.0)


def extremos(y, desde=0.0):
    """Instantes [ms desde t0] de los extremos de y en [t0+desde, t0+WIN)."""
    s = np.sign(np.diff(y)); i = np.where(s[:-1] != s[1:])[0] + 1; tt = t[i]
    return [round(float(x-T0)*1e3, 2) for x in tt[(tt >= T0+desde) & (tt < T0+WIN)]]


N_BASE = len(extremos(YB))          # igual con las dos ventanas: v1 no tiene extremos en t0..t0+0.1 ms
assert N_BASE == len(extremos(YB, SALTO))


def extra_extremos(phi, con_salto=False):
    """Extremos adicionales respecto de v1 sin transitorio, para la fase inicial phi [rad]."""
    return len(extremos(YB + ENV*np.sin(2*np.pi*FE*(t-T0) + phi), 0.0 if con_salto else SALTO)) - N_BASE


def resumen(rows, k):
    return dict(n_fases_con_extra=sum(r[k] > 0 for r in rows),
                histograma={str(v): sum(r[k] == v for r in rows) for v in sorted({r[k] for r in rows})},
                fases_con_extra_deg=[r['fase_deg'] for r in rows if r[k] > 0])


if __name__ == '__main__':
    rows = []
    for deg in range(0, 360, 5):
        y = YB + ENV*np.sin(2*np.pi*FE*(t-T0) + np.deg2rad(deg)); e = extremos(y)
        rows.append(dict(fase_deg=deg, extra_tras_inicio=len(extremos(y, SALTO))-N_BASE, extra_con_salto=len(e)-N_BASE, instantes_ms=e))
    out = dict(rejilla_Hz=FS2, ventana_ms=1e3*WIN, exclusion_salto_ms=1e3*SALTO, extremos_v1_sin_transitorio=N_BASE,
               fase_fundamental_en_t0_deg=round(float(np.degrees((2*np.pi*f1*T0+p1) % (2*np.pi))), 2), n_fases=len(rows),
               tras_inicio=resumen(rows, 'extra_tras_inicio'), con_salto=resumen(rows, 'extra_con_salto'), filas=rows)
    json.dump(out, open(os.path.join(HERE, 'd_fases.json'), 'w'), indent=1)
    print(f"extremos de v1 sin transitorio en la ventana: {N_BASE}; fase de la fundamental en t0: {out['fase_fundamental_en_t0_deg']} grados")
    for k in ('tras_inicio', 'con_salto'):
        print(f"{k}: fases con extremos adicionales {out[k]['n_fases_con_extra']} de {len(rows)}; histograma (extra: n fases) {out[k]['histograma']}")
    print(f"fase 0 (la del articulo): extra = {rows[0]['extra_con_salto']}, instantes {rows[0]['instantes_ms']} ms")
