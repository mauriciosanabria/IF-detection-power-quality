"""p_frontera_rf.py -- Escenario P: frontera refinada de Rilling y Flandrin (2008) frente a la
frontera de pendiente (11-nov-2026). No corre ninguna descomposicion: usa res_mc6_sweep.json (EMD)
y p_ceemdan30.json (CEEMDAN, 30 ensayos, 20 realizaciones).

Rilling y Flandrin (IEEE Trans. Signal Process. 56(1), 2008), con su modelo
cos 2 pi t + a cos(2 pi f t + phi), f < 1, es decir a = a1/ak y f = f1/fk = 1/r:
  - p. 87-88: frecuencia de corte f_c ~ 0.67 (10 iteraciones de cribado), por encima de la cual los
    dos tonos no se separan con ninguna amplitud; aqui r > 1/0.67 ~ 1.5.
  - Sec. IV-C, p. 91-92: para f < 1/3 la curva af = 1 es holgada y la frontera mas ajustada es
    a f sin(3 pi f / 2) = 1; aqui rho > sin(3 pi / (2 r)) para r > 3.
Reglas comparadas (prediccion de exito, tasa >= 50 %):
  pendiente : rho > 1
  RF        : r >= 1/0.67 y (rho > 1 si r <= 3; rho > sin(3 pi/(2 r)) si r > 3)
Salida: p_frontera_rf.json junto a este archivo."""
import os, json, math
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
d = json.load(open(os.path.join(R, 'experiments', 'journal', 'res_mc6_sweep.json')))
ce = {(c['r'], c['f2']): c['tasa'] for c in json.load(open(os.path.join(HERE, 'p_ceemdan30.json')))['celdas']}
FC = 0.67


def pendiente(r_amp, f2): return r_amp*f2/50 > 1 + 1e-12


def rf(r_amp, f2):
    r = f2/50; rho = r_amp*r
    if 1/r > FC: return False
    return rho > (1 + 1e-12 if r <= 3 + 1e-12 else math.sin(3*math.pi/(2*r)))


out = {'f_c': FC, 'r_minimo': 1/FC}
for nombre, tasa in (('EMD', lambda c: c['EMD']), ('CEEMDAN', lambda c: ce[(c['r'], c['f2'])])):
    out[nombre] = {}
    for regla, fn in (('pendiente', pendiente), ('RF', rf)):
        mal = [dict(r=c['r'], f2=c['f2'], rho=c['r']*c['f2']/50, tasa=tasa(c), prediccion=bool(fn(c['r'], c['f2']))) for c in d if fn(c['r'], c['f2']) != (tasa(c) >= 0.5)]
        out[nombre][regla] = dict(aciertos=len(d) - len(mal), de=len(d), fallos=mal)
        print(f"{nombre:8s} regla {regla:9s}: acierta {len(d)-len(mal)} de {len(d)}; falla en {[(100*m['r'], m['f2'], m['tasa']) for m in mal]}")
c70 = [c for c in d if c['f2'] == 70]
out['celdas_70Hz'] = dict(n=len(c70), EMD_max=max(c['EMD'] for c in c70), CEEMDAN_max=max(ce[(c['r'], c['f2'])] for c in c70), VMD_max=max(c['VMD'] for c in c70))
print('70 Hz (r = 1.4, bajo el corte):', out['celdas_70Hz'])
json.dump(out, open(os.path.join(HERE, 'p_frontera_rf.json'), 'w'), indent=1)
