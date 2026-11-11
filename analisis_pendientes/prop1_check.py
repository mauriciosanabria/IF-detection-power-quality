"""prop1_check.py -- Verificacion numerica de la Proposicion 1 (Ec. 3 del articulo; suplementario S6).

Para y = a1 sin(2 pi f1 t) + ak sin(2 pi fk t + phi), los extremos son los ceros de
g(t) = a1 f1 cos(2 pi f1 t) + ak fk cos(2 pi fk t + phi). Se cuentan los cambios de signo de g
en 20 s sobre una rejilla de 1 MHz y se comparan con
  caso (i)   rho < 1/r      : 2 f1
  caso (ii)  1/r < rho < 1  : Ec. (3), N = 2 f1 [1 + (2/pi)((r-1) arctan(rho sin a/(1-rho cos a)) - a)],
                              cos a = (1 + r rho^2)/((1+r) rho)
  caso (iii) rho > 1        : 2 fk
con rho = ak fk/(a1 f1) y r = fk/f1. Las frecuencias se eligen con relacion no conmensurable.

Uso: python prop1_check.py . No necesita el resto del repositorio. Salida: prop1_check.json."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
F1, T, FS2, PHI = 50.0, 20.0, 1e6, 0.7


def nbar(rho, r):
    if rho >= 1: return 2*F1*r
    if rho <= 1/r: return 2*F1
    al = math.acos((1 + r*rho**2)/((1 + r)*rho))
    return 2*F1*(1 + (2/math.pi)*((r - 1)*math.atan2(rho*math.sin(al), 1 - rho*math.cos(al)) - al))


def conteo(eta, fk, phi=PHI, T=T):
    """Extremos por segundo: cambios de signo de g en T segundos (por bloques de 1 s)."""
    n = 0; last = None
    for b in range(int(T)):
        t = b + np.arange(int(FS2))/FS2
        s = np.sign(F1*np.cos(2*np.pi*F1*t) + eta*fk*np.cos(2*np.pi*fk*t + phi))
        if last is not None and s[0] != last: n += 1
        n += int(np.count_nonzero(s[:-1] != s[1:])); last = s[-1]
    return n/T


# Caso (ii): 14 casos con r = fk/f1 lejos de racionales de orden bajo (parte fraccionaria tipo
# razon aurea), como pide la hipotesis de la Ec. (3). El primero es el citado en S6.
CASOS = {
    'ii': [(0.30, 120.37), (0.45, 91.23), (0.35, 91.23), (0.20, 161.08), (0.12, 161.08), (0.25, 157.08),
           (0.10, 230.90), (0.05, 230.90), (0.19, 230.90), (0.05, 369.45), (0.10, 369.45), (0.02, 665.83),
           (0.05, 665.83), (0.02, 970.71)],
    'i': [(2.56/127.84, 91.23), (0.019/127.84, 1000.22), (0.62/122.44, 161.08), (0.065/127.71, 123.45)],
    'iii': [(0.50, 150.31), (0.20, 401.37), (0.10, 1000.22), (0.40, 250.19)],
}
# Relaciones conmensurables o casi (r racional o cerca de un entero): el conteo depende de phi, y la
# Ec. (3) es su promedio sobre phi (ultima frase del caso ii de la Proposicion 1). Se promedian 60
# fases sorteadas (semilla 7), 5 s cada una.
CASI = [(0.10, 250.19), (0.19, 250.19), (0.05, 401.37), (0.02, 1000.22), (0.30, 120.00)]

if __name__ == '__main__':
    out = {}
    for caso, lista in CASOS.items():
        rows = []
        for eta, fk in lista:
            r = fk/F1; rho = eta*r; pred = nbar(rho, r); obs = conteo(eta, fk)
            assert (caso == 'i' and rho < 1/r) or (caso == 'ii' and 1/r < rho < 1) or (caso == 'iii' and rho > 1)
            rows.append(dict(ak_a1=eta, fk=fk, r=r, rho=rho, q=(pred/(2*F1) - 1)/(r - 1), N_ec3=pred, N_conteo=obs, dif=obs - pred))
            print(f"caso ({caso}) ak/a1={eta:.4f} fk={fk:8.2f} rho={rho:.3f}  Ec.(3)/regla {pred:8.2f}  conteo {obs:8.2f}  dif {obs-pred:+.2f}", flush=True)
        out[caso] = dict(n=len(rows), max_abs_dif=max(abs(x['dif']) for x in rows), filas=rows)
        print(f"  -> caso ({caso}): {len(rows)} casos, |dif| maxima {out[caso]['max_abs_dif']:.2f} extremos/s")
    rows = []
    for eta, fk in CASI:
        r = fk/F1; rho = eta*r; pred = nbar(rho, r); c = [conteo(eta, fk, ph, T=5.0) for ph in np.random.default_rng(7).uniform(0, 2*np.pi, 60)]
        rows.append(dict(ak_a1=eta, fk=fk, r=r, rho=rho, N_ec3=pred, N_min=min(c), N_max=max(c), N_media_phi=float(np.mean(c)), dif_media=float(np.mean(c)) - pred))
        print(f"casi conmensurable ak/a1={eta:.2f} fk={fk:8.2f} rho={rho:.3f}  Ec.(3) {pred:8.2f}  conteo segun phi {min(c):.2f}..{max(c):.2f}  media sobre phi {np.mean(c):8.2f}  dif {np.mean(c)-pred:+.2f}", flush=True)
    out['casi_conmensurables'] = dict(n_fases=60, T_s=5.0, semilla=7, filas=rows)
    # Relacion racional exacta (120 Hz, r = 12/5): media sobre 4000 fases sorteadas, contando en un
    # periodo comun exacto (0.1 s) sobre una rejilla de 2 MHz. Es la cifra citada en S6.
    eta, fk = 0.30, 120.0; tt = np.arange(int(2e5))/2e6; c = []
    for ph in np.random.default_rng(1).uniform(0, 2*np.pi, 4000):
        g = lambda x: F1*np.cos(2*np.pi*F1*x) + eta*fk*np.cos(2*np.pi*fk*x + ph)
        sg = np.sign(g(tt)); c.append(10*(int(np.count_nonzero(sg[:-1] != sg[1:])) + int(sg[-1] != np.sign(g(0.1)))))
    out['racional_120Hz'] = dict(ak_a1=eta, fk=fk, n_fases=4000, semilla=1, N_ec3=nbar(eta*fk/F1, fk/F1), N_media_phi=float(np.mean(c)),
                                 error_estandar=float(np.std(c)/np.sqrt(len(c))), N_min=int(min(c)), N_max=int(max(c)))
    print(f"120 Hz exactos: media sobre 4000 fases {np.mean(c):.2f} +- {np.std(c)/np.sqrt(len(c)):.2f}; Ec.(3) {out['racional_120Hz']['N_ec3']:.2f}; conteos entre {min(c)} y {max(c)}")
    out['ajustes'] = dict(f1=F1, T_s=T, rejilla_Hz=FS2, phi_rad=PHI)
    json.dump(out, open(os.path.join(HERE, 'prop1_check.json'), 'w'), indent=1)
