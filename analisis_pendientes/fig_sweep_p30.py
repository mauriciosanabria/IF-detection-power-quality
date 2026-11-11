"""fig_sweep_p30.py -- Fig. 2 del articulo (mapa del barrido P) con el panel de CEEMDAN de la corrida
p_ceemdan30.py (30 ensayos, 20 realizaciones por celda). Los otros cinco paneles salen de
res_mc6_sweep.json sin cambios. Mismo dibujo que experiments/journal/fig_sweep.py.

Uso: python fig_sweep_p30.py   (despues de p_ceemdan30.py)
Salida: res_mc6_sweep_p30.json (copia de res_mc6_sweep.json con la columna CEEMDAN nueva) y
sweep_map.pdf, junto a este archivo. No sobrescribe figures/sweep_map.pdf."""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
R_ = os.environ.get('M3_CODE', os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(R_, 'methods'))
from ieee_style import plt
import numpy as np
d = json.load(open(os.path.join(R_, 'experiments', 'journal', 'res_mc6_sweep.json')))
nueva = {(c['r'], c['f2']): c for c in json.load(open(os.path.join(HERE, 'p_ceemdan30.json')))['celdas']}
for x in d:
    c = nueva[(x['r'], x['f2'])]; assert c['n'] == 20
    x['CEEMDAN_20ensayos_5real'] = x['CEEMDAN']; x['CEEMDAN'] = c['tasa']
json.dump(d, open(os.path.join(HERE, 'res_mc6_sweep_p30.json'), 'w'), indent=1)
R = sorted({x['r'] for x in d}); F = sorted({x['f2'] for x in d})
meths = ['IpDFT', 'ESPRIT', 'PF+Hilbert', 'VMD', 'EMD', 'CEEMDAN']
fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.6), sharex=True, sharey=True)
for ax, m in zip(axs.flat, meths):
    Z = np.array([[next(x[m] for x in d if x['r'] == r and x['f2'] == f) for f in F] for r in R])
    ax.pcolormesh(np.arange(len(F)+1), np.arange(len(R)+1), Z, cmap='Greys', vmin=0, vmax=1, shading='flat')
    for i, r in enumerate(R):
        for j, f in enumerate(F):
            ax.text(j+0.5, i+0.5, f'{100*Z[i, j]:.0f}', ha='center', va='center', fontsize=5.2, color='white' if Z[i, j] > 0.55 else 'black')
    # rho = 1 boundary: r = 50/f
    xs = np.linspace(0, len(F), 200); fx = np.interp(xs-0.5, np.arange(len(F)), np.log(F)); rr = np.log(50/np.exp(fx))
    ys = np.interp(rr, np.log(R), np.arange(len(R))) + 0.5
    ok = (rr >= np.log(R[0])) & (rr <= np.log(R[-1])); ax.plot(xs[ok], ys[ok], color='#c0392b', lw=1.3, ls='--')
    ax.set_title(m, fontsize=8.5)
for ax in axs[-1]: ax.set_xticks(np.arange(len(F))+0.5); ax.set_xticklabels(F, fontsize=6.5); ax.set_xlabel('$f_2$ (Hz)', fontsize=8)
for ax in axs[:, 0]: ax.set_yticks(np.arange(len(R))+0.5); ax.set_yticklabels([f'{100*r:g}' for r in R], fontsize=6.5); ax.set_ylabel('$a_2/a_1$ (%)', fontsize=8)
for ax in axs.flat: ax.tick_params(direction='out', top=False, right=False)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'sweep_map.pdf')); plt.savefig(os.path.join(HERE, 'sweep_map_p30.png'), dpi=110)
print('ok')
