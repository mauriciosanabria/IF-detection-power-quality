"""vmd_n.py -- (d) mecanismo del error de VMD en el Escenario N (fundamental con muescas bajo rampa)."""
import sys, os, numpy as np, json
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.environ.get('M3_CODE', os.path.join(HERE, '..'))  # raiz del paquete de codigo del repositorio
sys.path.insert(0, os.path.join(R, 'experiments', 'journal')); sys.path.insert(0, R)
from mc_common import FS, rocof_signal, add_noise, CONFIGS
from methods.extra_estimators import vmd_modes, analytic_if
t, y0, ft = rocof_signal(); a1, _, p1 = CONFIGS['v1'][:3]
ph = np.mod(2*np.pi*np.cumsum(ft)/FS + p1, np.pi/3); h = (ph > 0.5) & (ph < 0.5 + 2*np.pi*50*0.0006)
yN = y0*(1 - 0.3*h)
def rm(f, w): m = (t >= w[0]) & (t <= w[1]); return float(np.sqrt(np.mean((f[m]-ft[m])**2)))
out = {}
for lab, base in [('R (sin muescas)', y0), ('N (con muescas)', yN)]:
    rng = np.random.default_rng(950); y = add_noise(base, 30, rng)
    u, c = vmd_modes(y, FS, 3); i = int(np.argmin(abs(c-50))); m1 = u[i]; f = analytic_if(m1, FS)[1]
    S = np.abs(np.fft.rfft(m1*np.hanning(len(m1))))**2; fr = np.fft.rfftfreq(len(m1), 1/FS)
    tot = S.sum(); fund = S[(fr > 40) & (fr < 60)].sum(); hi = S[fr > 100].sum()
    r = (t > .6) & (t < .8); d = f[r]-ft[r]; D = np.abs(np.fft.rfft(d-d.mean())); fd = np.fft.rfftfreq(len(d), 1/FS)
    # contenido de la entrada: energia de las muescas fuera de la fundamental y donde cae
    Sy = np.abs(np.fft.rfft((base)*np.hanning(len(base))))**2
    print(f"{lab}: centros VMD {np.round(np.sort(c),1)} Hz | modo fundamental: energia >100 Hz = {100*hi/tot:.2f} % | "
          f"error IF rampa {rm(f,(.4,.55)):.3f} Hz, asentado {rm(f,(.6,.8)):.3f} Hz | rizo del IF dominante a {fd[np.argmax(D[1:])+1]:.0f} Hz, p2p {np.ptp(d):.2f} Hz")
    # componentes de la entrada en armonicos de la muesca
    arm = {}
    for hh in (5, 7, 11, 13):
        b = (fr > hh*50-15) & (fr < hh*50+15); e_in = 100*np.sqrt(Sy[b].sum()/Sy[(fr>40)&(fr<60)].sum()); e_mo = 100*np.sqrt(S[b].sum()/fund)
        print(f"     entrada h{hh}: {e_in:.2f} % de la fundamental; en el modo: {e_mo:.2f} %")
        # rizo esperado del IF por un residuo de fraccion a'/a1 a la frecuencia fk: (a'/a1)(fk - f1), con f1 asentada en 49.4 Hz
        arm[f'h{hh}'] = dict(entrada_pct=float(e_in), en_modo_pct=float(e_mo), rizo_esperado_Hz=float(e_mo/100*(hh*50-49.4)))
    out[lab] = dict(centros=np.sort(c).tolist(), hi=float(hi/tot), ramp=rm(f, (.4, .55)), settled=rm(f, (.6, .8)), ripple=float(fd[np.argmax(D[1:])+1]),
                    rizo_p2p_Hz=float(np.ptp(d)), ventana_rizo_s=[0.6, 0.8], armonicos=arm)
json.dump(out, open(os.path.join(HERE, 'vmd_n.json'), 'w'), indent=1)
