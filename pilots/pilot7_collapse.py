"""Collapse test: does basin mismatch organise better against the energy-loss ratio D than against rho_inf?
Pools gen-alpha, HHT, WBZ, Newmark. D = median(dE_num/dE_ref) over ICs with 0.6<E0<2 after T=8."""
import numpy as np
from scipy.stats import spearmanr
from pilot5_families import step, param_for_rho, make_system, ref, run

def E(x, v): return 0.5*v**2 + 0.25*x**4 - 0.5*x**2 + 0.25

def ref_state(X, V, T, DC, hr=0.005):
    x, v = X.copy(), V.copy(); rhs = lambda x, v: (v, -DC*v + x - x**3)
    for _ in range(int(round(T/hr))):
        k1 = rhs(x, v); k2 = rhs(x+hr/2*k1[0], v+hr/2*k1[1]); k3 = rhs(x+hr/2*k2[0], v+hr/2*k2[1]); k4 = rhs(x+hr*k3[0], v+hr*k3[1])
        x = x+hr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v+hr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return x, v

def ga_state(X, V, T, h, p, DC):
    f, k = make_system(0.0); x, v = X.copy(), V.copy(); a = -DC*v - f(x)
    with np.errstate(all="ignore"):
        for _ in range(int(round(T/h))): x, v, a = step(x, v, a, h, f, k, DC, p)
    return x, v

n, tend, T = 100, 60.0, 8.0
xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs)
sel = (E(X, V) > 0.6) & (E(X, V) < 2.0); Xs, Vs = X[sel], V[sel]; E0 = E(Xs, Vs)
targets = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
for DC, h in [(0.15, 0.4), (0.05, 0.2)]:
    R = ref(X, V, tend, DC, 0.0); xr, vr = ref_state(Xs, Vs, T, DC); dref = E0 - E(xr, vr)
    rows = []
    for name in ["genalpha", "HHT", "WBZ", "Newmark"]:
        for t in targets:
            p, r = param_for_rho(name, t)
            if p is None: continue
            mis = 100*np.mean(run(X, V, tend, h, p, DC, 0.0) != R)
            xg, vg = ga_state(Xs, Vs, T, h, p, DC); D = np.nanmedian((E0 - E(xg, vg))/dref)
            rows.append((name, t, D, mis))
    print(f"\n=== c={DC}, h={h} ===\nfamily     rho_inf     D   mismatch%")
    for r in rows: print(f"{r[0]:<10} {r[1]:5.2f} {r[2]:6.2f} {r[3]:7.1f}")
    arr = np.array([(r[1], r[2], r[3]) for r in rows])
    print("Spearman(mismatch, |D-1|)      = %.2f" % spearmanr(np.abs(arr[:, 1]-1), arr[:, 2])[0])
    print("Spearman(mismatch, |rho-rho*|) uses best single-family rho*: n/a; Spearman(mismatch, 1-rho) = %.2f" % spearmanr(1-arr[:, 0], arr[:, 2])[0])
    # best-D-match per family
    for name in ["genalpha", "HHT", "WBZ", "Newmark"]:
        fr = [r for r in rows if r[0] == name]; i = int(np.argmin([abs(r[2]-1) for r in fr])); j = int(np.argmin([r[3] for r in fr]))
        print(f"{name:<9} D-closest-to-1: rho={fr[i][1]:.2f} (mis {fr[i][3]:.1f}%)   basin-best: rho={fr[j][1]:.2f} (mis {fr[j][3]:.1f}%)")
