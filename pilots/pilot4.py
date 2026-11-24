import numpy as np
import genalpha as GA
fA = lambda x: -x + x**3
kA = lambda x: -1 + 3*x**2
E = lambda x, v: 0.5*v**2 + 0.25*x**4 - 0.5*x**2 + 0.25

def ref(X, V, T, DC, hr=0.005):
    x, v = X.copy(), V.copy()
    rhs = lambda x, v: (v, -DC*v + x - x**3)
    for _ in range(int(round(T/hr))):
        k1 = rhs(x, v); k2 = rhs(x+hr/2*k1[0], v+hr/2*k1[1])
        k3 = rhs(x+hr/2*k2[0], v+hr/2*k2[1]); k4 = rhs(x+hr*k3[0], v+hr*k3[1])
        x = x+hr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v+hr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return x, v

def ga(X, V, T, h, rho, DC):
    x, v = X.copy(), V.copy(); a = -DC*v - fA(x); t = np.zeros_like(x)
    with np.errstate(all="ignore"):
        for _ in range(int(round(T/h))):
            x, v, a, _ = GA.step(x, v, a, t, h, fA, kA, DC, lambda tt: 0.0, rho, 1e-10, 30)
    return x, v

n = 70; xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs); X, V = X.ravel(), V.ravel()
E0 = E(X, V); sel = (E0 > 0.6) & (E0 < 2.0); X, V, E0 = X[sel], V[sel], E0[sel]
T = 8.0
rhos = np.round(np.arange(0, 1.001, 0.05), 2)
obs = {0.05: {0.2: .53, 0.3: .58, 0.4: .61, 0.5: .80, 0.6: .77}, 0.15: {0.2: .37, 0.3: .40, 0.4: .40, 0.5: .39, 0.6: .37}}
for DC in [0.05, 0.15]:
    xr, vr = ref(X, V, T, DC); dref = E0 - E(xr, vr)
    print(f"\n=== c={DC}  (N={len(X)} ICs, T={T}) ===")
    for h in [0.2, 0.3, 0.4, 0.5, 0.6]:
        D = []
        for r in rhos:
            xg, vg = ga(X, V, T, h, r, DC); dg = E0 - E(xg, vg)
            D.append(np.nanmedian(dg/dref))
        D = np.array(D); idx = np.where(np.diff(np.sign(D - 1)) != 0)[0]
        rs = float("nan")
        if len(idx):
            i = idx[0]; rs = rhos[i] + (1 - D[i])/(D[i+1]-D[i])*0.05
        print(f"h={h}: D(rho=0,.25,.5,.75,1)=", " ".join(f"{D[k]:.2f}" for k in (0, 5, 10, 15, 20)), f"| rho*(energy)={rs:.2f}  basin argmin={obs[DC][h]:.2f}")
