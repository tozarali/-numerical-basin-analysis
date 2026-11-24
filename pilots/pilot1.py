import numpy as np, json
import genalpha as GA
DC = 0.15
fA = lambda x: -x + x**3
kA = lambda x: -1 + 3*x**2

def ref(X, V, tend, hr=0.01):
    x, v = X.copy(), V.copy()
    rhs = lambda x, v: (v, -DC*v + x - x**3)
    for _ in range(int(tend/hr)):
        k1 = rhs(x, v); k2 = rhs(x+hr/2*k1[0], v+hr/2*k1[1])
        k3 = rhs(x+hr/2*k2[0], v+hr/2*k2[1]); k4 = rhs(x+hr*k3[0], v+hr*k3[1])
        x = x+hr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v+hr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return np.sign(x)

def ga(X, V, tend, h, rho):
    x, v = X.copy(), V.copy(); a = -DC*v - fA(x); t = np.zeros_like(x)
    with np.errstate(all="ignore"):
        for _ in range(int(round(tend/h))):
            x, v, a, _ = GA.step(x, v, a, t, h, fA, kA, DC, lambda tt: 0.0, rho, 1e-10, 30)
    o = np.sign(x); o[~np.isfinite(x)] = 0
    return o

n, tend = 100, 60.0
xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs)
R = ref(X, V, tend)
rhos = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
hs = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
M = np.zeros((len(hs), len(rhos)))
for i, h in enumerate(hs):
    for j, r in enumerate(rhos):
        M[i, j] = np.mean(ga(X, V, tend, h, r) != R)
print("rho:", rhos)
for h, row in zip(hs, M):
    print(f"h={h:<5}", " ".join(f"{100*m:5.1f}" for m in row), " argmin rho =", rhos[int(np.argmin(row))])
# robustness: tend and reference step
R2 = ref(X, V, 120.0)
print("ref disagreement tend 60 vs 120:", np.mean(R != R2))
