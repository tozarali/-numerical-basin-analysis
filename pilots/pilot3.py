import numpy as np
import genalpha as GA
fA = lambda x: -x + x**3
kA = lambda x: -1 + 3*x**2

def ref(X, V, tend, DC, hr=0.01):
    x, v = X.copy(), V.copy()
    rhs = lambda x, v: (v, -DC*v + x - x**3)
    for _ in range(int(tend/hr)):
        k1 = rhs(x, v); k2 = rhs(x+hr/2*k1[0], v+hr/2*k1[1])
        k3 = rhs(x+hr/2*k2[0], v+hr/2*k2[1]); k4 = rhs(x+hr*k3[0], v+hr*k3[1])
        x = x+hr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v+hr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return np.sign(x)

def ga(X, V, tend, h, rho, DC):
    x, v = X.copy(), V.copy(); a = -DC*v - fA(x); t = np.zeros_like(x)
    with np.errstate(all="ignore"):
        for _ in range(int(round(tend/h))):
            x, v, a, _ = GA.step(x, v, a, t, h, fA, kA, DC, lambda tt: 0.0, rho, 1e-10, 30)
    o = np.sign(x); o[~np.isfinite(x)] = 0
    return o

def decay_ratio(h, rho, DC, K=2.0):
    f = lambda x: K*x; k = lambda x: K*np.ones_like(x)
    M = np.zeros((3, 3))
    for j in range(3):
        s = np.zeros(3); s[j] = 1.0
        x, v, a, _ = GA.step(np.array(s[0]), np.array(s[1]), np.array(s[2]), 0.0, h, f, k, DC, lambda t: 0.0, rho, 1e-14, 5)
        M[:, j] = [float(x), float(v), float(a)]
    ev = np.linalg.eigvals(M); z = ev[np.argmax(np.abs(ev.imag))]
    return (-np.log(abs(z))/h) / (DC/2)

def rho_star(h, DC):
    rs = np.linspace(0, 1, 101); d = np.array([decay_ratio(h, r, DC) for r in rs])
    idx = np.where(np.diff(np.sign(d - 1)) != 0)[0]
    if len(idx) == 0: return float("nan")
    i = idx[0]; return rs[i] + (1 - d[i])/(d[i+1]-d[i])*0.01

n, tend = 100, 80.0
xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs)
rhos = np.round(np.arange(0, 1.01, 0.1), 1)
for DC in [0.05, 0.15, 0.4]:
    R = ref(X, V, tend, DC)
    print(f"\n=== c = {DC} ===  rho: " + " ".join(f"{r:5.1f}" for r in rhos))
    for h in [0.2, 0.3, 0.4, 0.5, 0.6]:
        m = np.array([100*np.mean(ga(X, V, tend, h, r, DC) != R) for r in rhos])
        # parabola-refined argmin around the discrete minimum
        j = int(np.argmin(m)); rb = rhos[j]
        if 0 < j < len(rhos)-1:
            y0, y1, y2 = m[j-1], m[j], m[j+1]; den = y0 - 2*y1 + y2
            rb = rhos[j] + 0.1*0.5*(y0 - y2)/den if den > 0 else rhos[j]
        print(f"h={h:<4}", " ".join(f"{v:5.1f}" for v in m), f" | argmin~{rb:.2f}  predicted rho*={rho_star(h, DC):.2f}")
