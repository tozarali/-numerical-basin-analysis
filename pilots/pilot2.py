import numpy as np
import genalpha as GA
DC, K = 0.15, 2.0   # linearisation at the wells (f'(+-1)=2)
f = lambda x: K*x
k = lambda x: K*np.ones_like(x)

def multiplier(h, rho):
    # build the one-step map on (x,v,a) by probing with unit vectors (linear problem)
    M = np.zeros((3, 3))
    for j in range(3):
        s = np.zeros(3); s[j] = 1.0
        x, v, a, _ = GA.step(np.array(s[0]), np.array(s[1]), np.array(s[2]), 0.0, h, f, k, DC, lambda t: 0.0, rho, 1e-14, 5)
        M[:, j] = [float(x), float(v), float(a)]
    ev = np.linalg.eigvals(M)
    ev = ev[np.argmax(np.abs(ev.imag))]  # dominant oscillatory pair
    return ev

w0 = np.sqrt(K); xi = DC/(2*w0); wd = w0*np.sqrt(1-xi**2)
rhos = np.round(np.arange(0, 1.01, 0.1), 1)
hs = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
print("exact: decay rate sigma=%.4f, wd=%.4f" % (xi*w0, wd))
for h in hs:
    rows = []
    for r in rhos:
        z = multiplier(h, r)
        sig = -np.log(abs(z))/h; w = abs(np.angle(z))/h
        rows.append((r, sig/(xi*w0), w/wd))
    print("h=%.1f  rho: decay-ratio / freq-ratio" % h)
    print("   " + "  ".join("%.1f:%5.2f/%5.3f" % t for t in rows))
