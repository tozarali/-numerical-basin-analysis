"""H5 test: same basin-error-vs-rho_inf curve across integrator families (gen-alpha, HHT, WBZ, Newmark),
matched by the numerically measured high-frequency spectral radius rho_inf. Also tilted double well (asymmetry)."""
import numpy as np

def step(x0, v0, a0, h, f, k, d, p, tol=1e-10, maxit=30, force_iters=None):
    am, af, g, b = p
    c1 = 1.0/(b*h*h); x1 = x0 + h*v0 + 0.5*h*h*a0; f0 = f(x0)
    for it in range(maxit):
        a1 = c1*(x1-x0-h*v0) - (0.5/b-1)*a0
        v1 = v0 + h*((1-g)*a0 + g*a1)
        R = (1-am)*a1 + am*a0 + d*((1-af)*v1 + af*v0) + (1-af)*f(x1) + af*f0
        if force_iters is None and np.all(np.abs(R) <= tol): break
        dR = (1-am)*c1 + d*(1-af)*g*h*c1 + (1-af)*k(x1)
        x1 = x1 - R/dR
        if force_iters is not None and it+1 >= force_iters: break
    a1 = c1*(x1-x0-h*v0) - (0.5/b-1)*a0
    v1 = v0 + h*((1-g)*a0 + g*a1)
    return x1, v1, a1

def spec_rad_inf(p):
    f = lambda x: x; k = lambda x: np.ones_like(x)
    M = np.zeros((3, 3))
    for j in range(3):
        s = np.zeros(3); s[j] = 1.0
        x, v, a = step(np.array(s[0]), np.array(s[1]), np.array(s[2]), 1e3, f, k, 0.0, p, force_iters=2)
        M[:, j] = [float(x), float(v), float(a)]
    return np.max(np.abs(np.linalg.eigvals(M)))

def fam(name, s):
    if name == "genalpha":
        r = s; am = (2*r-1)/(r+1); af = r/(r+1); return (am, af, 0.5-am+af, 0.25*(1-am+af)**2)
    if name == "HHT":      # s = a in [0,1/3]
        return (0.0, s, 0.5+s, 0.25*(1+s)**2)
    if name == "WBZ":      # s = m in [-1,0]
        return (s, 0.0, 0.5-s, 0.25*(1-s)**2)
    if name == "Newmark":  # s = delta >= 0
        return (0.0, 0.0, 0.5+s, 0.25*(1+s)**2)

RANGE = {"HHT": (0.0, 1/3), "WBZ": (0.0, -1.0), "Newmark": (0.0, 3.0)}
def param_for_rho(name, target):
    if name == "genalpha": return fam(name, target), target
    lo, hi = RANGE[name]; rlo, rhi = spec_rad_inf(fam(name, lo)), spec_rad_inf(fam(name, hi))
    if not (min(rlo, rhi)-1e-9 <= target <= max(rlo, rhi)+1e-9): return None, None
    for _ in range(60):
        mid = 0.5*(lo+hi); rm = spec_rad_inf(fam(name, mid))
        if (rm - target)*(rlo - target) > 0: lo, rlo = mid, rm
        else: hi = mid
    s = 0.5*(lo+hi); return fam(name, s), spec_rad_inf(fam(name, s))

def make_system(e):
    return (lambda x: x**3 - x - e), (lambda x: 3*x**2 - 1)

def ref(X, V, tend, DC, e, hr=0.01):
    x, v = X.copy(), V.copy(); rhs = lambda x, v: (v, -DC*v + x - x**3 + e)
    for _ in range(int(round(tend/hr))):
        k1 = rhs(x, v); k2 = rhs(x+hr/2*k1[0], v+hr/2*k1[1]); k3 = rhs(x+hr/2*k2[0], v+hr/2*k2[1]); k4 = rhs(x+hr*k3[0], v+hr*k3[1])
        x = x+hr/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]); v = v+hr/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
    return np.sign(x)

def run(X, V, tend, h, p, DC, e):
    f, k = make_system(e); x, v = X.copy(), V.copy(); a = -DC*v - f(x)
    with np.errstate(all="ignore"):
        for _ in range(int(round(tend/h))):
            x, v, a = step(x, v, a, h, f, k, DC, p)
    o = np.sign(x); o[~np.isfinite(x)] = 0; return o

if __name__ == "__main__":
    n, tend = 100, 60.0; xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs)
    targets = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    for (DC, h, e) in [(0.15, 0.4, 0.0), (0.05, 0.2, 0.0), (0.15, 0.4, 0.1)]:
        R = ref(X, V, tend, DC, e)
        print(f"\n=== c={DC} h={h} tilt e={e}  basin mismatch % vs rho_inf (family-matched) ===")
        print("rho_inf     " + " ".join(f"{t:5.1f}" for t in targets))
        for name in ["genalpha", "HHT", "WBZ", "Newmark"]:
            row = []
            for t in targets:
                p, r = param_for_rho(name, t)
                row.append("  n/a" if p is None else f"{100*np.mean(run(X, V, tend, h, p, DC, e) != R):5.1f}")
            print(f"{name:<11} " + " ".join(row))
