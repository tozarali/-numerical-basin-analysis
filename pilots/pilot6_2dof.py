"""2-DOF test: bistable mode x coupled to a stiff high-frequency mode y (the mode engineers damp with rho_inf).
V = (x^2-1)^2/4 + w^2 y^2/2 + g x^2 y^2/2. Wells at (+-1, 0). Basin slice: (x0, vx0), y0=Y0, vy0=0."""
import numpy as np
W, G, C1, C2, Y0 = 3.0, 1.0, 0.15, 0.30, 0.5

def fint(x, y): return x**3 - x + G*x*y*y, W*W*y + G*x*x*y
def hess(x, y): return 3*x*x - 1 + G*y*y, 2*G*x*y, W*W + G*x*x     # Kxx, Kxy, Kyy

def ref(X, V, tend, hr=0.005):
    x, y, vx, vy = X.copy(), np.full_like(X, Y0), V.copy(), np.zeros_like(X)
    def rhs(x, y, vx, vy):
        fx, fy = fint(x, y); return vx, vy, -C1*vx - fx, -C2*vy - fy
    for _ in range(int(round(tend/hr))):
        k1 = rhs(x, y, vx, vy)
        k2 = rhs(*(s + hr/2*q for s, q in zip((x, y, vx, vy), k1)))
        k3 = rhs(*(s + hr/2*q for s, q in zip((x, y, vx, vy), k2)))
        k4 = rhs(*(s + hr*q for s, q in zip((x, y, vx, vy), k3)))
        x, y, vx, vy = (s + hr/6*(a + 2*b + 2*c + d) for s, a, b, c, d in zip((x, y, vx, vy), k1, k2, k3, k4))
    return np.sign(x)

def ga(X, V, tend, h, rho, tol=1e-10, maxit=30):
    am = (2*rho-1)/(rho+1); af = rho/(rho+1); g = 0.5-am+af; b = 0.25*(1-am+af)**2; c1 = 1/(b*h*h)
    x, y, vx, vy = X.copy(), np.full_like(X, Y0), V.copy(), np.zeros_like(X)
    fx, fy = fint(x, y); ax, ay = -C1*vx - fx, -C2*vy - fy
    with np.errstate(all="ignore"):
        for _ in range(int(round(tend/h))):
            f0x, f0y = fx, fy
            x1, y1 = x + h*vx + 0.5*h*h*ax, y + h*vy + 0.5*h*h*ay
            for _ in range(maxit):
                a1x = c1*(x1-x-h*vx) - (0.5/b-1)*ax; a1y = c1*(y1-y-h*vy) - (0.5/b-1)*ay
                v1x = vx + h*((1-g)*ax + g*a1x); v1y = vy + h*((1-g)*ay + g*a1y)
                f1x, f1y = fint(x1, y1)
                Rx = (1-am)*a1x + am*ax + C1*((1-af)*v1x + af*vx) + (1-af)*f1x + af*f0x
                Ry = (1-am)*a1y + am*ay + C2*((1-af)*v1y + af*vy) + (1-af)*f1y + af*f0y
                if max(np.nanmax(np.abs(Rx)), np.nanmax(np.abs(Ry))) <= tol: break
                kxx, kxy, kyy = hess(x1, y1)
                Axx = (1-am)*c1 + C1*(1-af)*g*h*c1 + (1-af)*kxx; Ayy = (1-am)*c1 + C2*(1-af)*g*h*c1 + (1-af)*kyy; Axy = (1-af)*kxy
                det = Axx*Ayy - Axy*Axy
                x1 = x1 - (Ayy*Rx - Axy*Ry)/det; y1 = y1 - (-Axy*Rx + Axx*Ry)/det
            a1x = c1*(x1-x-h*vx) - (0.5/b-1)*ax; a1y = c1*(y1-y-h*vy) - (0.5/b-1)*ay
            v1x = vx + h*((1-g)*ax + g*a1x); v1y = vy + h*((1-g)*ay + g*a1y)
            fx, fy = fint(x1, y1)
            x, y, vx, vy, ax, ay = x1, y1, v1x, v1y, a1x, a1y
    o = np.sign(x); o[~np.isfinite(x)] = 0; return o

if __name__ == "__main__":
    n, tend = 100, 60.0; xs = np.linspace(-2, 2, n); X, V = np.meshgrid(xs, xs)
    R = ref(X, V, tend)
    rhos = [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
    print("2-DOF basin mismatch (percent) vs RK4 reference (slice x0,vx0; y0=%.1f)" % Y0)
    print("h \\ rho   " + " ".join(f"{r:5.1f}" for r in rhos))
    for h in [0.1, 0.2, 0.3, 0.4]:
        print(f"h={h:<6}  " + " ".join(f"{100*np.mean(ga(X.copy(), V.copy(), tend, h, r) != R):5.1f}" for r in rhos))
