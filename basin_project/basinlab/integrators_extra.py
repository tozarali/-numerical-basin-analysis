"""Extra SDOF schemes for the double well (vectorised over ICs).
Same signature idea as integrators.step: state arrays in, state arrays out.
Each scheme is a 'point' (no rho parameter) or has its own parameter."""
import numpy as np

from .model import force, tangent


def potential(x):
    return 0.25*(x*x - 1)**2


def discrete_gradient_step(x, v, h, damping, tol=1e-12, maxiter=30):
    """Gonzalez/Itoh-Abe discrete gradient (energy-consistent control):
    x1-x0 = h(v0+v1)/2,  v1-v0 = -h[c(v0+v1)/2 + (V(x1)-V(x0))/(x1-x0)].
    For the double well V=(x^2-1)^2/4 the difference quotient is exact and
    non-singular: dg = (x1+x0)(x1^2+x0^2-2)/4. With c=0 energy is conserved
    to solver tolerance."""
    x = np.asarray(x, float)
    v = np.asarray(v, float)
    x1 = x + h*v
    for _ in range(maxiter):
        dg = 0.25*(x1 + x)*(x1*x1 + x*x - 2)
        v1 = 2*(x1 - x)/h - v
        r = v1 - v + h*(damping*0.5*(v + v1) + dg)
        scale = 1 + np.abs(v) + h*np.abs(dg)
        if np.all(np.abs(r) <= tol*scale):
            return x1, v1
        ddg = 0.25*((x1*x1 + x*x - 2) + 2*x1*(x1 + x))
        x1 = x1 - r/(2/h + damping + h*ddg)
    raise RuntimeError("discrete gradient Newton failed")


def bdf2_integrate(x, v, end, h, damping, tol=1e-12, maxiter=30):
    """BDF2 on (x, v) with one trapezoidal start step; rho_inf = 0, O(h^2).
    Exact end time: requires end/h integer (keeps BDF2 coefficients fixed)."""
    n = round(end/h)
    if not np.isclose(n*h, end, atol=1e-10):
        raise ValueError("end must be a multiple of h for BDF2")
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    # start: trapezoidal (Newmark 1/4,1/2) single step via discrete Newton
    xm, vm = x.copy(), v.copy()
    x1 = x + h*v
    for _ in range(maxiter):
        v1 = 2*(x1 - x)/h - v
        r = (v1 - v) + h*0.5*(damping*(v+v1) + force(x) + force(x1))
        if np.all(np.abs(r) < tol):
            break
        x1 -= r / (2/h + h*0.5*(damping*2/h + tangent(x1)))
    x, v = x1, 2*(x1 - xm)/h - vm
    for _ in range(n - 1):
        xg = 2*x - xm
        for _ in range(maxiter):
            v1 = (3*xg - 4*x + xm) / (2*h)
            r = 3*v1 - 4*v + vm + 2*h*(damping*v1 + force(xg))
            if np.all(np.abs(r) < tol):
                break
            drdx = 9/(2*h) + 2*h*(damping*3/(2*h) + tangent(xg))
            xg = xg - r/drdx
        else:
            raise RuntimeError("BDF2 Newton failed")
        xm, vm, x, v = x, v, xg, (3*xg - 4*x + xm)/(2*h)
    return x, v


def bathe_rho_integrate(x, v, end, h, damping, rho_inf):
    """rho_inf-Bathe composite scheme (Noh & Bathe 2019).
    TODO: implement from the original paper; verify rho_inf by
    linear_analysis-style amplification of a linear SDOF before use.
    Do NOT write coefficients from memory."""
    raise NotImplementedError("Implement from Noh & Bathe (2019) and test first")
