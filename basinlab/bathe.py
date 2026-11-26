"""rho_inf-Bathe composite integrator (Noh & Bathe 2019, Comput. Struct. 212,
299-310; DOI 10.1016/j.compstruc.2018.11.001), nonlinear/MDOF, batched, with
per-IC failure masking.  M a + C v + f(q) = 0.

Sub-step 1 (t -> t+gamma h): trapezoidal rule.
Sub-step 2 (t -> t+h):
  q1 = q0 + h (c0 v0 + c1 v_g + c2 v1),   v1 = v0 + h (d0 a0 + d1 a_g + d2 a1)
with the paper's relations (Eq. 14, 18, 19, 21), written with the paper's q's
renamed c_i and s's renamed d_i:
  c1 = d1 = (rho + 1) / (2 gamma (rho - 1) + 4)
  c0 = d0 = (gamma - 1) c1 + 1/2,   c2 = d2 = 1/2 - gamma c1      (sum = 1)
  gamma = gamma0(rho) = (2 - sqrt(2 + 2 rho)) / (1 - rho)  (0.5 at rho = 1)
Checks (tests): rho=1 gives c = (1/4, 1/2, 1/4) = two trapezoidal sub-steps;
rho=0 gives gamma = 2 - sqrt 2; second-order accuracy; spectral radius at
large h/T equals rho.  The coefficient c2 = 1/2 - gamma c1 follows from the
consistency condition c0 + c1 + c2 = 1 and the standard Bathe method
(rho=0, gamma=1/2 gives c = 1/3 each)."""
import numpy as np


def gamma0(rho):
    if abs(rho - 1.0) < 1e-12:
        return 0.5
    return (2.0 - np.sqrt(2.0 + 2.0*rho))/(1.0 - rho)


def coefficients(rho, gamma=None):
    g = gamma0(rho) if gamma is None else gamma
    c1 = (rho + 1.0)/(2.0*g*(rho - 1.0) + 4.0)
    c0 = (g - 1.0)*c1 + 0.5
    c2 = 0.5 - g*c1
    return g, (c0, c1, c2)


def step_masked(q, v, a, h, M, C, model, rho, tol=1e-10, maxiter=60):
    g, (c0, c1, c2) = coefficients(rho)
    nb = q.shape[0]
    f0, _ = model.force_tangent(q)
    scale = np.maximum(1.0, np.max(np.abs(f0), axis=-1))
    # ---- sub-step 1: trapezoidal to t + g h --------------------------------
    hg = g*h
    qg = q + hg*v + 0.5*hg*hg*a
    active = np.ones(nb, bool)
    for _ in range(maxiter + 1):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        vg = 2*(qg[idx] - q[idx])/hg - v[idx]
        ag = 2*(vg - v[idx])/hg - a[idx]
        f, K = model.force_tangent(qg[idx])
        r = ag @ M.T + vg @ C.T + f
        ok = (np.max(np.abs(r), axis=-1) <= tol*scale[idx]) | ~np.all(np.isfinite(r), axis=-1)
        active[idx[ok]] = False
        work = ~ok
        if not np.any(work):
            continue
        J = 4*M/hg**2 + 2*C/hg + K[work]
        qg[idx[work]] -= np.linalg.solve(J, r[work][..., None])[..., 0]
    failed1 = active | ~np.all(np.isfinite(qg), axis=-1)
    vg = 2*(qg - q)/hg - v
    ag = 2*(vg - v)/hg - a
    # ---- sub-step 2 ----------------------------------------------------------
    q1 = q + h*v + 0.5*h*h*a
    active = ~failed1
    for _ in range(maxiter + 1):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        v1 = (q1[idx] - q[idx])/(h*c2) - (c0*v[idx] + c1*vg[idx])/c2
        a1 = (v1 - v[idx])/(h*c2) - (c0*a[idx] + c1*ag[idx])/c2
        f, K = model.force_tangent(q1[idx])
        r = a1 @ M.T + v1 @ C.T + f
        ok = (np.max(np.abs(r), axis=-1) <= tol*scale[idx]) | ~np.all(np.isfinite(r), axis=-1)
        active[idx[ok]] = False
        work = ~ok
        if not np.any(work):
            continue
        J = M/(h*h*c2*c2) + C/(h*c2) + K[work]
        q1[idx[work]] -= np.linalg.solve(J, r[work][..., None])[..., 0]
    failed = failed1 | active | ~np.all(np.isfinite(q1), axis=-1)
    v1 = (q1 - q)/(h*c2) - (c0*v + c1*vg)/c2
    a1 = (v1 - v)/(h*c2) - (c0*a + c1*ag)/c2
    return q1, v1, a1, failed


def integrate_to_masked(q, v, end, h, M, C, model, rho, maxiter=60):
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    f, _ = model.force_tangent(q)
    a = np.linalg.solve(M, (-(v @ C.T) - f)[..., None])[..., 0]
    alive = np.ones(q.shape[0], bool)
    nfull = int(np.floor(end/h + 1e-12))
    steps = [h]*nfull + ([end - nfull*h] if end - nfull*h > 1e-10 else [])
    for hs in steps:
        i = np.flatnonzero(alive)
        if i.size == 0:
            break
        q1, v1, a1, failed = step_masked(q[i], v[i], a[i], hs, M, C, model, rho, maxiter=maxiter)
        q[i], v[i], a[i] = q1, v1, a1
        alive[i[failed]] = False
    q[~alive] = np.nan
    v[~alive] = np.nan
    return q, v, int((~alive).sum())
