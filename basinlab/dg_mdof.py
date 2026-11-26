"""Energy-consistent (Gonzalez midpoint discrete gradient) integrator for
M q'' + C q' + grad V(q) = 0, batched over ICs, with per-IC failure masking.

  q1 - q0 = h vbar,  vbar = (v0 + v1)/2
  M (v1 - v0) = -h [ C vbar + G(q0, q1) ]
  G = g + (s/|dq|^2) dq,  g = grad V(qbar),  s = V(q1) - V(q0) - g.dq
so that G.dq = V(q1) - V(q0) exactly: with C = 0 the total energy is conserved
for any V (Gonzalez 1996). Needs model.potential and model.force_tangent.
In one dimension G is the exact difference quotient (V1-V0)/(q1-q0)."""
import numpy as np

_TINY = 1e-22


def _G_and_jac(model, q0, q1):
    qb = 0.5*(q0 + q1)
    dq = q1 - q0
    gb, Kb = model.force_tangent(qb)
    f1, _ = model.force_tangent(q1)
    s = model.potential(q1) - model.potential(q0) - np.sum(gb*dq, axis=-1)
    d2 = np.sum(dq*dq, axis=-1)
    small = d2 < _TINY
    d2s = np.where(small, 1.0, d2)
    G = gb + (s/d2s)[:, None]*dq
    # d s / d q1 = f1 - gb - Kb dq / 2
    ds = f1 - gb - 0.5*np.einsum("bij,bj->bi", Kb, dq)
    n = q0.shape[-1]
    J = (0.5*Kb + np.einsum("bi,bj->bij", dq, ds)/d2s[:, None, None]
         + (s/d2s)[:, None, None]*np.eye(n)
         - 2*(s/d2s**2)[:, None, None]*np.einsum("bi,bj->bij", dq, dq))
    G = np.where(small[:, None], gb, G)
    J = np.where(small[:, None, None], 0.5*Kb, J)
    return G, J


def step_masked(q, v, h, M, C, model, tol=1e-10, maxiter=60):
    q, v = np.asarray(q, float), np.asarray(v, float)
    nb, n = q.shape
    scale = np.maximum(1.0, np.max(np.abs(model.force_tangent(q)[0]), axis=-1))
    q1 = q + h*v
    active = np.ones(nb, bool)
    A = 2*M/h**2 + C/h
    for _ in range(maxiter + 1):
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        qi, q1i = q[idx], q1[idx]
        vb = (q1i - qi)/h
        G, dG = _G_and_jac(model, qi, q1i)
        r = (2*(vb - v[idx])/h) @ M.T + vb @ C.T + G
        ok = np.max(np.abs(r), axis=-1) <= tol*scale[idx]
        ok |= ~np.all(np.isfinite(r), axis=-1)
        active[idx[ok]] = False
        work = ~ok
        if not np.any(work):
            continue
        q1[idx[work]] -= np.linalg.solve(A + dG[work], r[work][..., None])[..., 0]
    v1 = 2*(q1 - q)/h - v
    failed = active | ~np.all(np.isfinite(q1), axis=-1)
    return q1, v1, failed


def integrate_to_masked(q, v, end, h, M, C, model, maxiter=60):
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    alive = np.ones(q.shape[0], bool)
    nfull = int(np.floor(end/h + 1e-12))
    steps = [h]*nfull + ([end - nfull*h] if end - nfull*h > 1e-10 else [])
    for hs in steps:
        i = np.flatnonzero(alive)
        if i.size == 0:
            break
        q1, v1, failed = step_masked(q[i], v[i], hs, M, C, model, maxiter=maxiter)
        q[i], v[i] = q1, v1
        alive[i[failed]] = False
    q[~alive] = np.nan
    v[~alive] = np.nan
    return q, v, int((~alive).sum())
