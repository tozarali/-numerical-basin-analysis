"""Generalized-alpha for M a + C v + f(q) = 0, batched over ICs.
model.force_tangent(q) -> (f [...,n], K [...,n,n]); M, C constant (n,n).
Used for truss, Beam2Mode and the FE arch (reduce batch size for FE)."""
import numpy as np


def step(q, v, a, h, M, C, model, p, tol=1e-10, maxiter=30):
    am, af, gamma, beta = p
    c1 = 1/(beta*h*h)
    f0, _ = model.force_tangent(q)
    q1 = q + h*v + 0.5*h*h*a
    for it in range(maxiter + 1):
        a1 = c1*(q1 - q - h*v) - (0.5/beta - 1)*a
        v1 = v + h*((1-gamma)*a + gamma*a1)
        f1, K1 = model.force_tangent(q1)
        r = ((1-am)*a1 + am*a) @ M.T + ((1-af)*v1 + af*v) @ C.T + (1-af)*f1 + af*f0
        scale = max(1.0, float(np.max(np.abs(f0))))
        if np.max(np.abs(r)) <= tol*scale:
            return q1, v1, a1
        if it == maxiter:
            break
        J = (1-am)*c1*M + (1-af)*gamma*h*c1*C + (1-af)*K1
        q1 = q1 - np.linalg.solve(J, r[..., None])[..., 0]
    raise RuntimeError("mdof generalized-alpha Newton failed")


def integrate_to(q, v, end, h, M, C, model, p):
    """Exact end time via one shortened final step (as energy_loss_retest)."""
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    f, _ = model.force_tangent(q)
    a = np.linalg.solve(M, (-(v @ C.T) - f)[..., None])[..., 0]
    n = int(np.floor(end/h + 1e-12))
    for _ in range(n):
        q, v, a = step(q, v, a, h, M, C, model, p)
    rem = end - n*h
    if rem > 1e-10:
        q, v, a = step(q, v, a, rem, M, C, model, p)
    return q, v


def step_masked(q, v, a, h, M, C, model, p, tol=1e-10, maxiter=60):
    """Per-IC Newton: only unconverged rows are updated; rows that do not
    converge within maxiter are returned in `failed` (caller marks them
    unresolved instead of aborting the whole batch)."""
    am, af, gamma, beta = p
    c1 = 1/(beta*h*h)
    f0, _ = model.force_tangent(q)
    scale = np.maximum(1.0, np.max(np.abs(f0), axis=-1))
    q1 = q + h*v + 0.5*h*h*a
    active = np.ones(q.shape[0], bool)
    J0 = (1-am)*c1*M + (1-af)*gamma*h*c1*C
    for _ in range(maxiter + 1):
        a1 = c1*(q1 - q - h*v) - (0.5/beta - 1)*a
        v1 = v + h*((1-gamma)*a + gamma*a1)
        idx = np.flatnonzero(active)
        if idx.size == 0:
            break
        f1, K1 = model.force_tangent(q1[idx])
        r = (((1-am)*a1[idx] + am*a[idx]) @ M.T + ((1-af)*v1[idx] + af*v[idx]) @ C.T
             + (1-af)*f1 + af*f0[idx])
        ok = np.max(np.abs(r), axis=-1) <= tol*scale[idx]
        ok |= ~np.all(np.isfinite(r), axis=-1)       # give up on non-finite rows
        active[idx[ok]] = False
        work = ~ok
        if not np.any(work):
            continue
        J = J0 + (1-af)*K1[work]
        q1[idx[work]] -= np.linalg.solve(J, r[work][..., None])[..., 0]
    a1 = c1*(q1 - q - h*v) - (0.5/beta - 1)*a
    v1 = v + h*((1-gamma)*a + gamma*a1)
    failed = active | ~np.all(np.isfinite(q1), axis=-1)
    return q1, v1, a1, failed


def integrate_to_masked(q, v, end, h, M, C, model, p, maxiter=60):
    """Like integrate_to, but failed ICs are frozen as NaN (-> label 0)."""
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    f, _ = model.force_tangent(q)
    a = np.linalg.solve(M, (-(v @ C.T) - f)[..., None])[..., 0]
    alive = np.ones(q.shape[0], bool)
    n = int(np.floor(end/h + 1e-12))
    steps = [h]*n + ([end - n*h] if end - n*h > 1e-10 else [])
    for hs in steps:
        i = np.flatnonzero(alive)
        if i.size == 0:
            break
        q1, v1, a1, failed = step_masked(q[i], v[i], a[i], hs, M, C, model, p, maxiter=maxiter)
        q[i], v[i], a[i] = q1, v1, a1
        alive[i[failed]] = False
    q[~alive] = np.nan
    v[~alive] = np.nan
    return q, v, int((~alive).sum())
