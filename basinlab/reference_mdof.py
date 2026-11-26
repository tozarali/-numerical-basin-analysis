"""Independent explicit RK4 reference for M q'' + C q' + f(q) = 0 with a
diagonal (lumped) mass, batched over ICs. Uses model.force if present
(no tangent needed), otherwise model.force_tangent(q)[0]."""
import numpy as np


def _force(model, q):
    if hasattr(model, "force"):
        return model.force(q)
    return model.force_tangent(q)[0]


def rk4(q, v, end, h, m_diag, C, model):
    n = round(end/h)
    if n < 1 or not np.isclose(n*h, end, rtol=0, atol=1e-9):
        raise ValueError("end must be a positive multiple of h")
    q, v = np.asarray(q, float).copy(), np.asarray(v, float).copy()
    minv = 1.0/np.asarray(m_diag, float)

    def acc(qq, vv):
        return -(vv @ C.T + _force(model, qq))*minv

    for _ in range(n):
        k1q, k1v = v, acc(q, v)
        k2q, k2v = v + h*k1v/2, acc(q + h*k1q/2, v + h*k1v/2)
        k3q, k3v = v + h*k2v/2, acc(q + h*k2q/2, v + h*k2v/2)
        k4q, k4v = v + h*k3v, acc(q + h*k3q, v + h*k3v)
        q = q + h*(k1q + 2*k2q + 2*k3q + k4q)/6
        v = v + h*(k1v + 2*k2v + 2*k3v + k4v)/6
    if not np.all(np.isfinite(q)):
        raise FloatingPointError("non-finite RK4 reference state")
    return q, v
