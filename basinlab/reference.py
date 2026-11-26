"""Independent fixed-step RK4 and adaptive SciPy reference trajectories."""
import numpy as np
from scipy.integrate import solve_ivp

from .model import force


def rk4(x, v, end, h, damping):
    steps = round(end / h)
    if damping <= 0 or h <= 0 or steps < 1 or not np.isclose(steps * h, end, rtol=0, atol=1e-10):
        raise ValueError("damping and h must be positive and end must be a positive multiple of h")
    x, v = np.asarray(x, dtype=float).copy(), np.asarray(v, dtype=float).copy()

    def rhs(q, w):
        return w, -damping * w - force(q)

    for _ in range(steps):
        q1, w1 = rhs(x, v)
        q2, w2 = rhs(x + h*q1/2, v + h*w1/2)
        q3, w3 = rhs(x + h*q2/2, v + h*w2/2)
        q4, w4 = rhs(x + h*q3, v + h*w3)
        x = x + h*(q1 + 2*q2 + 2*q3 + q4)/6
        v = v + h*(w1 + 2*w2 + 2*w3 + w4)/6
    if not np.all(np.isfinite(x) & np.isfinite(v)):
        raise FloatingPointError("Non-finite RK4 reference state")
    return x, v


def adaptive(x, v, end, damping, method="DOP853", rtol=1e-11, atol=1e-13):
    if end <= 0 or damping <= 0 or rtol <= 0 or atol <= 0:
        raise ValueError("end, damping and tolerances must be positive")
    x, v = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(v, dtype=float))
    out_x, out_v = np.empty(x.shape), np.empty(v.shape)
    for index in np.ndindex(x.shape):
        solution = solve_ivp(lambda t, y: (y[1], -damping*y[1] - force(y[0])),
                             (0, end), (x[index], v[index]), method=method,
                             rtol=rtol, atol=atol)
        if not solution.success:
            raise RuntimeError(f"{method} failed at {index}: {solution.message}")
        out_x[index], out_v[index] = solution.y[:, -1]
    return out_x, out_v
