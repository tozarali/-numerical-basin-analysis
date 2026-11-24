"""Exploratory energy-loss proxies and averaged linear theory.

The legacy ``d_sep`` name measures a fixed-time loss, not a verified orbit
or a horizon-free quantity. See experiments.dsep_exploratory for event checks.
"""
import numpy as np
from scipy.integrate import solve_ivp

from .integrators import step
from .model import energy, force
from .linear_analysis import principal_pair

BARRIER = 0.25  # double well V=(x^2-1)^2/4, saddle at x=0


def _gen_alpha_to(x, v, end, h, damping, p):
    """Scalar/array gen-alpha with an exact shortened last step."""
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    a = -damping*v - force(x)
    n = int(np.floor(end/h + 1e-12))
    for _ in range(n):
        x, v, a = step(x, v, a, h, damping, p)
    rem = end - n*h
    if rem > 1e-10:
        x, v, a = step(x, v, a, rem, damping, p)
    return x, v


def _reference_to(x, v, end, damping, rtol=1e-11, atol=1e-13):
    sol = solve_ivp(lambda t, y: (y[1], -damping*y[1] - force(y[0])),
                    (0, end), (x, v), method="DOP853", rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0, -1], sol.y[1, -1]


def conservative_period(E):
    """Period of the undamped orbit with energy E (double well)."""
    if E <= 0 or abs(E - BARRIER) < 1e-9:
        raise ValueError("E must be positive and not on the separatrix")
    if E < BARRIER:
        x0, v0 = 1.0, np.sqrt(2*E)          # well bottom, moving right
    else:
        x0, v0 = 0.0, np.sqrt(2*(E - BARRIER))  # saddle, moving right
    ev = lambda t, y: y[0] - x0
    ev.direction = 1
    sol = solve_ivp(lambda t, y: (y[1], -force(y[0])), (0, 1e4), (x0 - 1e-12, v0),
                    method="DOP853", rtol=1e-12, atol=1e-14, events=ev)
    times = sol.t_events[0]
    times = times[times > 1e-6]
    if not len(times):
        raise RuntimeError("period not found")
    return float(times[0])


def d_sep(p, h, damping, delta=0.05, window="half_period"):
    """Fixed-time energy-loss proxy starting at the saddle, not a true cycle.

    The historical D_sep return key is retained for compatibility. The default
    half-period is that of an undamped orbit and does not imply a damped orbit
    actually completes a crossing or excursion.

    Starting at the saddle
    with E0 = BARRIER*(1+delta). Window = half the conservative period of E0
    matched exactly in time for both solvers."""
    E0 = BARRIER * (1 + delta)
    x0, v0 = 0.0, np.sqrt(2*(E0 - BARRIER))
    tau = conservative_period(E0) / 2 if window == "half_period" else float(window)
    xr, vr = _reference_to(x0, v0, tau, damping)
    xn, vn = _gen_alpha_to(x0, v0, tau, h, damping, p)
    d_ref = E0 - energy(xr, vr)
    d_num = E0 - energy(xn, vn)
    return {"D_sep": float(d_num/d_ref), "E0": E0, "window": tau,
            "dE_ref": float(d_ref), "dE_num": float(d_num)}


def d_theory(p, h, damping, E):
    """Level-0 averaged theory: treat the orbit of energy E as harmonic with
    frequency w(E)=2*pi/T(E); ratio of numerical to physical energy loss per
    cycle from the exact linear amplification at (w h, xi=c/(2w)).
    KNOWN LIMITATION: near the separatrix the orbit is far from harmonic;
    level-1 (TODO) weights local losses by c*v^2 along the orbit."""
    w = 2*np.pi / conservative_period(E)
    xi = damping / (2*w)
    lam = principal_pair(p, w*h, xi)
    if lam is None:
        return np.nan
    steps_per_cycle = 2*np.pi / (w*h)
    retained_num = abs(lam) ** (2*steps_per_cycle)
    retained_phys = np.exp(-2*xi*w * 2*np.pi/w)   # = exp(-2*pi*c/w)
    return float((1 - retained_num) / (1 - retained_phys))


def crossing(rhos, values, target=1.0):
    """Linear interpolation of the first crossing of `target` (None if none)."""
    rhos, values = np.asarray(rhos, float), np.asarray(values, float)
    idx = np.flatnonzero((values[:-1]-target)*(values[1:]-target) <= 0)
    if not len(idx):
        return None
    i = int(idx[0])
    if values[i+1] == values[i]:
        return float((rhos[i]+rhos[i+1])/2)
    return float(rhos[i] + (target-values[i])*(rhos[i+1]-rhos[i])/(values[i+1]-values[i]))
