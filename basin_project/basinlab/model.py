"""Damped, unforced symmetric double-well oscillator."""
import numpy as np


def force(x):
    return x**3 - x


def tangent(x):
    return 3 * x**2 - 1


def energy(x, v):
    return 0.5 * v**2 + 0.25 * (x**2 - 1)**2


def classify(x, v):
    """Return +/-1 only when physical energy is strictly below the barrier."""
    x, v = np.broadcast_arrays(np.asarray(x), np.asarray(v))
    result = np.zeros(x.shape, dtype=np.int8)
    captured = np.isfinite(x) & np.isfinite(v) & (energy(x, v) < 0.25)
    result[captured & (x > 0)] = 1
    result[captured & (x < 0)] = -1
    return result
