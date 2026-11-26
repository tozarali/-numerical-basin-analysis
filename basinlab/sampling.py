"""Pre-registrable IC populations (fixed seed, recorded in summary.json)."""
import numpy as np
from scipy.stats import qmc


def sobol_box(n, lower, upper, seed=12345):
    """Scrambled Sobol points in a box; n is rounded up to a power of 2."""
    m = int(np.ceil(np.log2(n)))
    pts = qmc.Sobol(d=len(lower), scramble=True, seed=seed).random_base2(m)
    return qmc.scale(pts, lower, upper)


def energy_band(points, energy_fn, e_low, e_high):
    """Keep points with e_low < E0 <= e_high (e.g. (E_b, k*E_b])."""
    e = energy_fn(*points.T)
    keep = (e > e_low) & (e <= e_high)
    return points[keep], e[keep]
