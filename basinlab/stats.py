"""Continuous optimum location and bootstrap uncertainty over ICs."""
import numpy as np


def parabolic_minimum(rhos, values, half_width=3):
    """Fit a parabola to the points around the discrete minimum.
    Returns the vertex if it is a minimum inside the fitted window, else the
    discrete argmin (flagged)."""
    rhos, values = np.asarray(rhos, float), np.asarray(values, float)
    i = int(np.nanargmin(values))
    lo, hi = max(0, i-half_width), min(len(rhos), i+half_width+1)
    if hi - lo < 3:
        return float(rhos[i]), "edge"
    a, b, _ = np.polyfit(rhos[lo:hi], values[lo:hi], 2)
    if a <= 0:
        return float(rhos[i]), "not_convex"
    vertex = -b/(2*a)
    if not rhos[lo] <= vertex <= rhos[hi-1]:
        return float(rhos[i]), "outside"
    return float(vertex), "ok"


def bootstrap_optimum(rhos, errors, n_boot=1000, seed=0, half_width=3):
    """errors: bool array (n_rho, n_ic), True = label differs from reference
    for that IC. ICs are resampled jointly across rho (they are shared),
    so the CI reflects IC-sampling uncertainty of the optimum location."""
    errors = np.asarray(errors, bool)
    rng = np.random.default_rng(seed)
    n_ic = errors.shape[1]
    point, flag = parabolic_minimum(rhos, errors.mean(axis=1), half_width)
    samples = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n_ic, n_ic)
        samples[b] = parabolic_minimum(rhos, errors[:, idx].mean(axis=1), half_width)[0]
    lo, hi = np.percentile(samples, [2.5, 97.5])
    return {"rho_star": point, "flag": flag, "ci95": [float(lo), float(hi)],
            "boot_sd": float(samples.std(ddof=1))}


def proportion_ci(k, n, z=1.96):
    """Wilson interval for a mismatch fraction."""
    if n == 0:
        return (np.nan, np.nan)
    p = k/n
    den = 1 + z*z/n
    centre = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return (float(centre-half), float(centre+half))
