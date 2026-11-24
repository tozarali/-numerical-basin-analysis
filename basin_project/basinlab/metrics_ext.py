"""Additional, saturation-resistant basin measures (opus.md WP1-3)."""
import numpy as np

from .metrics import boundary_points, entropy


def net_basin_transfer(reference, candidate):
    """Signed area measure: fraction of resolved ICs moved -1 -> +1 minus
    +1 -> -1. Zero for symmetric errors; sign shows which basin grows."""
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    valid = (reference != 0) & (candidate != 0)
    n = int(valid.sum())
    if n == 0:
        return {"net": None, "to_plus": None, "to_minus": None, "valid": 0}
    to_plus = np.mean((reference == -1) & (candidate == 1) & valid) * valid.size / n
    to_minus = np.mean((reference == 1) & (candidate == -1) & valid) * valid.size / n
    return {"net": float(to_plus - to_minus), "to_plus": float(to_plus),
            "to_minus": float(to_minus), "valid": n}


def signed_boundary_shift(reference, candidate, spacing):
    """Signed displacement of the reference basin boundary.
    For each pair of adjacent reference cells with labels -1/+1, look at the
    candidate labels of the SAME two cells: both +1 -> the + basin grew there
    (sign +1); both -1 -> the - basin grew (sign -1); otherwise 0. The
    magnitude is the distance from the pair midpoint to the nearest candidate
    boundary midpoint. Returns mean signed and mean absolute shift (same
    units as `spacing`) and the fractions of growing/shrinking pairs."""
    from scipy.spatial import cKDTree
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    if np.any(reference == 0) or np.any(candidate == 0):
        raise ValueError("requires fully resolved grids")
    dy, dx = spacing
    pts, signs = [], []
    for axis in (0, 1):
        a = [slice(None), slice(None)]; b = [slice(None), slice(None)]
        a[axis], b[axis] = slice(None, -1), slice(1, None)
        ra, rb = reference[tuple(a)], reference[tuple(b)]
        ca, cb = candidate[tuple(a)], candidate[tuple(b)]
        rows, cols = np.nonzero(ra*rb == -1)
        sgn = np.where((ca[rows, cols] == 1) & (cb[rows, cols] == 1), 1,
                       np.where((ca[rows, cols] == -1) & (cb[rows, cols] == -1), -1, 0))
        if axis == 0:
            pts.append(np.column_stack(((rows + 0.5)*dy, cols*dx)))
        else:
            pts.append(np.column_stack((rows*dy, (cols + 0.5)*dx)))
        signs.append(sgn)
    a_pts = np.concatenate(pts); sign = np.concatenate(signs)
    b_pts = boundary_points(candidate, spacing)
    if not len(a_pts) or not len(b_pts):
        return {"mean_signed": None, "mean_abs": None, "points": int(len(a_pts))}
    dist, _ = cKDTree(b_pts).query(a_pts)
    return {"mean_signed": float(np.mean(sign*dist)), "mean_abs": float(dist.mean()),
            "frac_plus_grew": float(np.mean(sign == 1)), "frac_minus_grew": float(np.mean(sign == -1)),
            "points": int(len(a_pts))}


def entropy_box_sweep(labels, box_sizes=(2, 4, 5, 8, 10, 20)):
    """Sb/Sbb for several box widths (only those dividing the grid)."""
    n = np.asarray(labels).shape[0]
    return {int(b): entropy(labels, b) for b in box_sizes if n % b == 0}


def uncertainty_exponent(classify_ics, ics, eps_values, seed=0):
    """Grebogi-McDonald-Ott-Yorke: f(eps) = fraction of ICs whose outcome
    changes under a random +/-eps perturbation; f ~ eps^alpha.
    classify_ics: callable (N,2) -> labels in {-1,0,1}. Expensive: each eps
    needs two classifications of all ICs (reuse the base labels)."""
    rng = np.random.default_rng(seed)
    base = classify_ics(ics)
    out = []
    for eps in eps_values:
        direction = rng.normal(size=ics.shape)
        direction /= np.linalg.norm(direction, axis=1, keepdims=True)
        pert = classify_ics(ics + eps*direction)
        ok = (base != 0) & (pert != 0)
        out.append(float(np.mean(base[ok] != pert[ok])) if ok.any() else np.nan)
    f = np.asarray(out)
    good = f > 0
    alpha = (float(np.polyfit(np.log(np.asarray(eps_values)[good]), np.log(f[good]), 1)[0])
             if good.sum() >= 3 else None)
    return {"eps": list(map(float, eps_values)), "f": out, "alpha": alpha}
