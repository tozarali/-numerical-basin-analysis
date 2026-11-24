"""Grid-based basin metrics with unknown outcomes excluded explicitly."""
import numpy as np
from scipy.spatial import cKDTree


def _validate(labels):
    labels = np.asarray(labels)
    if labels.ndim != 2 or not np.all(np.isin(labels, [-1, 0, 1])):
        raise ValueError("labels must be a 2D array containing only -1, 0, 1")
    return labels


def mismatch(reference, candidate):
    reference, candidate = np.asarray(reference), np.asarray(candidate)
    if reference.shape != candidate.shape:
        raise ValueError("label array shapes must match")
    if (reference.ndim == 0 or not np.all(np.isin(reference, [-1, 0, 1]))
            or not np.all(np.isin(candidate, [-1, 0, 1]))):
        raise ValueError("labels must contain only -1, 0, 1")
    valid = (reference != 0) & (candidate != 0)
    count = int(valid.sum())
    return {"fraction": float(np.mean(reference[valid] != candidate[valid])) if count else None,
            "valid": count, "excluded": int(valid.size - count)}


def entropy(labels, box_cells):
    labels = _validate(labels)
    if box_cells < 2 or any(n % box_cells for n in labels.shape):
        raise ValueError("box_cells must be >= 2 and divide both dimensions")
    entropies = []
    excluded = 0
    for row in range(0, labels.shape[0], box_cells):
        for col in range(0, labels.shape[1], box_cells):
            box = labels[row:row+box_cells, col:col+box_cells]
            if np.any(box == 0):
                excluded += 1
                continue
            p = np.count_nonzero(box == 1) / box.size
            entropy_value = 0.0
            if 0 < p < 1:
                entropy_value = -p*np.log(p) - (1-p)*np.log(1-p)
            entropies.append(entropy_value)
    values = np.asarray(entropies)
    boundary = values[values > 0]
    return {"Sb": float(values.mean()) if values.size else None,
            "Sbb": float(boundary.mean()) if boundary.size else None,
            "boxes": int(values.size), "boundary_boxes": int(boundary.size),
            "excluded_boxes": excluded}


def boundary_points(labels, spacing):
    labels = _validate(labels)
    dy, dx = spacing
    if dx <= 0 or dy <= 0:
        raise ValueError("grid spacing must be positive")
    # Midpoints of adjacent, resolved, different labels; coordinates are relative.
    h = (labels[:, :-1] * labels[:, 1:] == -1)
    v = (labels[:-1, :] * labels[1:, :] == -1)
    rh, ch = np.nonzero(h)
    rv, cv = np.nonzero(v)
    return np.concatenate((np.column_stack((rh*dy, (ch+0.5)*dx)),
                           np.column_stack(((rv+0.5)*dy, cv*dx))))


def boundary_distance(reference, candidate, spacing):
    reference, candidate = _validate(reference), _validate(candidate)
    if reference.shape != candidate.shape:
        raise ValueError("grid shapes must match")
    if np.any(reference == 0) or np.any(candidate == 0):
        raise ValueError("boundary distances require fully resolved grids")
    a, b = boundary_points(reference, spacing), boundary_points(candidate, spacing)
    if not len(a) or not len(b):
        return {"mean_symmetric": None, "hausdorff": None,
                "reference_points": len(a), "candidate_points": len(b)}
    ab = cKDTree(b).query(a)[0]
    ba = cKDTree(a).query(b)[0]
    return {"mean_symmetric": float((ab.mean() + ba.mean()) / 2),
            "hausdorff": float(max(ab.max(), ba.max())),
            "reference_points": len(a), "candidate_points": len(b)}
