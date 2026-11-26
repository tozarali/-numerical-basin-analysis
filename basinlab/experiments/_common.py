"""Shared helpers for the new experiments (quick mode, cached family maps,
fixed-end-time SDOF label runners, part-file checkpoints)."""
import functools
import json
import os
from pathlib import Path

import numpy as np

from ..integrators import parameters_for_spectral_radius, step
from ..integrators_extra import bdf2_integrate, discrete_gradient_step
from ..model import classify, force

ROOT = Path(__file__).resolve().parents[2]


def quick():
    """BASIN_QUICK=1 -> smoke sizes; outputs go to data/<name>_quick."""
    return os.environ.get("BASIN_QUICK", "0") == "1"


def out_name(name):
    return name + "_quick" if quick() else name


@functools.lru_cache(maxsize=None)
def family_parameters(family, rho):
    """Cached (p, measured_rho) or None if the family cannot attain rho."""
    try:
        return parameters_for_spectral_radius(family, float(rho))
    except ValueError:
        return None


FAMILY_GRIDS = {  # admissible targets are filtered by family_parameters
    "generalized-alpha": (0.0, 1.0),
    "HHT": (0.5, 1.0),
    "WBZ": (0.0, 1.0),
    "Newmark": (0.6, 1.0),
}


def rho_grid(family, step_size):
    lo, hi = FAMILY_GRIDS[family]
    grid = np.round(np.arange(lo, hi + 1e-9, step_size), 4)
    return [float(r) for r in grid if family_parameters(family, float(r)) is not None]


def alpha_state(x, v, end, h, damping, p):
    """Alpha-family state at exactly `end` (shortened last step)."""
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    a = -damping * v - force(x)
    n = int(np.floor(end / h + 1e-12))
    for _ in range(n):
        x, v, a = step(x, v, a, h, damping, p)
    rem = end - n * h
    if rem > 1e-10:
        x, v, a = step(x, v, a, rem, damping, p)
    return x, v


def discrete_gradient_state(x, v, end, h, damping):
    x, v = np.asarray(x, float).copy(), np.asarray(v, float).copy()
    n = int(np.floor(end / h + 1e-12))
    for _ in range(n):
        x, v = discrete_gradient_step(x, v, h, damping)
    rem = end - n * h
    if rem > 1e-10:
        x, v = discrete_gradient_step(x, v, rem, damping)
    return x, v


def control_labels(name, x, v, end, h, damping):
    if name == "BDF2":
        return classify(*bdf2_integrate(x, v, end, h, damping))
    if name == "discrete-gradient":
        return classify(*discrete_gradient_state(x, v, end, h, damping))
    raise ValueError(name)


def grid(n, lo=-2.0, hi=2.0):
    q = lo + (np.arange(n) + 0.5) * (hi - lo) / n
    return np.meshgrid(q, q)


class Parts:
    """Resumable part files: data/<name>/parts/<key>.npz (+ .json)."""

    def __init__(self, out_dir):
        self.dir = Path(out_dir) / "parts"
        self.dir.mkdir(parents=True, exist_ok=True)

    def has(self, key):
        return (self.dir / f"{key}.json").exists()

    def save(self, key, meta, **arrays):
        if arrays:
            np.savez_compressed(self.dir / f"{key}.npz", **arrays)
        tmp = self.dir / f"{key}.json.tmp"
        tmp.write_text(json.dumps(meta, indent=1, default=float), encoding="utf-8")
        tmp.replace(self.dir / f"{key}.json")      # json last = part complete

    def load(self, key):
        meta = json.loads((self.dir / f"{key}.json").read_text(encoding="utf-8"))
        npz = self.dir / f"{key}.npz"
        arrays = dict(np.load(npz)) if npz.exists() else {}
        return meta, arrays


def saturated(mismatch_fraction, sb=None, threshold=0.30):
    """Saturation flag used throughout (mismatch >= 30 % or high basin entropy)."""
    flag = mismatch_fraction is not None and mismatch_fraction >= threshold
    if sb is not None:
        flag = flag or sb >= 0.9 * np.log(2)
    return bool(flag)
