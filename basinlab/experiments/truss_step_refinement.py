"""Time-step refinement against the saved structural reference labels."""
import json
import time
from pathlib import Path

import numpy as np

from ..metrics import mismatch
from ..truss import Truss, generalized_alpha
from .truss_basins import grid


def main():
    started = time.perf_counter()
    out = Path(__file__).resolve().parents[2] / "data" / "truss_basins"
    source = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    cfg = source["configuration"]
    model = Truss(cfg["half_span"], cfg["rise"], cfg["axial_stiffness_EA_over_L0"],
                  cfg["viscous_damping_per_coordinate"])
    results = []
    for name, horizontal, constrained in (("1d_constrained", 0.0, True),
                                           ("2d_perturbed", .6, False)):
        q, v = grid(96, horizontal)
        with np.load(out / f"{name}_n96.npz") as data:
            reference = data["reference"]
            if not np.array_equal(data["y"], q[..., 1]) or not np.array_equal(data["vy"], v[..., 1]):
                raise ValueError(f"Stored {name} grid differs from expected grid")
        for h in (.2, .1, .05):
            candidate = generalized_alpha(model, q, v, 60, h, .4, constrained)
            labels = model.classify(*candidate)
            result = {"model": name, "n": 96, "h": h, "rho": .4,
                      "mismatch": mismatch(reference, labels)}
            results.append(result)
            print(f"{name} h={h}: {result['mismatch']}", flush=True)
    (out / "step_refinement.json").write_text(
        json.dumps({"results": results, "runtime_seconds": round(time.perf_counter()-started, 3)},
                   indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
