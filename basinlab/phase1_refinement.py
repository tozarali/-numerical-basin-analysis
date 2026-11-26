"""Time-step refinement using the previously saved Phase 1 reference grid."""
import json
import time
from pathlib import Path

import numpy as np

from .integrators import integrate, parameters
from .metrics import mismatch
from .model import classify


def main():
    start = time.perf_counter()
    out = Path(__file__).resolve().parents[1] / "data"
    configuration = json.loads((out / "phase1_results.json").read_text(
        encoding="utf-8"))["configuration"]
    if configuration["damping"] != .15 or configuration["end"] != 80.0:
        raise ValueError("Phase 1 reference configuration does not match refinement")
    with np.load(out / "phase1_n48.npz") as grid:
        x, v, reference = grid["x"], grid["v"], grid["reference"]
    results = []
    for h in (.4, .2, .1, .05):
        state = integrate(x, v, 80, h, .15, parameters("generalized-alpha", .4))
        results.append({"h": h, "mismatch": mismatch(reference, classify(*state))})
    result = {"grid_n": 48, "damping": .15, "end": 80, "rho": .4,
              "steps": results, "runtime_seconds": round(time.perf_counter() - start, 3)}
    (out / "phase1_time_refinement.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
