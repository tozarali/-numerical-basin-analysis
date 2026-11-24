"""Run the small, reproducible Phase 1 convergence experiment."""
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from .integrators import integrate, parameters
from .metrics import boundary_distance, entropy, mismatch
from .model import classify
from .reference import adaptive, rk4


def grid(n):
    q = -2 + (np.arange(n) + 0.5) * 4 / n
    return np.meshgrid(q, q)


def compare_states(first, second):
    return float(np.max(np.abs(np.asarray(first) - np.asarray(second))))


def main():
    start = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    out = root / "data"
    out.mkdir(exist_ok=True)
    damping, end = 0.15, 80.0
    report = {
        "configuration": {"domain": [-2, 2, -2, 2], "grid": "cell centers",
                          "damping": damping, "end": end, "box_width": 1.0,
                          "candidate": "generalized-alpha, rho=0.4, h=0.4",
                          "python": platform.python_version(), "numpy": np.__version__,
                          "scipy": scipy.__version__},
        "independent_reference": {}, "grids": []}

    x, v = grid(16)
    low = adaptive(x, v, end, damping, rtol=1e-10, atol=1e-12)
    high = adaptive(x, v, end, damping, rtol=1e-12, atol=1e-14)
    other = adaptive(x, v, end, damping, method="RK45", rtol=1e-11, atol=1e-13)
    fixed = rk4(x, v, end, .005, damping)
    labels = classify(*high)
    report["independent_reference"] = {
        "grid_n": 16,
        "DOP853_1e-10_vs_1e-12": mismatch(labels, classify(*low)),
        "DOP853_vs_RK45": mismatch(labels, classify(*other)),
        "DOP853_vs_RK4_h_0.005": mismatch(labels, classify(*fixed)),
        "max_state_difference_DOP853_RK45": compare_states(high, other),
        "max_state_difference_DOP853_RK4": compare_states(high, fixed)}

    for n in (24, 48, 96):
        x, v = grid(n)
        ref_coarse = rk4(x, v, end, .01, damping)
        ref = rk4(x, v, end, .005, damping)
        ref_long = rk4(x, v, 120.0, .01, damping)
        candidate = integrate(x, v, end, .4, damping,
                              parameters("generalized-alpha", .4))
        candidate_long = integrate(x, v, 120.0, .4, damping,
                                   parameters("generalized-alpha", .4))
        r, rc, rl = classify(*ref), classify(*ref_coarse), classify(*ref_long)
        z, zl = classify(*candidate), classify(*candidate_long)
        np.savez_compressed(out / f"phase1_n{n}.npz",
                            x=x, v=v, reference=r, candidate=z,
                            reference_coarse=rc, reference_long=rl, candidate_long=zl)
        spacing = (4/n, 4/n)
        row = {"n": n, "points": n*n,
               "reference_h_convergence": mismatch(r, rc),
               "reference_time_convergence": mismatch(r, rl),
               "candidate_time_convergence": mismatch(z, zl),
               "candidate_mismatch": mismatch(r, z),
               "reference_entropy": entropy(r, n//4),
               "candidate_entropy": entropy(z, n//4)}
        if np.all(r != 0) and np.all(z != 0):
            row["boundary_distance"] = boundary_distance(r, z, spacing)
        else:
            row["boundary_distance"] = None
        report["grids"].append(row)
        print(f"n={n}: mismatch={row['candidate_mismatch']}, "
              f"Sb={row['reference_entropy']['Sb']}, "
              f"boundary={row['boundary_distance']}", flush=True)

    report["runtime_seconds"] = round(time.perf_counter() - start, 3)
    (out / "phase1_results.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Results written to {out}; runtime {report['runtime_seconds']} s", flush=True)


if __name__ == "__main__":
    main()
