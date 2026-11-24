"""Re-test the D=1 prediction with exact common end times and basin labels."""
import csv
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from ..integrators import parameters, step
from ..metrics import mismatch
from ..model import classify, energy, force
from ..reference import rk4


def exact_step_grid(x, v, end, h, damping, p):
    """Use h-sized steps and one shorter last step so all schemes end at end."""
    if h <= 0 or end <= 0:
        raise ValueError("h and end must be positive")
    x, v = np.asarray(x, dtype=float).copy(), np.asarray(v, dtype=float).copy()
    a = -damping * v - force(x)
    full_steps = int(np.floor(end / h + 1e-12))
    for _ in range(full_steps):
        x, v, a = step(x, v, a, h, damping, p)
    remainder = end - full_steps * h
    if remainder > 1e-10:
        x, v, a = step(x, v, a, remainder, damping, p)
    return x, v


def interpolate_crossing(rhos, values, target=1.0):
    crossings = np.flatnonzero((values[:-1] - target) * (values[1:] - target) <= 0)
    if not len(crossings):
        return None
    i = int(crossings[0])
    if values[i + 1] == values[i]:
        return float((rhos[i] + rhos[i + 1]) / 2)
    return float(rhos[i] + (target - values[i]) *
                 (rhos[i + 1] - rhos[i]) / (values[i + 1] - values[i]))


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "energy_loss_retest"
    out.mkdir(parents=True, exist_ok=True)

    damping_values = (.05, .15)
    h_values = (.2, .3, .4, .5, .6)
    rhos = np.round(np.arange(0, 1.0001, .05), 2)
    energy_end, basin_end, reference_h, n = 8.0, 60.0, .005, 100
    q = -2 + (np.arange(n) + .5) * 4 / n
    X, V = np.meshgrid(q, q)
    E0 = energy(X, V)
    selected = (E0 > .6) & (E0 < 2.0)
    XS, VS, E0S = X[selected], V[selected], E0[selected]

    rows = []
    summary = {"configuration": {
        "model": "x'' + c*x' - x + x^3 = 0",
        "domain": [-2, 2, -2, 2], "grid": "100x100 cell centers",
        "energy_initial_subset": ".6 < E0 < 2.0", "D_end_time": energy_end,
        "basin_end_time": basin_end, "reference": "vectorized RK4",
        "reference_h": reference_h, "rho_step": .05,
        "last_step": "shortened to hit the exact requested end time",
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__},
        "cases": []}

    for damping in damping_values:
        ref_energy = rk4(XS, VS, energy_end, reference_h, damping)
        delta_ref = E0S - energy(*ref_energy)
        valid_d = np.isfinite(delta_ref) & (delta_ref > 1e-12)
        ref_basin = rk4(X, V, basin_end, reference_h, damping)
        labels_ref = classify(*ref_basin)
        labels_long = classify(*rk4(X, V, 2 * basin_end, reference_h, damping))
        ref_time_check = mismatch(labels_ref, labels_long)

        for h in h_values:
            case = {"damping": damping, "h": h, "reference_time_check": ref_time_check,
                    "rho": []}
            d_medians = []
            mismatch_values = []
            for rho in rhos:
                p = parameters("generalized-alpha", float(rho))
                energy_state = exact_step_grid(XS, VS, energy_end, h, damping, p)
                delta_num = E0S - energy(*energy_state)
                ratio = delta_num[valid_d] / delta_ref[valid_d]
                basin_state = exact_step_grid(X, V, basin_end, h, damping, p)
                labels = classify(*basin_state)
                mis = mismatch(labels_ref, labels)
                row = {
                    "damping": damping, "h": h, "rho": float(rho),
                    "D_median": float(np.median(ratio)),
                    "D_q10": float(np.quantile(ratio, .1)),
                    "D_q90": float(np.quantile(ratio, .9)),
                    "D_valid_ICs": int(valid_d.sum()),
                    "basin_mismatch": mis["fraction"],
                    "basin_valid": mis["valid"],
                    "basin_excluded": mis["excluded"]}
                rows.append(row)
                case["rho"].append(row)
                d_medians.append(row["D_median"])
                mismatch_values.append(row["basin_mismatch"])

            d_medians = np.asarray(d_medians)
            mismatch_values = np.asarray(mismatch_values)
            d_star = interpolate_crossing(rhos, d_medians)
            j_d = int(np.argmin(np.abs(d_medians - 1)))
            j_b = int(np.argmin(mismatch_values))
            case["D_crossing_rho"] = d_star
            case["D_closest_rho"] = float(rhos[j_d])
            case["basin_best_rho"] = float(rhos[j_b])
            case["rho_difference_crossing_minus_basin_best"] = (
                None if d_star is None else float(d_star - rhos[j_b]))
            summary["cases"].append(case)
            print(f"c={damping:.2f} h={h:.1f}: D*= {d_star}; "
                  f"nearest={rhos[j_d]:.2f}, basin best={rhos[j_b]:.2f}, "
                  f"mismatch={mismatch_values[j_b]:.4f}", flush=True)

    with (out / "rho_sweep.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary["runtime_seconds"] = round(time.perf_counter() - started, 3)
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} rows; runtime {summary['runtime_seconds']} s", flush=True)


if __name__ == "__main__":
    main()
