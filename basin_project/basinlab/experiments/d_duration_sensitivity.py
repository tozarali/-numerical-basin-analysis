"""Check whether the empirical D=1 crossing depends on its chosen horizon."""
import csv
import json
import time
from pathlib import Path

import numpy as np

from ..integrators import parameters
from ..model import energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid, interpolate_crossing


CASES = ((.05, .4), (.05, .5), (.05, .6),
         (.15, .4), (.15, .5), (.15, .6))
HORIZONS = (2.0, 4.0, 8.0, 12.0, 16.0)


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "d_duration_sensitivity"
    out.mkdir(parents=True, exist_ok=True)
    n, reference_h = 100, .005
    q = -2 + (np.arange(n) + .5) * 4 / n
    X, V = np.meshgrid(q, q)
    E0 = energy(X, V)
    selected = (E0 > .6) & (E0 < 2)
    X, V, E0 = X[selected], V[selected], E0[selected]
    rhos = np.round(np.arange(0, 1.0001, .05), 2)
    rows, summaries = [], []

    for damping, h in CASES:
        for horizon in HORIZONS:
            ref_state = rk4(X, V, horizon, reference_h, damping)
            delta_ref = E0 - energy(*ref_state)
            valid = np.isfinite(delta_ref) & (delta_ref > 1e-12)
            D = []
            for rho in rhos:
                state = exact_step_grid(X, V, horizon, h, damping,
                                        parameters("generalized-alpha", float(rho)))
                ratio = (E0 - energy(*state))[valid] / delta_ref[valid]
                median = float(np.median(ratio))
                D.append(median)
                rows.append({"damping": damping, "h": h, "horizon": horizon,
                             "rho": float(rho), "D_median": median,
                             "D_q10": float(np.quantile(ratio, .1)),
                             "D_q90": float(np.quantile(ratio, .9)),
                             "valid_ICs": int(valid.sum())})
            crossing = interpolate_crossing(rhos, np.asarray(D))
            summaries.append({"damping": damping, "h": h, "horizon": horizon,
                              "D_crossing_rho": crossing,
                              "valid_ICs": int(valid.sum())})
            print(f"c={damping:.2f} h={h:.1f} T={horizon:g}: D*={crossing}",
                  flush=True)

    with (out / "D_by_horizon.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {"initial_condition_region": ".6 < E0 < 2", "grid_n": n,
              "reference": "RK4, h=.005", "rho_step": .05,
              "results": summaries, "runtime_seconds": round(time.perf_counter() - started, 3)}
    (out / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} points; runtime {result['runtime_seconds']} s",
          flush=True)


if __name__ == "__main__":
    main()
