"""Refine rho minima where matched-time D and basin optima separated."""
import csv
import json
import time
from pathlib import Path

import numpy as np

from ..integrators import parameters
from ..metrics import boundary_distance, entropy, mismatch
from ..model import classify, energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid, interpolate_crossing


CASES = (
    (.05, .5, np.round(np.arange(.50, .851, .025), 3)),
    (.05, .6, np.round(np.arange(.50, .851, .025), 3)),
    (.15, .4, np.round(np.arange(.25, .601, .025), 3)),
    (.15, .5, np.round(np.arange(.25, .601, .025), 3)),
    (.15, .6, np.round(np.arange(.25, .601, .025), 3)),
)


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "focused_rho_retest"
    out.mkdir(parents=True, exist_ok=True)
    n, energy_end, basin_end, reference_h = 100, 8.0, 60.0, .005
    q = -2 + (np.arange(n) + .5) * 4 / n
    X, V = np.meshgrid(q, q)
    E0 = energy(X, V)
    selected = (E0 > .6) & (E0 < 2.0)
    XS, VS, E0S = X[selected], V[selected], E0[selected]
    rows = []
    summary = {"configuration": {"grid_n": n, "damping_end": energy_end,
                                 "basin_end": basin_end, "reference_h": reference_h,
                                 "rho_refinement": "0.025 near observed minima",
                                 "mismatch_subset": ".6 < E0 < 2.0",
                                 "entropy_box_cells": 10},
               "cases": []}

    for damping, h, rhos in CASES:
        ref_e = rk4(XS, VS, energy_end, reference_h, damping)
        delta_ref = E0S - energy(*ref_e)
        valid_d = np.isfinite(delta_ref) & (delta_ref > 1e-12)
        ref_state = rk4(X, V, basin_end, reference_h, damping)
        ref_labels = classify(*ref_state)
        ref_subset = ref_labels[selected]
        ref_entropy = entropy(ref_labels, 10)
        case_rows = []
        for rho in rhos:
            p = parameters("generalized-alpha", float(rho))
            e_state = exact_step_grid(XS, VS, energy_end, h, damping, p)
            delta_num = E0S - energy(*e_state)
            ratio = delta_num[valid_d] / delta_ref[valid_d]
            state = exact_step_grid(X, V, basin_end, h, damping, p)
            labels = classify(*state)
            d_full, d_subset = mismatch(ref_labels, labels), mismatch(ref_subset, labels[selected])
            bd = boundary_distance(ref_labels, labels, (4/n, 4/n))
            ent = entropy(labels, 10)
            row = {"damping": damping, "h": h, "rho": float(rho),
                   "D_median": float(np.median(ratio)),
                   "D_q10": float(np.quantile(ratio, .1)),
                   "D_q90": float(np.quantile(ratio, .9)),
                   "basin_mismatch_full": d_full["fraction"],
                   "basin_mismatch_energy_subset": d_subset["fraction"],
                   "boundary_mean_symmetric": bd["mean_symmetric"],
                   "boundary_hausdorff": bd["hausdorff"],
                   "reference_Sb": ref_entropy["Sb"], "candidate_Sb": ent["Sb"],
                   "reference_Sbb": ref_entropy["Sbb"], "candidate_Sbb": ent["Sbb"]}
            rows.append(row)
            case_rows.append(row)

        rho_values = np.asarray([r["rho"] for r in case_rows])
        d_values = np.asarray([r["D_median"] for r in case_rows])
        mis_values = np.asarray([r["basin_mismatch_energy_subset"] for r in case_rows])
        d_star = interpolate_crossing(rho_values, d_values)
        best = int(np.argmin(mis_values))
        case = {"damping": damping, "h": h, "D_crossing_rho": d_star,
                "subset_basin_best_rho": float(rho_values[best]),
                "subset_basin_minimum": float(mis_values[best]),
                "rho_points": case_rows}
        summary["cases"].append(case)
        print(f"c={damping:.2f} h={h:.1f}: D*={d_star}, "
              f"same-subset basin best={rho_values[best]:.3f}, "
              f"mismatch={mis_values[best]:.4f}", flush=True)

    with (out / "rho_refinement.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary["runtime_seconds"] = round(time.perf_counter() - started, 3)
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} points; runtime {summary['runtime_seconds']} s", flush=True)


if __name__ == "__main__":
    main()
