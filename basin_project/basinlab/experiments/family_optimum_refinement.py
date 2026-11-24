"""Fine rho refinement around each integrator family's D and basin optima."""
import csv
import json
import time
from pathlib import Path

import numpy as np

from ..integrators import parameters_for_spectral_radius
from ..metrics import boundary_distance, mismatch
from ..model import classify, energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid, interpolate_crossing


TARGETS = {
    "generalized-alpha": np.round(np.arange(.35, .6001, .025), 3),
    "HHT": np.round(np.arange(.51, .7001, .025), 3),
    "WBZ": np.round(np.arange(.50, .9001, .025), 3),
    "Newmark": np.round(np.arange(.80, 1.0001, .025), 3)}


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "family_optimum_refinement"
    out.mkdir(parents=True, exist_ok=True)
    damping, h, energy_end, basin_end = .15, .4, 4.0, 60.0
    n, reference_h = 100, .005
    q = -2 + (np.arange(n) + .5) * 4 / n
    X, V = np.meshgrid(q, q)
    E0 = energy(X, V)
    selected = (E0 > .6) & (E0 < 2)
    XS, VS, E0S = X[selected], V[selected], E0[selected]
    ref_e = rk4(XS, VS, energy_end, reference_h, damping)
    delta_ref = E0S - energy(*ref_e)
    valid_d = np.isfinite(delta_ref) & (delta_ref > 1e-12)
    ref_labels = classify(*rk4(X, V, basin_end, reference_h, damping))
    rows, group_stats = [], []

    for family, targets in TARGETS.items():
        family_rows = []
        for target in targets:
            p, actual_rho = parameters_for_spectral_radius(family, float(target))
            e_state = exact_step_grid(XS, VS, energy_end, h, damping, p)
            ratio = ((E0S - energy(*e_state))[valid_d] / delta_ref[valid_d])
            state = exact_step_grid(X, V, basin_end, h, damping, p)
            labels = classify(*state)
            full = mismatch(ref_labels, labels)
            subset = mismatch(ref_labels[selected], labels[selected])
            bd = boundary_distance(ref_labels, labels, (4/n, 4/n))
            row = {"family": family, "rho": float(target), "rho_measured": actual_rho,
                   "D_median": float(np.median(ratio)),
                   "mismatch_subset": subset["fraction"],
                   "mismatch_full": full["fraction"],
                   "boundary_mean_symmetric": bd["mean_symmetric"],
                   "boundary_hausdorff": bd["hausdorff"]}
            rows.append(row)
            family_rows.append(row)
        ds = np.asarray([r["D_median"] for r in family_rows])
        errors = np.asarray([r["mismatch_subset"] for r in family_rows])
        star = interpolate_crossing(targets, ds)
        best = int(np.argmin(errors))
        group_stats.append({"family": family, "D_crossing_rho": star,
                            "basin_best_rho": float(targets[best]),
                            "minimum_mismatch": float(errors[best]),
                            "crossing_minus_basin_best": (
                                None if star is None else float(star - targets[best]))})
        print(f"{family}: D*={star}, basin best={targets[best]:.3f}, "
              f"mismatch={errors[best]:.4f}", flush=True)

    with (out / "rho_refinement.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {"configuration": {"damping": damping, "h": h, "energy_end": energy_end,
                                "basin_end": basin_end, "grid_n": n,
                                "reference_h": reference_h,
                                "same_energy_subset_for_D_and_mismatch": True},
              "families": group_stats,
              "runtime_seconds": round(time.perf_counter() - started, 3)}
    (out / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} points; runtime {result['runtime_seconds']} s",
          flush=True)


if __name__ == "__main__":
    main()
