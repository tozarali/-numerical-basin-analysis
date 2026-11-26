"""Re-test family collapse with matched end times and common IC populations."""
import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from ..integrators import parameters_for_spectral_radius
from ..metrics import boundary_distance, mismatch
from ..model import classify, energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid, interpolate_crossing


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "family_collapse_retest"
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
    targets_by_family = {
        "generalized-alpha": np.round(np.arange(0, 1.0001, .05), 2),
        "HHT": np.round(np.arange(.55, 1.0001, .05), 2),
        "WBZ": np.round(np.arange(.05, 1.0001, .05), 2),
        "Newmark": np.round(np.arange(.6, 1.0001, .05), 2)}
    rows, group_stats = [], []

    for family, targets in targets_by_family.items():
        family_rows = []
        for target in targets:
            p, actual_rho = parameters_for_spectral_radius(family, float(target))
            e_state = exact_step_grid(XS, VS, energy_end, h, damping, p)
            ratio = ((E0S - energy(*e_state))[valid_d] / delta_ref[valid_d])
            basin_state = exact_step_grid(X, V, basin_end, h, damping, p)
            labels = classify(*basin_state)
            full = mismatch(ref_labels, labels)
            subset = mismatch(ref_labels[selected], labels[selected])
            bdist = boundary_distance(ref_labels, labels, (4/n, 4/n))
            row = {"family": family, "rho_target": float(target),
                   "rho_measured": actual_rho, "D_median": float(np.median(ratio)),
                   "D_q10": float(np.quantile(ratio, .1)),
                   "D_q90": float(np.quantile(ratio, .9)),
                   "mismatch_full": full["fraction"],
                   "mismatch_energy_subset": subset["fraction"],
                   "boundary_mean_symmetric": bdist["mean_symmetric"],
                   "boundary_hausdorff": bdist["hausdorff"]}
            rows.append(row)
            family_rows.append(row)
        ds = np.asarray([r["D_median"] for r in family_rows])
        mismatches = np.asarray([r["mismatch_energy_subset"] for r in family_rows])
        group_stats.append({
            "family": family, "D_crossing_rho": interpolate_crossing(targets, ds),
            "basin_best_rho": float(targets[int(np.argmin(mismatches))]),
            "basin_minimum": float(np.min(mismatches)),
            "spearman_mismatch_abs_D_error": float(spearmanr(np.abs(ds - 1), mismatches).statistic),
            "spearman_mismatch_one_minus_rho": float(spearmanr(1-targets, mismatches).statistic)})
        print(f"{family}: D*={group_stats[-1]['D_crossing_rho']}, "
              f"basin best={group_stats[-1]['basin_best_rho']}, "
              f"min subset mismatch={group_stats[-1]['basin_minimum']:.4f}",
              flush=True)

    arr_d = np.asarray([abs(r["D_median"] - 1) for r in rows])
    arr_e = np.asarray([r["mismatch_energy_subset"] for r in rows])
    arr_r = np.asarray([1-r["rho_target"] for r in rows])
    result = {"configuration": {"model": "x'' + c*x' - x + x^3 = 0",
                                "damping": damping, "h": h,
                                "energy_end": energy_end, "basin_end": basin_end,
                                "grid_n": n, "reference_h": reference_h,
                                "rho_ranges": {name: [float(values[0]), float(values[-1])]
                                               for name, values in targets_by_family.items()},
                                "same_energy_subset_for_D_and_mismatch": True},
              "families": group_stats,
              "pooled": {
                  "spearman_mismatch_abs_D_error": float(spearmanr(arr_d, arr_e).statistic),
                  "spearman_mismatch_one_minus_rho": float(spearmanr(arr_r, arr_e).statistic)},
              "runtime_seconds": round(time.perf_counter() - started, 3)}
    with (out / "rho_sweep.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (out / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} matched-family points; runtime {result['runtime_seconds']} s",
          flush=True)


if __name__ == "__main__":
    main()
