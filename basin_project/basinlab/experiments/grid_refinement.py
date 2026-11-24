"""Recheck selected D and basin optima on a denser initial-condition grid."""
import json
import time
from pathlib import Path

import numpy as np

from ..integrators import parameters
from ..metrics import boundary_distance, entropy, mismatch
from ..model import classify, energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid, interpolate_crossing


CASES = ((.15, .6, (.375, .40, .425, .45, .475, .50, .525)),
         (.05, .5, (.60, .65, .70, .75, .80)))


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "grid_refinement"
    out.mkdir(parents=True, exist_ok=True)
    n, reference_h, basin_end, energy_end = 200, .005, 60.0, 8.0
    q = -2 + (np.arange(n) + .5) * 4 / n
    X, V = np.meshgrid(q, q)
    E0 = energy(X, V)
    selected = (E0 > .6) & (E0 < 2)
    XS, VS, E0S = X[selected], V[selected], E0[selected]
    summaries = []

    for damping, h, rhos in CASES:
        ref_e = rk4(XS, VS, energy_end, reference_h, damping)
        delta_ref = E0S - energy(*ref_e)
        valid_d = np.isfinite(delta_ref) & (delta_ref > 1e-12)
        ref_labels = classify(*rk4(X, V, basin_end, reference_h, damping))
        ref_long = classify(*rk4(X, V, 2*basin_end, reference_h, damping))
        time_check = mismatch(ref_labels, ref_long)
        rows = []
        for rho in rhos:
            p = parameters("generalized-alpha", rho)
            e_state = exact_step_grid(XS, VS, energy_end, h, damping, p)
            ratio = ((E0S - energy(*e_state))[valid_d] / delta_ref[valid_d])
            state = exact_step_grid(X, V, basin_end, h, damping, p)
            labels = classify(*state)
            full = mismatch(ref_labels, labels)
            subset = mismatch(ref_labels[selected], labels[selected])
            bdist = boundary_distance(ref_labels, labels, (4/n, 4/n))
            row = {"rho": rho, "D_median": float(np.median(ratio)),
                   "D_q10": float(np.quantile(ratio, .1)),
                   "D_q90": float(np.quantile(ratio, .9)),
                   "full_mismatch": full["fraction"],
                   "subset_mismatch": subset["fraction"],
                   "boundary_mean_symmetric": bdist["mean_symmetric"],
                   "boundary_hausdorff": bdist["hausdorff"]}
            rows.append(row)
            print(f"n={n} c={damping:.2f} h={h:.1f} rho={rho:.3f}: "
                  f"D={row['D_median']:.3f}, subset mismatch={row['subset_mismatch']:.4f}",
                  flush=True)
        d_star = interpolate_crossing(np.asarray(rhos),
                                      np.asarray([r["D_median"] for r in rows]))
        best = min(rows, key=lambda row: row["subset_mismatch"])
        summaries.append({"damping": damping, "h": h, "grid_n": n,
                          "reference_T60_T120_mismatch": time_check,
                          "D_crossing_rho": d_star,
                          "subset_basin_best_rho": best["rho"],
                          "points": rows})

    result = {"grid": "200x200 cell centers", "reference_h": reference_h,
              "cases": summaries, "runtime_seconds": round(time.perf_counter() - started, 3)}
    (out / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved dense-grid comparison; runtime {result['runtime_seconds']} s",
          flush=True)


if __name__ == "__main__":
    main()
