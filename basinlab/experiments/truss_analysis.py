"""Reanalyse stored truss labels on the same initial-condition set used for D."""
import json
from pathlib import Path

import numpy as np

from ..experiments.energy_loss_retest import interpolate_crossing
from ..metrics import mismatch
from ..truss import Truss


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "truss_basins"
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    configuration = summary["configuration"]
    model = Truss(half_span=configuration["half_span"],
                  rise=configuration["rise"],
                  axial_stiffness=configuration["axial_stiffness_EA_over_L0"],
                  damping=configuration["viscous_damping_per_coordinate"])
    cases = []
    for case in summary["cases"]:
        n, name = case["n"], case["model"]
        with np.load(out / f"{name}_n{n}.npz") as arrays:
            y, vy = arrays["y"], arrays["vy"]
            q = np.stack((np.full_like(y, case["horizontal_initial"]), y), axis=-1)
            v = np.stack((np.zeros_like(vy), vy), axis=-1)
            e0 = model.energy(q, v)
            selected = (e0 > 1.5*model.barrier) & (e0 < 8*model.barrier)
            reference = arrays["reference"][selected]
            for h in (.2, .4, .6):
                group = [r for r in case["runs"] if r["h"] == h]
                errors = []
                for row in group:
                    labels = arrays[f"h{h}_rho{row['rho']}"][selected]
                    errors.append(mismatch(reference, labels)["fraction"])
                rhos = np.asarray([row["rho"] for row in group])
                crossing = interpolate_crossing(
                    rhos, np.asarray([row["D_median"] for row in group]))
                best = int(np.argmin(errors))
                cases.append({"model": name, "n": n, "h": h,
                              "D_crossing_rho": crossing,
                              "D_time": group[0]["D_horizon"],
                              "same_subset_basin_best_rho": float(rhos[best]),
                              "same_subset_min_mismatch": float(errors[best]),
                              "same_subset_mismatches": errors,
                              "rho": rhos.tolist(),
                              "D_selection_count": int(selected.sum())})
                print(f"{name} n={n} h={h}: D*={crossing}, "
                      f"subset optimum={rhos[best]}, mismatch={errors[best]:.4f}")
    (out / "same_population_analysis.json").write_text(
        json.dumps({"cases": cases}, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
