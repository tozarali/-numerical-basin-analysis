"""WP2 (exploratory): basin mismatch over (c, h, rho, family) on a fixed
100x100 grid, plus BDF2 and discrete-gradient controls. Resumable: one part
file per (c, h). Labels are kept (int8) for diagnostic_eval.
Not preregistered; asymmetric well (epsilon) is NOT included because the
asymmetric model/barrier is not yet implemented and validated."""
import numpy as np

from ..metrics import entropy, mismatch
from ..metrics_ext import net_basin_transfer
from ..model import classify
from ..reference import rk4
from ..runmeta import experiment
from ._common import (FAMILY_GRIDS, Parts, alpha_state, control_labels, family_parameters,
                      grid, out_name, quick, rho_grid, saturated)

T, T_CHECK, REF_H = 60.0, 120.0, 0.005
CONTROLS = ["BDF2", "discrete-gradient"]


def config():
    if quick():
        return {"n": 24, "c": [0.15], "h": [0.2, 0.4], "rho_step": 0.25,
                "families": list(FAMILY_GRIDS)}
    return {"n": 100, "c": [0.05, 0.1, 0.15, 0.25, 0.4],
            "h": [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8], "rho_step": 0.05,
            "families": list(FAMILY_GRIDS)}


def cases(cfg):
    out = []
    for fam in cfg["families"]:
        out += [(fam, r) for r in rho_grid(fam, cfg["rho_step"])]
    return out + [(c, None) for c in CONTROLS]


def main():
    cfg = dict(config(), T=T, T_check=T_CHECK, reference_h=REF_H, box=[-2, 2],
               saturation="mismatch>=0.30 or candidate Sb>=0.9 ln2 (box 0.4)",
               candidate_time_check_rhos=[0.0, 0.5, 1.0])
    case_list = cases(cfg)
    with experiment(out_name("regime_map"), cfg) as run:
        parts = Parts(run["out_dir"])
        n = cfg["n"]
        X, V = grid(n)
        box = max(2, int(round(0.4*n/4)))
        while n % box:
            box += 1
        cells = []
        for c in cfg["c"]:
            rkey = f"ref_c{c}"
            if parts.has(rkey):
                rmeta, arr = parts.load(rkey); ref = arr["labels"]
            else:
                ref = classify(*rk4(X, V, T, REF_H, c))
                ref2 = classify(*rk4(X, V, T_CHECK, REF_H, c))
                rmeta = {"c": c, "time_check_changes": int((ref != ref2).sum()),
                         "unresolved": int((ref == 0).sum())}
                parts.save(rkey, rmeta, labels=ref)
            for h in cfg["h"]:
                key = f"c{c}_h{h}"
                if parts.has(key):
                    cells += parts.load(key)[0]["cells"]; continue
                labels = np.zeros((len(case_list), n, n), np.int8)
                rows = []
                for i, (fam, rho) in enumerate(case_list):
                    if rho is None:
                        lab = control_labels(fam, X, V, T, h, c)
                        measured = None
                    else:
                        p, measured = family_parameters(fam, rho)
                        lab = classify(*alpha_state(X, V, T, h, c, p))
                    labels[i] = lab
                    mm = mismatch(ref, lab)
                    ent = entropy(lab, box)
                    row = {"c": c, "h": h, "family": fam, "rho": rho, "rho_measured": measured,
                           "mismatch": mm["fraction"], "valid": mm["valid"],
                           "excluded": mm["excluded"], "Sb": ent["Sb"],
                           "net_transfer": net_basin_transfer(ref, lab)["net"],
                           "reference_time_check_changes": rmeta["time_check_changes"]}
                    if rho is not None and any(abs(rho - r) < 1e-9 for r in cfg["candidate_time_check_rhos"]):
                        lab2 = classify(*alpha_state(X, V, T_CHECK, h, c, p))
                        row["candidate_time_check_changes"] = int((lab != lab2).sum())
                    row["saturated"] = saturated(row["mismatch"], row["Sb"])
                    rows.append(row)
                parts.save(key, {"cells": rows, "cases": [[f, r] for f, r in case_list]},
                           labels=labels)
                cells += rows
                best = min((r for r in rows if r["rho"] is not None and r["mismatch"] is not None),
                           key=lambda r: r["mismatch"])
                print(f"c={c} h={h}: best {best['family']} rho={best['rho']} "
                      f"mismatch={best['mismatch']:.3f}", flush=True)
        # per (c,h,family) discrete optimum among non-saturated cells
        optima = []
        for c in cfg["c"]:
            for h in cfg["h"]:
                for fam in cfg["families"]:
                    fr = [r for r in cells if r["c"] == c and r["h"] == h and r["family"] == fam
                          and r["mismatch"] is not None]
                    if not fr:
                        continue
                    b = min(fr, key=lambda r: r["mismatch"])
                    optima.append({"c": c, "h": h, "family": fam, "rho_best": b["rho"],
                                   "mismatch_best": b["mismatch"],
                                   "mismatch_rho1": next((r["mismatch"] for r in fr if r["rho"] == 1.0), None),
                                   "best_cell_saturated": b["saturated"],
                                   "all_saturated": all(r["saturated"] for r in fr)})
        run["results"] = {"n_cells": len(cells), "cells": cells, "family_optima": optima,
                          "note": "discrete grid optima, no uncertainty; exploratory only"}
    return run


if __name__ == "__main__":
    main()
