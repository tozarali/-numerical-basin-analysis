"""WP1: grid-resolution convergence of basin measures (exploratory).
c=0.15, h in {0.2,0.4,0.6}, gen-alpha rho in {0,0.4,1}, n in {100,200,400}.
Entropy is evaluated at FIXED phase-space box widths (0.2 and 0.4)."""
import numpy as np

from ..metrics import boundary_distance, entropy, mismatch
from ..metrics_ext import net_basin_transfer, signed_boundary_shift
from ..model import classify
from ..reference import rk4
from ..runmeta import experiment
from ..stats import proportion_ci
from ._common import Parts, alpha_state, family_parameters, grid, out_name, quick

DAMPING, T, T_CHECK, REF_H = 0.15, 60.0, 120.0, 0.005


def config():
    if quick():
        return {"n": [24, 48], "h": [0.4], "rho": [0.0, 0.4, 1.0]}
    return {"n": [100, 200, 400], "h": [0.2, 0.4, 0.6], "rho": [0.0, 0.4, 1.0]}


def main():
    cfg = dict(config(), damping=DAMPING, T=T, T_check=T_CHECK, reference_h=REF_H,
               box=[-2, 2], entropy_box_widths=[0.2, 0.4])
    with experiment(out_name("grid_convergence"), cfg) as run:
        parts = Parts(run["out_dir"])
        rows, checks = [], []
        for n in cfg["n"]:
            X, V = grid(n)
            spacing = (4/n, 4/n)
            key = f"ref_n{n}"
            if parts.has(key):
                meta, arr = parts.load(key)
                ref = arr["labels"]
            else:
                ref = classify(*rk4(X, V, T, REF_H, DAMPING))
                ref_h2 = classify(*rk4(X, V, T, REF_H/2, DAMPING))
                ref_t2 = classify(*rk4(X, V, T_CHECK, REF_H, DAMPING))
                meta = {"n": n, "step_check_changes": int((ref != ref_h2).sum()),
                        "time_check_changes": int((ref != ref_t2).sum()),
                        "unresolved": int((ref == 0).sum())}
                parts.save(key, meta, labels=ref)
            checks.append(meta)
            for h in cfg["h"]:
                for rho in cfg["rho"]:
                    key = f"n{n}_h{h}_rho{rho}"
                    if parts.has(key):
                        rows.append(parts.load(key)[0]); continue
                    p, _ = family_parameters("generalized-alpha", rho)
                    lab = classify(*alpha_state(X, V, T, h, DAMPING, p))
                    lab2 = classify(*alpha_state(X, V, T_CHECK, h, DAMPING, p))
                    mm = mismatch(ref, lab)
                    k = int(round(mm["fraction"]*mm["valid"])) if mm["fraction"] is not None else 0
                    row = {"n": n, "h": h, "rho": rho, "mismatch": mm["fraction"],
                           "mismatch_ci95": proportion_ci(k, mm["valid"]),
                           "valid": mm["valid"], "excluded": mm["excluded"],
                           "candidate_time_check_changes": int((lab != lab2).sum()),
                           "net_transfer": net_basin_transfer(ref, lab)["net"]}
                    if not np.any(ref == 0) and not np.any(lab == 0):
                        bd = boundary_distance(ref, lab, spacing)
                        sb = signed_boundary_shift(ref, lab, spacing)
                        row.update(boundary_mean=bd["mean_symmetric"], hausdorff=bd["hausdorff"],
                                   signed_shift=sb["mean_signed"])
                    for w in cfg["entropy_box_widths"]:
                        cells = int(round(w*n/4))
                        if cells >= 2 and n % cells == 0:
                            e = entropy(lab, cells)
                            row[f"Sb_w{w}"], row[f"Sbb_w{w}"] = e["Sb"], e["Sbb"]
                    parts.save(key, row)
                    rows.append(row)
                    print(f"n={n} h={h} rho={rho} mismatch={row['mismatch']:.4f}", flush=True)
        # convergence summary: change of mismatch between successive n
        conv = []
        for h in cfg["h"]:
            for rho in cfg["rho"]:
                seq = [r for r in rows if r["h"] == h and r["rho"] == rho]
                seq.sort(key=lambda r: r["n"])
                conv.append({"h": h, "rho": rho,
                             "mismatch_by_n": {r["n"]: r["mismatch"] for r in seq},
                             "last_change_pct_points": (100*abs(seq[-1]["mismatch"] - seq[-2]["mismatch"])
                                                        if len(seq) > 1 else None)})
        run["results"] = {"reference_checks": checks, "cells": rows, "convergence": conv,
                          "gate_rule_opus_WP1": "change n=200->400 <= 0.5 percentage points"}
    return run


if __name__ == "__main__":
    main()
