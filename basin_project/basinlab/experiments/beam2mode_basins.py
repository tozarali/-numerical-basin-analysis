"""WP4 (exploratory): second independent structural model - two-mode
buckled beam (beam2mode). With q2(0)=0 the model reduces EXACTLY to the
SDOF double well (sanity check, reported). Coupling enters through an
initial second-mode displacement q20 > 0; its effect on the reference
basins is reported. IC plane: (q1, q1') on [-2,2]^2, q2' = 0."""
import numpy as np

from ..beam2mode import Beam2Mode
from ..mdof_genalpha import integrate_to_masked
from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ._common import Parts, family_parameters, grid, out_name, quick, rho_grid

T, REF_H = 60.0, 0.005


def config():
    if quick():
        return {"n": 16, "q20": [0.0, 0.3], "alpha": [1.0], "h": [0.4],
                "families": ["generalized-alpha"], "rho_step": 0.5}
    return {"n": 60, "q20": [0.0, 0.3], "alpha": [1.0, 3.0], "h": [0.2, 0.4, 0.6],
            "families": ["generalized-alpha", "WBZ", "Newmark"], "rho_step": 0.1}


def ics(n, q20):
    X, V = grid(n)
    q = np.column_stack((X.ravel(), np.full(X.size, q20)))
    v = np.column_stack((V.ravel(), np.zeros(X.size)))
    return q, v


def main():
    cfg = dict(config(), T=T, reference_h=REF_H, P=2.0, damping=[0.15, 0.15])
    with experiment(out_name("beam2mode_basins"), cfg) as run:
        parts = Parts(run["out_dir"])
        n, rows, refs = cfg["n"], [], {}
        for alpha in cfg["alpha"]:
            model = Beam2Mode(P=2.0, alpha=alpha, damping=(0.15, 0.15))
            M, C = np.eye(2), np.diag(model.damping)
            for q20 in cfg["q20"]:
                q0, v0 = ics(n, q20)
                rkey = f"ref_a{alpha}_q{q20}"
                if parts.has(rkey):
                    meta, arr = parts.load(rkey); ref = arr["labels"]
                else:
                    ref = model.classify(*rk4(q0, v0, T, REF_H, np.ones(2), C, model))
                    ref2 = model.classify(*rk4(q0, v0, T, REF_H/2, np.ones(2), C, model))
                    meta = {"alpha": alpha, "q20": q20, "step_check_changes": int((ref != ref2).sum()),
                            "unresolved": int((ref == 0).sum())}
                    parts.save(rkey, meta, labels=ref)
                refs[(alpha, q20)] = (meta, ref)
                for h in cfg["h"]:
                    for fam in cfg["families"]:
                        for rho in rho_grid(fam, cfg["rho_step"]):
                            key = f"a{alpha}_q{q20}_h{h}_{fam}_{rho}"
                            if parts.has(key):
                                rows.append(parts.load(key)[0]); continue
                            p, _ = family_parameters(fam, rho)
                            qT, vT, nfail = integrate_to_masked(q0, v0, T, h, M, C, model, p)
                            lab = model.classify(qT, vT)
                            mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
                            row = {"alpha": alpha, "q20": q20, "h": h, "family": fam, "rho": rho,
                                   "mismatch": mm["fraction"], "valid": mm["valid"],
                                   "excluded": mm["excluded"], "newton_failures": nfail}
                            parts.save(key, row)
                            rows.append(row)
                print(f"alpha={alpha} q20={q20} done", flush=True)
        coupling = []
        for alpha in cfg["alpha"]:
            a0, a3 = refs.get((alpha, 0.0)), refs.get((alpha, 0.3))
            if a0 and a3:
                ok = (a0[1] != 0) & (a3[1] != 0)
                coupling.append({"alpha": alpha,
                                 "reference_label_change_q20_0_vs_0.3": float(np.mean(a0[1][ok] != a3[1][ok]))})
        run["results"] = {"reference_checks": [m for m, _ in refs.values()],
                          "cells": rows, "coupling_effect": coupling,
                          "sanity": "q20=0 reduces to the SDOF double well (compare with regime_map c=0.15)"}
    return run


if __name__ == "__main__":
    main()
