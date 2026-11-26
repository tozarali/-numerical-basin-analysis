"""Exploratory: same IC grids/references as beam2mode_basins, with rho-Bathe
and the Gonzalez discrete gradient added as controls (generalized-alpha is
rerun at the same rho grid so all methods share one table)."""
import numpy as np

from ..beam2mode import Beam2Mode
from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ._common import Parts, grid, out_name, quick
from ._ext import integrate, key_of, method_list, summarize

T, REF_H = 60.0, 0.005


def config():
    if quick():
        return {"n": 12, "q20": [0.3], "alpha": [1.0], "h": [0.4], "rho_step": 0.5}
    return {"n": 60, "q20": [0.0, 0.3], "alpha": [1.0, 3.0], "h": [0.2, 0.4, 0.6], "rho_step": 0.1}


def main():
    cfg = config()
    with experiment(out_name("beam2mode_controls"), dict(cfg, T=T, reference_h=REF_H, P=2.0)) as run:
        parts = Parts(run["out_dir"])
        n = cfg["n"]
        X, V = grid(n)
        rows, refs = [], []
        for alpha in cfg["alpha"]:
            model = Beam2Mode(P=2.0, alpha=alpha, damping=(0.15, 0.15))
            M, C = np.eye(2), np.diag(model.damping)
            for q20 in cfg["q20"]:
                q0 = np.column_stack((X.ravel(), np.full(X.size, q20)))
                v0 = np.column_stack((V.ravel(), np.zeros(X.size)))
                rkey = f"ref_a{alpha}_q{q20}"
                if parts.has(rkey):
                    meta, arr = parts.load(rkey); ref = arr["labels"]
                else:
                    ref = model.classify(*rk4(q0, v0, T, REF_H, np.ones(2), C, model))
                    ref2 = model.classify(*rk4(q0, v0, T, REF_H/2, np.ones(2), C, model))
                    meta = {"alpha": alpha, "q20": q20, "step_check_changes": int((ref != ref2).sum()),
                            "unresolved": int((ref == 0).sum())}
                    parts.save(rkey, meta, labels=ref)
                refs.append(meta)
                for h in cfg["h"]:
                    for name, kind, rho in method_list(cfg["rho_step"]):
                        key = f"a{alpha}_q{q20}_h{h}_{key_of(kind, rho)}"
                        if parts.has(key):
                            rows.append(parts.load(key)[0]); continue
                        qT, vT, nfail = integrate(kind, rho, q0, v0, T, h, M, C, model)
                        lab = model.classify(qT, vT)
                        mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
                        row = {"alpha": alpha, "q20": q20, "h": h, "family": name, "kind": kind,
                               "rho": rho, "mismatch": mm["fraction"], "valid": mm["valid"],
                               "excluded": mm["excluded"], "newton_failures": nfail}
                        parts.save(key, row)
                        rows.append(row)
                print(f"alpha={alpha} q20={q20} done", flush=True)
        run["results"] = {"reference_checks": refs, "cells": rows,
                          "summary": summarize(rows, ["alpha", "q20", "h"])}
    return run


if __name__ == "__main__":
    main()
