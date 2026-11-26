"""Exploratory: tilted (asymmetric) double well V = (x^2-1)^2/4 - eps x, with
all integrator families on the same IC grid: tuned generalized-alpha rho,
rho-Bathe (Noh & Bathe 2019), trapezoid (= gen-alpha rho=1) and the
Gonzalez discrete gradient. eps = 0 reproduces the symmetric benchmark
(with Bathe added), so this one experiment also serves as the
'bathe_regime' comparison. Reference: RK4 h=0.005, T=60 (checked with h/2)."""
import numpy as np

from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ..well import Well
from ._common import Parts, grid, out_name, quick
from ._ext import integrate, key_of, method_list, summarize

T, REF_H = 60.0, 0.005


def config():
    if quick():
        return {"n": 16, "eps": [0.0, 0.1], "c": [0.15], "h": [0.4], "rho_step": 0.5}
    return {"n": 64, "eps": [0.0, 0.05, 0.1], "c": [0.05, 0.15, 0.30],
            "h": [0.2, 0.4, 0.6], "rho_step": 0.1}


def main():
    cfg = config()
    with experiment(out_name("asym_well"), dict(cfg, T=T, reference_h=REF_H)) as run:
        parts = Parts(run["out_dir"])
        n = cfg["n"]
        X, V = grid(n)
        q0, v0 = X.reshape(-1, 1), V.reshape(-1, 1)
        M = np.eye(1)
        rows, refs = [], []
        for eps in cfg["eps"]:
            model = Well(eps)
            for c in cfg["c"]:
                C = np.array([[c]])
                rkey = f"ref_e{eps}_c{c}"
                if parts.has(rkey):
                    meta, arr = parts.load(rkey); ref = arr["labels"]
                else:
                    ref = model.classify(*rk4(q0, v0, T, REF_H, np.ones(1), C, model))
                    ref2 = model.classify(*rk4(q0, v0, T, REF_H/2, np.ones(1), C, model))
                    meta = {"eps": eps, "c": c, "step_check_changes": int((ref != ref2).sum()),
                            "unresolved": int((ref == 0).sum()),
                            "fraction_right_well": float(np.mean(ref[ref != 0] == 1)),
                            "saddle": model.saddle, "barrier": model.barrier}
                    parts.save(rkey, meta, labels=ref)
                refs.append(meta)
                for h in cfg["h"]:
                    for name, kind, rho in method_list(cfg["rho_step"]):
                        key = f"e{eps}_c{c}_h{h}_{key_of(kind, rho)}"
                        if parts.has(key):
                            rows.append(parts.load(key)[0]); continue
                        qT, vT, nfail = integrate(kind, rho, q0, v0, T, h, M, C, model)
                        lab = model.classify(qT, vT)
                        mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
                        row = {"eps": eps, "c": c, "h": h, "family": name, "kind": kind, "rho": rho,
                               "mismatch": mm["fraction"], "valid": mm["valid"],
                               "excluded": mm["excluded"], "newton_failures": nfail}
                        parts.save(key, row)
                        rows.append(row)
                print(f"eps={eps} c={c} done", flush=True)
        run["results"] = {"reference_checks": refs, "cells": rows,
                          "summary": summarize(rows, ["eps", "c", "h"])}
    return run


if __name__ == "__main__":
    main()
