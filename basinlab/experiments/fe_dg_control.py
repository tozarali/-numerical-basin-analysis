"""Exploratory: FE buckled clamped beam with the energy-consistent controls.
Same ICs, reference (re-used from data/fe_arch_basins if present, otherwise
recomputed) and step sizes as fe_arch_basins; adds the Gonzalez discrete
gradient (exact energy conservation for the nonlinear FE potential when
undamped) and rho-Bathe, and reruns generalized-alpha on the same rho grid."""
import numpy as np

from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ._common import ROOT, Parts, out_name, quick
from ._ext import integrate, key_of, method_list, summarize
from .fe_arch_basins import setup


def config():
    if quick():
        return {"n_el": 8, "n": 8, "omega_h": [0.57], "rho_step": 0.5, "periods": 10}
    return {"n_el": 16, "n": 20, "omega_h": [0.28, 0.57, 0.85], "rho_step": 0.2, "periods": 15}


def main():
    cfg = config()
    with experiment(out_name("fe_dg_control"), cfg) as run:
        parts = Parts(run["out_dir"])
        a, q0, v0, barrier, w_all, lam = setup(cfg["n_el"], cfg["n"])
        M, C, md = a.lumped_mass(), a.damping_matrix(), a.mass_diag()
        T = cfg["periods"]*2*np.pi/w_all[0]
        h_ref = 0.5*2.785/w_all[-1]
        h_ref = T/np.ceil(T/h_ref)
        n = cfg["n"]
        old = ROOT/"data"/out_name("fe_arch_basins")/"parts"/"ref.npz"
        if parts.has("ref"):
            meta, arr = parts.load("ref"); ref = arr["labels"]
        elif old.exists():
            ref = np.load(old)["labels"]; meta = {"reused_from": old.relative_to(ROOT).as_posix()}
            parts.save("ref", meta, labels=ref)
        else:
            ref = a.classify(*rk4(q0, v0, T, h_ref, md, C, a), barrier)
            meta = {"unresolved": int((ref == 0).sum())}
            parts.save("ref", meta, labels=ref)
        rows = []
        for oh in cfg["omega_h"]:
            h = oh/w_all[0]
            for name, kind, rho in method_list(cfg["rho_step"]):
                key = f"oh{oh}_{key_of(kind, rho)}"
                if parts.has(key):
                    rows.append(parts.load(key)[0]); continue
                qT, vT, nfail = integrate(kind, rho, q0, v0, T, h, M, C, a)
                lab = a.classify(qT, vT, barrier)
                mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
                row = {"omega_h": oh, "h": h, "family": name, "kind": kind, "rho": rho,
                       "mismatch": mm["fraction"], "valid": mm["valid"],
                       "excluded": mm["excluded"], "newton_failures": nfail}
                parts.save(key, row)
                rows.append(row)
                print(f"Omega h={oh} {kind} rho={rho}: {row['mismatch']} fails={nfail}", flush=True)
        run["results"] = {"reference": meta, "cells": rows, "summary": summarize(rows, ["omega_h"])}
    return run


if __name__ == "__main__":
    main()
