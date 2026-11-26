"""Exploratory: equal-cost comparison with the rho_inf-Bathe method.

rho_inf-Bathe solves two implicit sub-steps per step, so at the same step h
it costs about twice as much as generalized-alpha or the discrete gradient
(one implicit solve per step). Here every single-stage method is run with
step h/2, so that all methods perform the same number of implicit solves
per unit time as rho_inf-Bathe at step h. Same initial conditions and
references as the base experiments:

  sdof : symmetric double well (eps = 0), 64x64 grid, as in asym_well
  beam : two-mode beam, 60x60 grid, as in beam2mode_controls
  fe   : finite-element buckled beam, 20x20 grid, as in fe_dg_control

Wall-clock time per cell is recorded as a secondary cost measure."""
import time

import numpy as np

from ..beam2mode import Beam2Mode
from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ..well import Well
from ._common import ROOT, Parts, grid, out_name, quick
from ._ext import integrate, key_of

T, REF_H = 60.0, 0.005


def config():
    if quick():
        return {"models": ["sdof"], "sdof": {"n": 16, "c": [0.15], "H": [0.4]}, "rho_step": 0.5}
    return {"models": ["sdof", "beam", "fe"],
            "sdof": {"n": 64, "c": [0.05, 0.15, 0.30], "H": [0.2, 0.4, 0.6]},
            "beam": {"n": 60, "alpha": [1.0, 3.0], "q20": [0.0, 0.3], "H": [0.2, 0.4, 0.6]},
            "fe": {"n_el": 16, "n": 20, "omega_H": [0.28, 0.57, 0.85], "periods": 15},
            "rho_step": 0.1, "fe_rho_step": 0.2}


def methods(rho_step):
    """(family, kind, rho, step divisor). Bathe at H, single-stage at H/2."""
    rhos = [round(float(r), 4) for r in np.arange(0.0, 1.0 + 1e-9, rho_step)]
    out = [("rho-Bathe", "bathe", r, 1) for r in rhos]
    out += [("generalized-alpha", "ga", r, 2) for r in rhos]
    out += [("discrete-gradient", "dg", None, 2)]
    return out


def run_cells(parts, tag, group, q0, v0, Tend, H, M, C, model, ref, n, classify, rho_step):
    rows = []
    for fam, kind, rho, div in methods(rho_step):
        key = f"{tag}_H{H}_{key_of(kind, rho)}_d{div}"
        if parts.has(key):
            rows.append(parts.load(key)[0]); continue
        t0 = time.perf_counter()
        qT, vT, nfail = integrate(kind, rho, q0, v0, Tend, H/div, M, C, model)
        sec = time.perf_counter() - t0
        lab = classify(qT, vT)
        mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
        row = dict(group, H=H, h=H/div, family=fam, kind=kind, rho=rho,
                   implicit_solves_per_H=2, mismatch=mm["fraction"], valid=mm["valid"],
                   excluded=mm["excluded"], newton_failures=nfail, seconds=sec)
        parts.save(key, row)
        rows.append(row)
    return rows


def summarize(rows, keys):
    from collections import defaultdict
    g = defaultdict(list)
    for r in rows:
        g[tuple(r[k] for k in keys)].append(r)
    out = []
    for k, rs in g.items():
        def best(kind):
            c = [r for r in rs if r["kind"] == kind]
            return min(c, key=lambda r: r["mismatch"])
        bb, bg = best("bathe"), best("ga")
        at = lambda kind, rho: [r for r in rs if r["kind"] == kind and r["rho"] == rho][0]["mismatch"]
        dg = [r for r in rs if r["kind"] == "dg"][0]
        out.append(dict(zip(keys, k), best_bathe=bb["mismatch"], best_bathe_rho=bb["rho"],
                        bathe_rho0=at("bathe", 0.0), bathe_rho1=at("bathe", 1.0),
                        best_ga_half=bg["mismatch"], best_ga_half_rho=bg["rho"],
                        trapezoid_half=at("ga", 1.0), dg_half=dg["mismatch"],
                        seconds_bathe=float(np.mean([r["seconds"] for r in rs if r["kind"] == "bathe"])),
                        seconds_ga_half=float(np.mean([r["seconds"] for r in rs if r["kind"] == "ga"])),
                        seconds_dg_half=dg["seconds"]))
    return out


def main():
    cfg = config()
    with experiment(out_name("cost_matched_controls"), dict(cfg, T=T, reference_h=REF_H)) as run:
        parts = Parts(run["out_dir"])
        res = {}
        if "sdof" in cfg["models"]:
            s = cfg["sdof"]; n = s["n"]
            X, V = grid(n)
            q0, v0 = X.reshape(-1, 1), V.reshape(-1, 1)
            model, M = Well(0.0), np.eye(1)
            rows = []
            for c in s["c"]:
                C = np.array([[c]])
                rkey = f"sdof_ref_c{c}"
                if parts.has(rkey):
                    ref = parts.load(rkey)[1]["labels"]
                else:
                    ref = model.classify(*rk4(q0, v0, T, REF_H, np.ones(1), C, model))
                    parts.save(rkey, {"c": c, "unresolved": int((ref == 0).sum())}, labels=ref)
                for H in s["H"]:
                    rows += run_cells(parts, f"sdof_c{c}", {"c": c}, q0, v0, T, H, M, C, model,
                                      ref, n, model.classify, cfg["rho_step"])
                print(f"sdof c={c} done", flush=True)
            res["sdof"] = {"cells": rows, "summary": summarize(rows, ["c", "H"])}
        if "beam" in cfg["models"]:
            s = cfg["beam"]; n = s["n"]
            X, V = grid(n)
            rows = []
            for alpha in s["alpha"]:
                model = Beam2Mode(P=2.0, alpha=alpha, damping=(0.15, 0.15))
                M, C = np.eye(2), np.diag(model.damping)
                for q20 in s["q20"]:
                    q0 = np.column_stack((X.ravel(), np.full(X.size, q20)))
                    v0 = np.column_stack((V.ravel(), np.zeros(X.size)))
                    rkey = f"beam_ref_a{alpha}_q{q20}"
                    if parts.has(rkey):
                        ref = parts.load(rkey)[1]["labels"]
                    else:
                        ref = model.classify(*rk4(q0, v0, T, REF_H, np.ones(2), C, model))
                        parts.save(rkey, {"alpha": alpha, "q20": q20}, labels=ref)
                    for H in s["H"]:
                        rows += run_cells(parts, f"beam_a{alpha}_q{q20}", {"alpha": alpha, "q20": q20},
                                          q0, v0, T, H, M, C, model, ref, n, model.classify, cfg["rho_step"])
                    print(f"beam alpha={alpha} q20={q20} done", flush=True)
            res["beam"] = {"cells": rows, "summary": summarize(rows, ["alpha", "q20", "H"])}
        if "fe" in cfg["models"]:
            from .fe_arch_basins import setup
            s = cfg["fe"]; n = s["n"]
            a, q0, v0, barrier, w_all, lam = setup(s["n_el"], n)
            M, C, md = a.lumped_mass(), a.damping_matrix(), a.mass_diag()
            Tend = s["periods"]*2*np.pi/w_all[0]
            old = ROOT/"data"/out_name("fe_arch_basins")/"parts"/"ref.npz"
            if parts.has("fe_ref"):
                ref = parts.load("fe_ref")[1]["labels"]
            elif old.exists():
                ref = np.load(old)["labels"]
                parts.save("fe_ref", {"reused_from": "data/fe_arch_basins/parts/ref.npz"}, labels=ref)
            else:
                h_ref = 0.5*2.785/w_all[-1]
                h_ref = Tend/np.ceil(Tend/h_ref)
                ref = a.classify(*rk4(q0, v0, Tend, h_ref, md, C, a), barrier)
                parts.save("fe_ref", {"unresolved": int((ref == 0).sum())}, labels=ref)
            rows = []
            cls = lambda q, v: a.classify(q, v, barrier)
            for oH in s["omega_H"]:
                H = oH/w_all[0]
                rows += run_cells(parts, f"fe_oh{oH}", {"omega_H": oH}, q0, v0, Tend, H, M, C, a,
                                  ref, n, cls, cfg["fe_rho_step"])
                print(f"fe omega_H={oH} done", flush=True)
            res["fe"] = {"cells": rows, "summary": summarize(rows, ["omega_H"])}
        run["results"] = res
    return run


if __name__ == "__main__":
    main()
