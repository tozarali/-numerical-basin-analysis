"""Confirmatory four-family comparison implementing protocol/PREREG.md
(draft of 2026-11-24). LOCKED: the full run refuses to start unless
  1. protocol/PREREG.md contains the line 'Status: APPROVED' and no 'DRAFT', and
  2. that file is tracked by Git and unmodified relative to HEAD.
The file hash and commit are stored in the output. Quick mode
(BASIN_QUICK=1) runs a small *smoke* version without the lock; its output is
data/confirmatory_family_comparison_smoke and is never confirmatory.

Design choices the draft left open are implemented as PROPOSED defaults
(see PROPOSED below); they must be copied into PREREG.md before approval.
"""
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from ..model import classify, energy
from ..reference import rk4
from ..runmeta import experiment
from ..sampling import sobol_box
from ..stats import parabolic_minimum
from ._common import ROOT, Parts, alpha_state, family_parameters, quick

PREREG = ROOT / "protocol" / "PREREG.md"
DAMPING, H = 0.15, 0.4
SOBOL_SEED, SOBOL_LOG2, N_TARGET = 20261124, 17, 32768
BOOT_SEED, N_BOOT = 20261125, 2000
BASELINE_RHO = 0.6


def rho_grids():
    r = lambda a, b: [round(x, 3) for x in np.arange(a, b + 1e-9, 0.025)]
    return {"generalized-alpha": r(0.0, 1.0),
            "HHT": [0.51] + r(0.525, 1.0),
            "WBZ": r(0.05, 1.0),
            "Newmark": r(0.6, 1.0)}


PROPOSED = {
    "quadratic_fit_half_width_points": 3,          # +-0.075 on the 0.025 grid
    "duplicate_minima": "smallest rho among equal discrete minima",
    "primary_observed_optimum": "local quadratic vertex (fallback discrete argmin when the fit is flagged)",
    "family_reliable": "reference checks pass AND every cell within +-0.1 of the discrete "
                       "optimum is reliable (also reporting the strict all-cells rule)",
    "bootstrap_missing_or_multiple_crossing": "prediction error set to the worst possible "
                                              "error on that family's rho grid (conservative)",
}


def lock_status():
    text = PREREG.read_text(encoding="utf-8") if PREREG.exists() else ""
    approved = bool(re.search(r"^Status:\s*APPROVED\b", text, re.M)) and "DRAFT" not in text
    sha = hashlib.sha256(text.encode("utf-8")).hexdigest() if text else None
    git = {"tracked": False, "clean": False, "commit": None}
    try:
        repo = ROOT.parent
        rel = PREREG.relative_to(repo).as_posix()
        git["tracked"] = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=repo,
                                        capture_output=True).returncode == 0
        git["clean"] = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=repo,
                                      capture_output=True).returncode == 0
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True)
        git["commit"] = out.stdout.strip() if out.returncode == 0 else None
    except (OSError, ValueError):
        pass
    return {"approved_text": approved, "sha256": sha, **git,
            "ok": approved and git["tracked"] and git["clean"] and git["commit"] is not None}


def population(n_target):
    pts = sobol_box(2**SOBOL_LOG2, [-2, -2], [2, 2], seed=SOBOL_SEED)
    e = energy(pts[:, 0], pts[:, 1])
    idx = np.flatnonzero((e > 0.6) & (e < 2.0))
    if len(idx) < n_target:
        raise RuntimeError(f"only {len(idx)} Sobol points in band; protocol says stop")
    idx = idx[:n_target]
    return pts[idx], idx


def unique_crossing(rhos, values, target=1.0):
    v = np.asarray(values) - target
    sign_change = np.flatnonzero(v[:-1]*v[1:] < 0)
    exact = np.flatnonzero(v == 0)
    if len(sign_change) + len(exact) != 1:
        return None
    if len(exact):
        return float(rhos[exact[0]])
    i = sign_change[0]
    return float(rhos[i] + (0 - v[i])*(rhos[i+1] - rhos[i])/(v[i+1] - v[i]))


def optimum(rhos, curve):
    rhos = np.asarray(rhos)
    i = int(np.flatnonzero(curve == np.nanmin(curve))[0])
    vertex, flag = parabolic_minimum(rhos, curve, PROPOSED["quadratic_fit_half_width_points"])
    primary = vertex if flag == "ok" else float(rhos[i])
    return {"discrete": float(rhos[i]), "quadratic": vertex, "fit_flag": flag, "primary": primary}


def analyse(grids, err, resolved, ratio, dvalid, cell_ok, ref_ok, idx=None):
    """One evaluation (full sample if idx is None, else a bootstrap sample)."""
    out = {}
    for fam, (lo, hi) in grids.items():
        rh = np.asarray(fam_rhos[fam])
        e = err[lo:hi] if idx is None else err[lo:hi][:, idx]
        r = resolved[lo:hi] if idx is None else resolved[lo:hi][:, idx]
        curve = e.sum(axis=1)/np.maximum(r.sum(axis=1), 1)
        dv = dvalid if idx is None else dvalid[idx]
        rr = ratio[lo:hi] if idx is None else ratio[lo:hi][:, idx]
        dmed = np.median(rr[:, dv], axis=1)
        opt = optimum(rh, curve)
        cross = unique_crossing(rh, dmed)
        near = np.abs(rh - opt["discrete"]) <= 0.1 + 1e-9
        reliable = bool(ref_ok and np.all(cell_ok[lo:hi][near]))
        worst = max(abs(rh[0] - opt["primary"]), abs(rh[-1] - opt["primary"]))
        err_pred = abs(cross - opt["primary"]) if cross is not None else worst
        out[fam] = {"optimum": opt, "D_crossing": cross, "reliable": reliable,
                    "reliable_strict": bool(ref_ok and np.all(cell_ok[lo:hi])),
                    "prediction_error": err_pred, "crossing_missing_or_multiple": cross is None,
                    "baseline_error": abs(BASELINE_RHO - opt["primary"]),
                    "curve": curve.tolist(), "D_median": dmed.tolist()}
    return out


fam_rhos = {}


def main():
    smoke = quick()
    name = "confirmatory_family_comparison_smoke" if smoke else "confirmatory_family_comparison"
    lock = lock_status()
    if not smoke and not lock["ok"]:
        raise SystemExit(
            "LOCKED: protocol/PREREG.md must contain 'Status: APPROVED', no 'DRAFT', and be "
            f"committed and unmodified in Git. Current: {json.dumps(lock)}")
    n_target = 512 if smoke else N_TARGET
    n_boot = 100 if smoke else N_BOOT
    grids = rho_grids()
    if smoke:
        grids = {f: g[::8] for f, g in grids.items()}
    fam_rhos.clear(); fam_rhos.update(grids)
    cfg = {"damping": DAMPING, "h": H, "sobol_seed": SOBOL_SEED, "sobol_points": 2**SOBOL_LOG2,
           "n_target": n_target, "energy_band": [0.6, 2.0], "rho_grids": grids,
           "reference": {"method": "RK4", "h": 0.005, "T": 60, "checks": ["h=0.0025", "T=120"]},
           "D": {"T_D": 4.0, "statistic": "median over valid ICs"},
           "bootstrap": {"seed": BOOT_SEED, "replicates": n_boot}, "baseline_rho": BASELINE_RHO,
           "proposed_resolutions": PROPOSED, "prereg_lock": lock, "smoke": smoke}
    with experiment(name, cfg, status="exploratory" if smoke else "confirmatory") as run:
        parts = Parts(run["out_dir"])
        pts, sel = population(n_target)
        x, v = pts[:, 0], pts[:, 1]
        if parts.has("reference"):
            rmeta, ra = parts.load("reference")
            ref, dref = ra["labels"], ra["dref"]
        else:
            ref = classify(*rk4(x, v, 60.0, 0.005, DAMPING))
            ref_h = classify(*rk4(x, v, 60.0, 0.0025, DAMPING))
            ref_t = classify(*rk4(x, v, 120.0, 0.005, DAMPING))
            dref = energy(x, v) - energy(*rk4(x, v, 4.0, 0.005, DAMPING))
            rmeta = {"step_changes": int((ref != ref_h).sum()), "time_changes": int((ref != ref_t).sum()),
                     "unresolved": int((ref == 0).sum())}
            parts.save("reference", rmeta, labels=ref, dref=dref, selection=sel)
        ref_ok = rmeta["step_changes"] == 0 and rmeta["time_changes"] == 0 and rmeta["unresolved"] == 0
        dvalid = np.isfinite(dref) & (dref > 1e-12)
        all_cases = [(f, r) for f, g in grids.items() for r in g]
        err = np.zeros((len(all_cases), len(x)), bool)
        resolved = np.zeros_like(err)
        ratio = np.full((len(all_cases), len(x)), np.nan, np.float32)
        cell_ok = np.zeros(len(all_cases), bool)
        cells = []
        for i, (fam, rho) in enumerate(all_cases):
            key = f"{fam}_{rho}"
            if not parts.has(key):
                p, measured = family_parameters(fam, rho)
                lab60 = classify(*alpha_state(x, v, 60.0, H, DAMPING, p))
                lab120 = classify(*alpha_state(x, v, 120.0, H, DAMPING, p))
                e4 = energy(*alpha_state(x, v, 4.0, H, DAMPING, p))
                rr = ((energy(x, v) - e4)/np.where(dvalid, dref, 1)).astype(np.float32)
                ok = (ref != 0) & (lab60 != 0)
                mism = float(np.mean(ref[ok] != lab60[ok])) if ok.any() else None
                meta = {"family": fam, "rho": rho, "rho_measured": measured, "mismatch": mism,
                        "unresolved": int((lab60 == 0).sum()),
                        "time_check_changes": int((lab60 != lab120).sum()),
                        "D_median": float(np.median(rr[dvalid]))}
                meta["reliable"] = (meta["unresolved"] == 0 and meta["time_check_changes"] == 0
                                    and mism is not None and mism < 0.30)
                parts.save(key, meta, lab60=lab60, ratio=rr)
                print(f"{fam} rho={rho}: mismatch={mism}", flush=True)
            meta, arr = parts.load(key)
            lab60 = arr["lab60"]
            resolved[i] = (ref != 0) & (lab60 != 0)
            err[i] = resolved[i] & (ref != lab60)
            ratio[i] = arr["ratio"]
            cell_ok[i] = meta["reliable"]
            cells.append(meta)
        bounds, k = {}, 0
        for fam, g in grids.items():
            bounds[fam] = (k, k + len(g)); k += len(g)
        point = analyse(bounds, err, resolved, ratio, dvalid, cell_ok, ref_ok)
        rng = np.random.default_rng(BOOT_SEED)
        diffs, boot_opt = [], {f: [] for f in grids}
        for b in range(n_boot):
            idx = rng.integers(0, len(x), len(x))
            res = analyse(bounds, err, resolved, ratio, dvalid, cell_ok, ref_ok, idx)
            rel = [f for f in grids if point[f]["reliable"]]
            if rel:
                diffs.append(np.mean([res[f]["prediction_error"] - res[f]["baseline_error"] for f in rel]))
            for f in grids:
                boot_opt[f].append(res[f]["optimum"]["primary"])
        reliable = [f for f in grids if point[f]["reliable"]]
        hits = [f for f in reliable if not point[f]["crossing_missing_or_multiple"]
                and point[f]["prediction_error"] <= 0.1]
        ci = np.percentile(diffs, [2.5, 97.5]).tolist() if diffs else None
        if len(reliable) < 3:
            verdict = "inconclusive (fewer than three reliable families)"
        elif len(hits)/len(reliable) >= 0.8 and ci is not None and ci[1] < 0:
            verdict = "passed"
        else:
            verdict = "failed"
        for f in grids:
            point[f]["optimum_bootstrap_ci95"] = np.percentile(boot_opt[f], [2.5, 97.5]).tolist()
        run["results"] = {"reference_checks": rmeta, "reference_ok": ref_ok, "families": point,
                          "reliable_families": reliable, "hits_within_0.1": hits,
                          "mean_error_difference_vs_baseline_ci95": ci,
                          "gate_verdict": verdict if not smoke else f"SMOKE ONLY ({verdict})",
                          "cells": cells}
        print("Gate:", run["results"]["gate_verdict"], flush=True)
    return run


if __name__ == "__main__":
    main()
