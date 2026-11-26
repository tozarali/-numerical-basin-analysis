"""Short exploratory timing probe for the proposed family comparison."""
import json
import os
import platform
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

from ..integrators import parameters_for_spectral_radius
from ..model import classify, energy
from ..reference import rk4
from .energy_loss_retest import exact_step_grid


FAMILY_CASE_COUNTS = {"generalized-alpha": 41, "HHT": 21, "WBZ": 39, "Newmark": 17}
PROBE_RHOS = (0.6, 0.8, 1.0)
TARGET_N = 32768


def smoke_population():
    q = -2 + (np.arange(24) + 0.5) * 4 / 24
    x, v = np.meshgrid(q, q)
    selected = (energy(x, v) > 0.6) & (energy(x, v) < 2)
    return x[selected], v[selected]


def _timed(function, *args):
    start = time.perf_counter()
    result = function(*args)
    return result, time.perf_counter() - start


def _differences(a, b):
    return {"different_labels": int(np.count_nonzero(a != b)),
            "unresolved_a": int(np.count_nonzero(a == 0)),
            "unresolved_b": int(np.count_nonzero(b == 0))}


def main():
    root = Path(__file__).resolve().parents[2]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    out = root / "data" / "confirmatory_preparation" / (stamp + "_" + uuid.uuid4().hex[:8])
    out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    x, v = smoke_population()
    n = len(x)
    damping, h = 0.15, 0.4
    reference = {}
    labels = {}
    try:
        for name, end, ref_h in (("T4_h005", 4.0, 0.005),
                                  ("T60_h005", 60.0, 0.005),
                                  ("T60_h0025", 60.0, 0.0025),
                                  ("T120_h005", 120.0, 0.005)):
            print(f"Reference: {name}", flush=True)
            state, duration = _timed(rk4, x, v, end, ref_h, damping)
            reference[name] = duration
            labels[name] = classify(*state)
            if name == "T4_h005":
                dref = energy(x, v) - energy(*state)
        valid = np.isfinite(dref) & (dref > 1e-12)
        if not np.any(valid):
            raise ArithmeticError("No valid reference energy losses")
        candidate = []
        for family in FAMILY_CASE_COUNTS:
            for rho in PROBE_RHOS:
                print(f"Candidate: {family}, rho={rho}", flush=True)
                p, mapped = parameters_for_spectral_radius(family, rho)
                e_state, te = _timed(exact_step_grid, x, v, 4.0, h, damping, p)
                s60, t60 = _timed(exact_step_grid, x, v, 60.0, h, damping, p)
                s120, t120 = _timed(exact_step_grid, x, v, 120.0, h, damping, p)
                ratio = (energy(x, v)[valid] - energy(*e_state)[valid]) / dref[valid]
                lab60, lab120 = classify(*s60), classify(*s120)
                candidate.append({"family": family, "rho": rho, "rho_mapping": mapped,
                                  "seconds": {"T4": te, "T60": t60, "T120": t120},
                                  "D_median": float(np.median(ratio)),
                                  "D_valid": int(np.count_nonzero(valid)),
                                  "reference_vs_T60": _differences(labels["T60_h005"], lab60),
                                  "T60_vs_T120": _differences(lab60, lab120)})
        average_by_family = {
            family: float(np.mean([sum(row["seconds"].values())
                                   for row in candidate if row["family"] == family]))
            for family in FAMILY_CASE_COUNTS}
        scale = TARGET_N / n
        estimate = scale * (sum(reference.values()) +
                            sum(FAMILY_CASE_COUNTS[f] * average_by_family[f]
                                for f in FAMILY_CASE_COUNTS))
        summary = {
            "status": "exploratory_smoke_not_confirmatory",
            "configuration": {"grid": "24x24 cell centers",
                              "energy_band": "0.6 < E0 < 2", "n_selected": n,
                              "damping": damping, "h": h,
                              "rho_samples": list(PROBE_RHOS),
                              "proposed_case_counts": FAMILY_CASE_COUNTS,
                              "target_n_sobol": TARGET_N},
            "environment": {"python": platform.python_version(),
                            "numpy": np.__version__, "scipy": scipy.__version__,
                            "conda_prefix": os.environ.get("CONDA_PREFIX")},
            "reference_seconds": reference,
            "reference_step_check": _differences(labels["T60_h005"], labels["T60_h0025"]),
            "reference_time_check": _differences(labels["T60_h005"], labels["T120_h005"]),
            "candidate_samples": candidate,
            "estimated_seconds_linear_scaling_compute_only": estimate,
            "estimate_warning": "Rough linear extrapolation from 24x24; excludes Sobol generation, "
                                "bootstrap, checkpoints, I/O, memory effects and environment changes. "
                                "Not a runtime guarantee or a substitute for a full-size benchmark.",
            "runtime_seconds": time.perf_counter() - started,
            "outcome": "completed"}
        (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n",
                                          encoding="utf-8")
        print(f"Saved: {out}", flush=True)
        print(f"Rough compute-only estimate: {estimate / 60:.1f} minutes", flush=True)
    except Exception as error:
        failure = {"status": "exploratory_smoke_not_confirmatory", "outcome": "failed",
                   "error": f"{type(error).__name__}: {error}"}
        (out / "summary.json").write_text(json.dumps(failure, indent=2) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    main()
