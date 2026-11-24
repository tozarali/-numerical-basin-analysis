"""WP5 step 2 (exploratory): basin mismatch for the FE buckled clamped beam.
IC plane analogous to the SDOF (x, v): midspan deflection w0 (other DOFs
relaxed to the constrained equilibrium) and a velocity w0' along the
buckling shape phi = (q_up - q_flat)/w_m, both scaled by w_m.
Physical damping: mass-proportional, xi = 0.15/(2 sqrt 2) at the well
frequency (same modal ratio as the SDOF c=0.15). Step sizes are given as
Omega h with Omega the lowest well frequency (SDOF h=0.2,0.4,0.6 <-> 0.28,
0.57,0.85). Reference: explicit RK4 with h below the stability limit of the
highest mode; checked with h/2 and 2T. Barrier: flat compressed state
(verify with fe_arch_validation before interpreting)."""
import numpy as np
from scipy.linalg import eigh

from ..fe_arch import Arch
from ..mdof_genalpha import integrate_to_masked
from ..metrics import mismatch
from ..reference_mdof import rk4
from ..runmeta import experiment
from ._common import Parts, family_parameters, out_name, quick

FACTOR, XI = 1.6, 0.15/(2*np.sqrt(2))


def config():
    if quick():
        return {"n_el": 8, "n": 8, "omega_h": [0.57], "rho": [0.0, 1.0], "periods": 10}
    return {"n_el": 16, "n": 20, "omega_h": [0.28, 0.57, 0.85],
            "rho": [0.0, 0.2, 0.4, 0.6, 0.8, 1.0], "periods": 15}


def setup(n_el, n):
    a = Arch(n_el=n_el)
    a.end_shortening = FACTOR*a.euler_shortening
    up = a.buckled_states()[1][0]
    flat = a.flat_state()[0]
    wm = float(up[a.mid_index])
    M = a.lumped_mass()
    _, K = a.force_tangent(up)
    w_all = np.sqrt(np.abs(eigh(K, M, eigvals_only=True)))
    a.rayleigh_mass = 2*XI*w_all[0]
    s = -2 + (np.arange(n) + 0.5)*4/n
    Q, lam = a.constrained_midspan(s*wm)
    phi = (up - flat)/wm
    S, SD = np.meshgrid(np.arange(n), s)        # rows: velocity, cols: deflection
    q0 = Q[S.ravel()]
    v0 = (SD.ravel()*wm*w_all[0])[:, None]*phi   # velocity scaled like v = omega*x
    return a, q0, v0, float(a.potential(flat)), w_all, lam


def main():
    cfg = config()
    with experiment(out_name("fe_arch_basins"), dict(cfg, shortening_factor=FACTOR, xi=XI)) as run:
        parts = Parts(run["out_dir"])
        a, q0, v0, barrier, w_all, lam = setup(cfg["n_el"], cfg["n"])
        M, C = a.lumped_mass(), a.damping_matrix()
        md = a.mass_diag()
        T = cfg["periods"]*2*np.pi/w_all[0]
        h_ref = 0.5*2.785/w_all[-1]                      # half the RK4 stability limit
        h_ref = T/np.ceil(T/h_ref)
        if parts.has("ref"):
            meta, arr = parts.load("ref"); ref = arr["labels"]
        else:
            ref = a.classify(*rk4(q0, v0, T, h_ref, md, C, a), barrier)
            ref_h2 = a.classify(*rk4(q0, v0, T, h_ref/2, md, C, a), barrier)
            ref_t2 = a.classify(*rk4(q0, v0, 2*T, h_ref, md, C, a), barrier)
            meta = {"T": T, "h_ref": h_ref, "step_check_changes": int((ref != ref_h2).sum()),
                    "time_check_changes": int((ref != ref_t2).sum()),
                    "unresolved": int((ref == 0).sum()), "label_counts": np.bincount(ref + 1, minlength=3).tolist()}
            parts.save("ref", meta, labels=ref)
        print("reference:", meta, flush=True)
        rows = []
        n = cfg["n"]
        for oh in cfg["omega_h"]:
            h = oh/w_all[0]
            for rho in cfg["rho"]:
                key = f"oh{oh}_rho{rho}"
                if parts.has(key):
                    rows.append(parts.load(key)[0]); continue
                p, _ = family_parameters("generalized-alpha", rho)
                qT, vT, nfail = integrate_to_masked(q0, v0, T, h, M, C, a, p)
                lab = a.classify(qT, vT, barrier)
                mm = mismatch(ref.reshape(n, n), lab.reshape(n, n))
                row = {"omega_h": oh, "h": h, "rho": rho, "mismatch": mm["fraction"],
                       "valid": mm["valid"], "excluded": mm["excluded"], "newton_failures": nfail}
                parts.save(key, row, labels=lab)
                rows.append(row)
                print(f"Omega h={oh} rho={rho}: mismatch={row['mismatch']}, fails={nfail}", flush=True)
        run["results"] = {"reference": meta, "cells": rows, "barrier": barrier,
                          "frequencies_lowest3": w_all[:3].tolist(), "omega_max": float(w_all[-1]),
                          "constrained_ic_unstable_points": int((lam < 0).sum())}
    return run


if __name__ == "__main__":
    main()
