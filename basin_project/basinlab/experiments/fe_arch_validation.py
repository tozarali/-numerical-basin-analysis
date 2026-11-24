"""WP5 step 1 (exploratory verification of the FE model, no basins):
  * Euler buckling: lowest eigenvalue of the flat compressed state changes
    sign at Delta_cr = 4 pi^2 EI L / (EA L^2)  (closed form);
  * post-buckling amplitude vs first-order theory;
  * mesh convergence (n_el = 8, 16, 32) of amplitude, barrier, frequencies;
  * index-1 saddle search: Newton from perturbed flat states; reports every
    equilibrium found with its number of negative eigenvalues, so we know
    whether the flat state is the lowest saddle (barrier for classification)."""
import numpy as np
from scipy.linalg import eigh

from ..fe_arch import Arch
from ..runmeta import experiment
from ._common import out_name, quick

SHORTENING_FACTOR = 1.6


def buckling_crossing(n_el):
    a = Arch(n_el=n_el)
    dcr = a.euler_shortening
    facs = np.linspace(0.9, 1.1, 21)
    lam = []
    for f in facs:
        a.end_shortening = f*dcr
        lam.append(float(a.flat_state()[2][0]))
    lam = np.array(lam)
    i = int(np.flatnonzero((lam[:-1] > 0) & (lam[1:] <= 0))[0])
    f_cross = facs[i] + lam[i]*(facs[i+1] - facs[i])/(lam[i] - lam[i+1])
    return float(f_cross)


def saddle_search(a, n_tries=12, seed=0):
    rng = np.random.default_rng(seed)
    flat, _, _ = a.flat_state()
    found = []
    s = a.X[:, 0]
    for t in range(n_tries):
        d = np.zeros(a.ndof)
        k = rng.integers(1, 4)
        amp = rng.uniform(-0.04, 0.04)
        d[1::3] = amp*np.sin(k*np.pi*s/a.span)        # k=2: antisymmetric guess
        q0 = flat + d[a.free] + rng.normal(scale=1e-4, size=len(a.free))
        q, conv, ev = a.equilibrium(q0, maxiter=80)
        if not conv:
            continue
        e = float(a.potential(q))
        key = (round(e, 7), int((ev < 0).sum()))
        if key not in [(round(f["energy"], 7), f["negative_eigs"]) for f in found]:
            found.append({"energy": e, "negative_eigs": int((ev < 0).sum()),
                          "w_mid": float(q[a.mid_index])})
    return sorted(found, key=lambda f: f["energy"])


def main():
    meshes = [8, 16] if quick() else [8, 16, 32]
    cfg = {"meshes": meshes, "shortening_factor": SHORTENING_FACTOR, "EA": 1e4, "EI": 1.0,
           "span": 1.0, "imperfection": 0.0}
    with experiment(out_name("fe_arch_validation"), cfg) as run:
        rows = []
        for n in meshes:
            a = Arch(n_el=n)
            a.end_shortening = SHORTENING_FACTOR*a.euler_shortening
            st = a.buckled_states()
            flat, _, evf = a.flat_state()
            up, conv_up, ev_up = st[1]
            down, conv_dn, ev_dn = st[-1]
            wm = float(up[a.mid_index])
            excess = a.end_shortening - a.euler_shortening
            M = a.lumped_mass()
            _, K = a.force_tangent(up)
            w = np.sqrt(eigh(K, M, eigvals_only=True)[:3])
            saddles = saddle_search(a, 6 if quick() else 16)
            index1 = [f for f in saddles if f["negative_eigs"] == 1]
            row = {"n_el": n, "buckling_factor_crossing": buckling_crossing(n),
                   "up_converged": conv_up, "down_converged": conv_dn,
                   "up_min_eig": float(ev_up[0]), "down_min_eig": float(ev_dn[0]),
                   "w_mid_up": wm, "w_mid_down": float(down[a.mid_index]),
                   "w_mid_first_order": float(np.sqrt(4*a.span*excess/np.pi**2)),
                   "flat_negative_eigs": int((evf < 0).sum()),
                   "barrier_flat_minus_up": float(a.potential(flat) - a.potential(up)),
                   "well_frequencies": w.tolist(),
                   "equilibria_found": saddles,
                   "lowest_index1_energy_minus_up": (float(index1[0]["energy"] - a.potential(up))
                                                     if index1 else None)}
            rows.append(row)
            print(f"n_el={n}: Pcr factor={row['buckling_factor_crossing']:.4f}, "
                  f"w_m={wm:.5f} (1st order {row['w_mid_first_order']:.5f}), "
                  f"barrier={row['barrier_flat_minus_up']:.5f}", flush=True)
        checks = {
            "euler_load_within_1pct_finest": abs(rows[-1]["buckling_factor_crossing"] - 1) < 0.01,
            "mirror_states_stable": all(r["up_min_eig"] > 0 and r["down_min_eig"] > 0 for r in rows),
            "flat_is_index1": all(r["flat_negative_eigs"] == 1 for r in rows),
            "flat_is_lowest_index1_found": all(
                r["lowest_index1_energy_minus_up"] is None or
                r["lowest_index1_energy_minus_up"] >= r["barrier_flat_minus_up"] - 1e-9 for r in rows)}
        run["results"] = {"meshes": rows, "checks": checks,
                          "open_items": ["literature dynamic benchmark", "openseespy cross-check",
                                         "rotary inertia justification", "saddle search is not exhaustive"]}
    return run


if __name__ == "__main__":
    main()
