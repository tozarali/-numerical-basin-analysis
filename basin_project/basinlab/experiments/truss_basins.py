"""Exploratory 1-DOF and 2-DOF nonlinear von Mises truss basin experiments."""
import csv
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy

from ..metrics import boundary_distance, mismatch
from ..truss import Truss, dop853, generalized_alpha, rk4


def grid(n, horizontal):
    y = -1.2 + (np.arange(n)+0.5)*2.4/n
    vy = -1.0 + (np.arange(n)+0.5)*2/n
    Y, V = np.meshgrid(y, vy)
    q = np.stack((np.full_like(Y, horizontal), Y), axis=-1)
    v = np.stack((np.zeros_like(V), V), axis=-1)
    return q, v


def coupling_diagnostic(model):
    """Compare nonlinear vertical forces on reference 2-DOF trajectories."""
    q, v = grid(20, .6)
    sum_delta, sum_baseline, count = 0.0, 0.0, 0
    for index in range(41):
        coupled, _ = model.force_tangent(q)
        constrained_q = q.copy()
        constrained_q[..., 0] = 0
        uncoupled, _ = model.force_tangent(constrained_q)
        delta = coupled[..., 1] - uncoupled[..., 1]
        sum_delta += float(np.sum(delta**2))
        sum_baseline += float(np.sum(uncoupled[..., 1]**2))
        count += delta.size
        if index != 40:
            q, v = rk4(model, q, v, .5, .01)
    return {"sample_grid_n": 20, "times": "0:0.5:20",
            "samples": count, "rms_force_change": float(np.sqrt(sum_delta/count)),
            "rms_constrained_force": float(np.sqrt(sum_baseline/count)),
            "relative_rms_force_change": float(np.sqrt(sum_delta/sum_baseline))}


def reference_samples(model, q, v, reference_q, reference_v, constrained):
    flatq, flatv = q.reshape(-1, 2), v.reshape(-1, 2)
    flat_refq, flat_refv = reference_q.reshape(-1, 2), reference_v.reshape(-1, 2)
    indices = np.linspace(0, len(flatq)-1, 24, dtype=int)
    dop, tight = [], []
    for index in indices:
        dop.append(dop853(model, flatq[index], flatv[index], 60,
                          constrained=constrained))
        tight.append(dop853(model, flatq[index], flatv[index], 60,
                            rtol=1e-12, atol=1e-14, constrained=constrained))
    dop, tight = np.asarray(dop), np.asarray(tight)
    ref_labels = model.classify(flat_refq[indices], flat_refv[indices])
    dop_labels = model.classify(dop[:, 0], dop[:, 1])
    tight_labels = model.classify(tight[:, 0], tight[:, 1])
    return {
        "samples": len(indices),
        "indices": indices.tolist(),
        "RK4_vs_DOP853": mismatch(ref_labels, dop_labels),
        "DOP853_tolerance_comparison": mismatch(dop_labels, tight_labels),
        "max_state_error_RK4_vs_DOP853": float(np.max(np.abs(
            np.stack((flat_refq[indices], flat_refv[indices]), axis=1)-dop))),
        "max_state_error_DOP853_tolerances": float(np.max(np.abs(dop-tight)))}


def main():
    started = time.perf_counter()
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "truss_basins"
    out.mkdir(parents=True, exist_ok=True)
    model = Truss()
    rhos = np.round(np.arange(0, 1.001, .1), 1)
    hs = (.2, .4, .6)
    ns = (48, 96)
    ref_step, end = .005, 60.0
    summary = {"configuration": {
        "model": "two equal linear elastic axial bars, geometrically nonlinear",
        "half_span": model.half_span, "rise": model.rise,
        "axial_stiffness_EA_over_L0": model.axial_stiffness,
        "apex_mass": 1.0, "viscous_damping_per_coordinate": model.damping,
        "barrier": model.barrier,
        "grid": "cell centers; y in [-1.2,1.2], vertical velocity in [-1,1]",
        "reference": "RK4 h=0.005, compared to h=0.01 and DOP853",
        "basin_end": end, "candidate": "generalized-alpha",
        "D_window": 8.0, "D_population": "1.5 < initial E/barrier < 8",
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__},
        "coupling_diagnostic": coupling_diagnostic(model), "cases": []}
    rows = []
    base_references = {}

    for n in ns:
        for horizontal, constrained in ((0.0, True), (.6, False)):
            name = "1d_constrained" if constrained else "2d_perturbed"
            q, v = grid(n, horizontal)
            e0 = model.energy(q, v)
            selected = (e0 > 1.5*model.barrier) & (e0 < 8*model.barrier)
            qd, vd, ed = q[selected], v[selected], e0[selected]
            r_q, r_v = rk4(model, q, v, end, ref_step, constrained)
            c_q, c_v = rk4(model, q, v, end, .01, constrained)
            l_q, l_v = rk4(model, q, v, 2*end, .01, constrained)
            labels = model.classify(r_q, r_v)
            labels_coarse = model.classify(c_q, c_v)
            labels_long = model.classify(l_q, l_v)
            base_references[(n, name)] = labels
            ref_e_q, ref_e_v = rk4(model, qd, vd, 8, ref_step, constrained)
            d_ref = ed-model.energy(ref_e_q, ref_e_v)
            valid_d = np.isfinite(d_ref) & (d_ref > 1e-12)
            ref_e_78 = rk4(model, qd, vd, 7.8, ref_step, constrained)
            d_ref_78 = ed-model.energy(*ref_e_78)
            valid_d_78 = np.isfinite(d_ref_78) & (d_ref_78 > 1e-12)
            if not np.any(valid_d):
                raise RuntimeError(f"No valid energy-loss samples for {name} n={n}")
            independent = reference_samples(model, q, v, r_q, r_v, constrained) if n == 48 else None
            case = {"model": name, "horizontal_initial": horizontal, "n": n,
                    "reference_step_check": mismatch(labels, labels_coarse),
                    "reference_horizon_check": mismatch(labels, labels_long),
                    "reference_unresolved_60": int(np.count_nonzero(labels == 0)),
                    "reference_unresolved_120": int(np.count_nonzero(labels_long == 0)),
                    "D_selected": int(selected.sum()), "D_valid": int(valid_d.sum()),
                    "independent_reference": independent, "runs": []}
            arrays = {"y": q[..., 1], "vy": v[..., 1],
                      "reference": labels, "reference_h01": labels_coarse,
                      "reference_t120": labels_long}
            for h in hs:
                for rho in rhos:
                    cstate = generalized_alpha(model, q, v, end, h, float(rho), constrained)
                    candidate = model.classify(*cstate)
                    e_state = generalized_alpha(model, qd, vd, 8, h, float(rho),
                                                constrained) if h != .6 else None
                    # 8 is not a multiple of 0.6: compare the common horizon 7.8.
                    if e_state is None:
                        e_state = generalized_alpha(model, qd, vd, 7.8, h, float(rho), constrained)
                        loss_ref, valid, d_horizon = d_ref_78, valid_d_78, 7.8
                    else:
                        loss_ref, valid, d_horizon = d_ref, valid_d, 8.0
                    ratios = (ed-model.energy(*e_state))[valid]/loss_ref[valid]
                    mis = mismatch(labels, candidate)
                    dist = (boundary_distance(labels, candidate, (2/n, 2.4/n))
                            if np.all(labels != 0) and np.all(candidate != 0) else None)
                    row = {"model": name, "n": n, "h": h, "rho": float(rho),
                           "mismatch": mis["fraction"], "valid": mis["valid"],
                           "excluded": mis["excluded"],
                           "candidate_unresolved": int(np.count_nonzero(candidate == 0)),
                           "D_median": float(np.median(ratios)),
                           "D_q10": float(np.quantile(ratios, .1)),
                           "D_q90": float(np.quantile(ratios, .9)),
                           "D_horizon": d_horizon, "D_valid": int(valid.sum()),
                           "boundary_mean_symmetric": None if dist is None else dist["mean_symmetric"],
                           "boundary_hausdorff": None if dist is None else dist["hausdorff"]}
                    if rho in (0.0, .4, 1.0):
                        long_state = generalized_alpha(model, q, v, 120, h,
                                                       float(rho), constrained)
                        long_labels = model.classify(*long_state)
                        row["candidate_time_check"] = mismatch(candidate, long_labels)
                        row["candidate_unresolved_120"] = int(np.count_nonzero(long_labels == 0))
                    rows.append(row)
                    case["runs"].append(row)
                    arrays[f"h{h}_rho{rho}"] = candidate
                print(f"{name} n={n} h={h}: " + " ".join(
                    f"{r['rho']:.1f}:{100*r['mismatch']:.2f}%" for r in case["runs"]
                    if r["h"] == h), flush=True)
            np.savez_compressed(out / f"{name}_n{n}.npz", **arrays)
            summary["cases"].append(case)
        compare = mismatch(base_references[(n, "1d_constrained")],
                           base_references[(n, "2d_perturbed")])
        summary.setdefault("reference_perturbation_effect", []).append(
            {"n": n, "horizontal_offset": .6, "comparison": compare})
        print(f"n={n}: 1D vs 2D reference outcome difference={compare}", flush=True)

    csv_rows = [{k: v for k, v in row.items()
                 if k not in ("candidate_time_check", "candidate_unresolved_120")}
                for row in rows]
    with (out / "rho_sweep.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    summary["runtime_seconds"] = round(time.perf_counter()-started, 3)
    (out / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(rows)} runs; runtime {summary['runtime_seconds']} s", flush=True)


if __name__ == "__main__":
    main()
