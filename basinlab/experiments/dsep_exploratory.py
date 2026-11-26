"""Explore near-barrier losses; distinguish real cycles from fixed windows."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ..dtheory import BARRIER, _gen_alpha_to, conservative_period
from ..integrators import parameters_for_spectral_radius
from ..linear_analysis import spectral_radius
from ..model import energy, force
from ..runmeta import experiment


def reference_events(delta, damping, rtol=1e-11, atol=1e-13):
    if not np.isfinite(delta) or not np.isfinite(damping) or delta <= 0 or damping <= 0:
        raise ValueError("delta and damping must be positive")
    e0 = BARRIER * (1 + delta)
    v0 = np.sqrt(2 * (e0 - BARRIER))

    def rhs(t, y):
        return y[1], -damping*y[1] - force(y[0]), damping*y[1]**2

    def first_turn(t, y):
        return y[1]

    def barrier(t, y):
        return energy(y[0], y[1]) - BARRIER

    def positive_crossing(t, y):
        return y[0]

    first_turn.direction = -1
    barrier.direction = -1
    positive_crossing.direction = 1
    sol = solve_ivp(rhs, (0, 60), (0, v0, 0), method="DOP853",
                    rtol=rtol, atol=atol, dense_output=True,
                    events=(first_turn, barrier, positive_crossing))
    if not sol.success:
        raise RuntimeError(sol.message)
    if not len(sol.t_events[0]):
        raise RuntimeError("Reference did not reach a turning point within 60 time units")
    turn = float(sol.t_events[0][0])
    capture = float(sol.t_events[1][0]) if len(sol.t_events[1]) else None
    returns = sol.t_events[2][sol.t_events[2] > 1e-8]
    cycle = float(returns[0]) if len(returns) else None
    if cycle is not None and capture is not None and capture < cycle:
        raise RuntimeError("Reference returned across the saddle after physical capture")
    if cycle is None and capture is None:
        cycle_status = "unresolved_before_horizon"
    elif cycle is None:
        cycle_status = "no_cycle_captured_below_barrier"
    else:
        cycle_status = "completed"
    return sol, {"e0": float(e0), "v0": float(v0), "turn_time": turn,
                 "capture_time": capture, "cycle_time": cycle,
                 "cycle_status": cycle_status,
                 "final_energy_balance_residual": float(
                     e0 - energy(*sol.y[:2, -1]) - sol.y[2, -1])}


def loss_ratio(sol, e0, v0, end, h, damping, p):
    ref = sol.sol(end)
    loss_ref = e0 - energy(ref[0], ref[1])
    if loss_ref <= 1e-12:
        raise ArithmeticError("Reference loss is too small to form a ratio")
    xn, vn = _gen_alpha_to(0.0, v0, end, h, damping, p)
    return float((e0 - energy(xn, vn)) / loss_ref)


def crossings(rhos, values):
    result = []
    for i in range(len(values)-1):
        a, b = values[i]-1, values[i+1]-1
        if a*b < 0:
            result.append(float(rhos[i] - a*(rhos[i+1]-rhos[i])/(b-a)))
        elif a == 0:
            result.append(float(rhos[i]))
    if values[-1] == 1:
        result.append(float(rhos[-1]))
    return result


def main():
    root = Path(__file__).resolve().parents[2]
    ranges = {"generalized-alpha": (0.0, 1.0), "HHT": (0.51, 1.0),
              "WBZ": (0.05, 1.0), "Newmark": (0.6, 1.0)}
    config = {"damping": 0.15, "h": [0.4, 0.5, 0.6],
              "delta": [0.02, 0.05, 0.1], "seed": None,
              "rho_ranges": ranges, "rho_points_per_family": 21,
              "initial_condition": "x=0, v=+sqrt(2*(E0-Eb))",
              "cycle": "next positive-direction x=0 crossing after t=0",
              "matching": "numerical and reference energy at reference event time; "
                          "shortened final numerical step",
              "reference": "DOP853 rtol=1e-11 atol=1e-13, checked against 1e-12/1e-14",
              "horizon": 60, "n_initial_conditions_per_delta": 1}
    with experiment("dsep_exploratory", config, status="exploratory") as run:
        out = Path(run["out_dir"])
        baseline_path = root / "data" / "family_optimum_refinement" / "summary.json"
        baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
        best = {item["family"]: item for item in baseline["families"]}
        all_rows, cases, references = [], [], []
        params = {}
        for family, (low, high) in ranges.items():
            grid = np.linspace(low, high, 21)
            params[family] = [(float(rho), *parameters_for_spectral_radius(family, float(rho)))
                              for rho in grid]

        for delta in config["delta"]:
            sol, ref = reference_events(delta, config["damping"])
            tight_sol, tight = reference_events(delta, config["damping"], 1e-12, 1e-14)
            half_period = conservative_period(ref["e0"]) / 2
            ref["conservative_half_period"] = half_period
            ref["tight_turn_time_difference"] = abs(ref["turn_time"] - tight["turn_time"])
            ref["tight_capture_time_difference"] = (
                abs(ref["capture_time"] - tight["capture_time"])
                if ref["capture_time"] is not None and tight["capture_time"] is not None
                else None)
            ref["tight_energy_difference_at_half_period"] = float(
                abs(energy(*sol.sol(half_period)[:2]) - energy(*tight_sol.sol(half_period)[:2])))
            ref["delta"] = delta
            references.append(ref)
            windows = {"first_turn_matched_time": ref["turn_time"],
                       "conservative_half_period": half_period}
            if ref["cycle_time"] is not None:
                windows["actual_cycle_matched_time"] = ref["cycle_time"]
            for h in config["h"]:
                for family, family_params in params.items():
                    family_rows = []
                    for rho, p, measured in family_params:
                        row = {"delta": delta, "h": h, "family": family,
                               "rho": rho, "rho_parameter_mapping": measured,
                               "rho_at_1000": spectral_radius(p, 1000.0),
                               "cycle_status": ref["cycle_status"]}
                        for label, end in windows.items():
                            row[label] = loss_ratio(sol, ref["e0"], ref["v0"],
                                                    end, h, config["damping"], p)
                        row.setdefault("actual_cycle_matched_time", None)
                        all_rows.append(row)
                        family_rows.append(row)
                    case = {"delta": delta, "h": h, "family": family,
                            "actual_cycle_status": ref["cycle_status"]}
                    for label in windows:
                        case[label + "_D1_crossings"] = crossings(
                            [row["rho"] for row in family_rows],
                            [row[label] for row in family_rows])
                    if h == 0.4:
                        case["existing_basin_best_rho"] = best[family]["basin_best_rho"]
                        case["existing_finite_window_D_crossing"] = best[family]["D_crossing_rho"]
                    cases.append(case)
        with (out / "ratios.csv").open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=list(all_rows[0]))
            writer.writeheader()
            writer.writerows(all_rows)
        run["results"] = {"reference_events": references, "cases": cases,
                          "interpretation": "Only actual_cycle_matched_time is a cycle measure. "
                          "The other ratios are finite-window exploratory proxies, not D_sep."}
        print(json.dumps({"reference_events": references, "cases_h_0.4":
                          [case for case in cases if case["h"] == 0.4]}, indent=2))


if __name__ == "__main__":
    main()
