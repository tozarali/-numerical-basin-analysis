"""Exploratory linear curves and descriptive fit to existing basin results."""
import csv
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from ..integrators import parameters_for_spectral_radius
from ..linear_analysis import dissipation_and_period, spectral_radius
from ..runmeta import experiment


def _fit(y, *features):
    x = np.column_stack((np.ones(len(y)), *features))
    coefficients, _, rank, singular = np.linalg.lstsq(x, y, rcond=None)
    if rank != x.shape[1]:
        raise ValueError("Regression design is rank deficient")
    residual = y - x @ coefficients
    r2 = 1 - np.dot(residual, residual) / np.dot(y-y.mean(), y-y.mean())
    return {"coefficients": coefficients.tolist(), "r_squared": float(r2),
            "condition_number": float(singular[0] / singular[-1])}


def _partial_spearman(x, y, control):
    x, y, control = map(rankdata, (x, y, control))
    design = np.column_stack((np.ones(len(x)), control))
    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    return float(np.corrcoef(xr, yr)[0, 1])


def main():
    root = Path(__file__).resolve().parents[2]
    source = root / "data" / "family_optimum_refinement" / "rho_refinement.csv"
    config = {"source": str(source.relative_to(root)), "h": 0.4, "damping": 0.15,
              "well_frequency": float(np.sqrt(2)), "curve_rho_targets": [0.6, 1.0],
              "curve_omega_h": [0.02, 3.0, 150],
              "seed": None, "population": "existing 100x100 grid, 0.6 < E0 < 2",
              "regression": "in-sample pooled OLS; not independent observations"}
    with experiment("linear_exploratory", config, status="exploratory") as run:
        out = Path(run["out_dir"])
        with source.open(newline="", encoding="utf-8") as file:
            original = list(csv.DictReader(file))
        p_cache = {}

        def params(family, rho):
            key = family, rho
            if key not in p_cache:
                p_cache[key] = parameters_for_spectral_radius(family, rho)
            return p_cache[key]

        omega = np.sqrt(2)
        xi = config["damping"] / (2 * omega)
        rows = []
        for row in original:
            family, rho = row["family"], float(row["rho"])
            p, _ = params(family, rho)
            linear = dissipation_and_period(p, omega * config["h"], xi)
            rows.append({"family": family, "rho": rho,
                         "D_median": float(row["D_median"]),
                         "abs_D_error": abs(float(row["D_median"]) - 1),
                         "period_error": linear["period_error"],
                         "xi_bar": linear["xi_bar"],
                         "mismatch_subset": float(row["mismatch_subset"])})

        with (out / "regression_points.csv").open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        x = np.array([row["abs_D_error"] for row in rows])
        z = np.abs(np.array([row["period_error"] for row in rows]))
        y = np.array([row["mismatch_subset"] for row in rows])
        results = {
            "n_rows": len(rows),
            "n_by_family": {family: sum(row["family"] == family for row in rows)
                            for family in sorted({row["family"] for row in rows})},
            "models": {"D_only": _fit(y, x), "period_only": _fit(y, z),
                       "D_and_period": _fit(y, x, z),
                       "D_period_interaction": _fit(y, x, z, x*z)},
            "partial_spearman_mismatch_D_given_period": _partial_spearman(x, y, z),
            "partial_spearman_mismatch_period_given_D": _partial_spearman(z, y, x),
            "limitations": "In-sample, one c and h, shared ICs and correlated rho points; "
                           "no uncertainty, causal or generalization claim."}

        curves = []
        for family in ("generalized-alpha", "HHT", "WBZ", "Newmark"):
            for rho in config["curve_rho_targets"]:
                p, measured = params(family, rho)
                for damping_ratio in (0.0, xi):
                    for wh in np.linspace(0.02, 3.0, 150):
                        values = dissipation_and_period(p, float(wh), damping_ratio)
                        curves.append({"family": family, "rho_target": rho,
                                       "rho_parameter_mapping": measured,
                                       "rho_at_1000": spectral_radius(p, 1000.0, 0.0),
                                       "omega_h": float(wh), "physical_xi": damping_ratio,
                                       "xi_bar": values["xi_bar"],
                                       "period_error": values["period_error"],
                                       "spectral_radius": spectral_radius(p, float(wh), damping_ratio)})
        with (out / "curves.csv").open("w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=list(curves[0]))
            writer.writeheader()
            writer.writerows(curves)
        run["results"] = results
        print(results, flush=True)


if __name__ == "__main__":
    main()
