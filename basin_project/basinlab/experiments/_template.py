"""TEMPLATE — copy to <name>.py. Run: python -m basinlab.experiments.<name>
Exploratory by default. Set STATUS='confirmatory' only if protocol/PREREG.md
is committed and this config matches it exactly."""
import numpy as np

from ..integrators import parameters_for_spectral_radius
from ..model import classify, energy
from ..reference import rk4
from ..runmeta import experiment
from ..sampling import energy_band, sobol_box
from ..stats import bootstrap_optimum
from .energy_loss_retest import exact_step_grid

NAME = "template_smoke"
STATUS = "exploratory"
CONFIG = {"family": "generalized-alpha", "damping": 0.15, "h": 0.4,
          "rhos": [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8], "T_basin": 60.0,
          "reference_h": 0.01, "n_sobol": 1024, "seed": 12345,
          "box": [[-2, -2], [2, 2]], "energy_band": [0.25, 2.0]}


def main(cfg=CONFIG):
    with experiment(NAME, cfg, status=STATUS) as run:
        pts = sobol_box(cfg["n_sobol"], *cfg["box"], seed=cfg["seed"])
        pts, _ = energy_band(pts, energy, *cfg["energy_band"])
        x, v = pts[:, 0], pts[:, 1]
        ref = classify(*rk4(x, v, cfg["T_basin"], cfg["reference_h"], cfg["damping"]))
        rhos = np.asarray(cfg["rhos"])
        errors = np.zeros((len(rhos), len(x)), bool)
        unresolved = []
        for i, rho in enumerate(rhos):
            p, _ = parameters_for_spectral_radius(cfg["family"], float(rho))
            lab = classify(*exact_step_grid(x, v, cfg["T_basin"], cfg["h"], cfg["damping"], p))
            ok = (ref != 0) & (lab != 0)
            errors[i] = ok & (ref != lab)
            unresolved.append(int((~ok).sum()))
        np.savez_compressed(f"{run['out_dir']}/errors.npz", rhos=rhos, errors=errors, ics=pts)
        run["results"] = {"n_ic": int(len(x)), "unresolved_per_rho": unresolved,
                          "mismatch": errors.mean(axis=1).tolist(),
                          "optimum": bootstrap_optimum(rhos, errors, n_boot=300)}
    return run


if __name__ == "__main__":
    print(main()["results"])
