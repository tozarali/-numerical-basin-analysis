"""Quantify intercoordinate work on sampled 2-DOF reference trajectories."""
import json
from pathlib import Path

import numpy as np

from ..truss import Truss, rk4
from .truss_basins import grid


def diagnostic(sample_step):
    model = Truss()
    q, v = grid(20, .6)
    vertical = np.zeros_like(q)
    vertical[..., 1] = q[..., 1]
    initial_y_energy = .5*v[..., 1]**2 + model.potential(vertical)

    def powers(position, velocity):
        coupled, _ = model.force_tangent(position)
        position0 = position.copy()
        position0[..., 0] = 0
        uncoupled, _ = model.force_tangent(position0)
        return (-(coupled[..., 1]-uncoupled[..., 1])*velocity[..., 1],
                model.damping*velocity[..., 1]**2)

    previous_work, previous_diss = powers(q, v)
    work = np.zeros(q.shape[:-1])
    absolute_work = np.zeros_like(work)
    dissipated = np.zeros_like(work)
    for _ in range(round(20/sample_step)):
        q, v = rk4(model, q, v, sample_step, .005)
        power, diss = powers(q, v)
        work += .5*sample_step*(previous_work+power)
        absolute_work += .5*sample_step*(np.abs(previous_work)+np.abs(power))
        dissipated += .5*sample_step*(previous_diss+diss)
        previous_work, previous_diss = power, diss

    vertical[..., 1] = q[..., 1]
    final_y_energy = .5*v[..., 1]**2 + model.potential(vertical)
    residual = final_y_energy-initial_y_energy-work+dissipated
    valid = dissipated > 1e-10
    if not np.any(valid):
        raise RuntimeError("No trajectories have measurable vertical damping loss")
    return {"sample_step": sample_step, "RK4_step": .005,
            "absolute_work_over_vertical_damping_median":
                float(np.median(absolute_work[valid]/dissipated[valid])),
            "absolute_work_over_vertical_damping_q90":
                float(np.quantile(absolute_work[valid]/dissipated[valid], .9)),
            "absolute_net_work_over_vertical_damping_median":
                float(np.median(np.abs(work[valid])/dissipated[valid])),
            "energy_balance_max_absolute_residual": float(np.max(np.abs(residual))),
            "energy_balance_rms_residual": float(np.sqrt(np.mean(residual**2))),
            "energy_balance_max_dissipated": float(np.max(dissipated)),
            "valid": int(valid.sum()), "excluded": int((~valid).sum())}


def main():
    out = Path(__file__).resolve().parents[2] / "data" / "truss_basins"
    out.mkdir(parents=True, exist_ok=True)
    results = [diagnostic(.1), diagnostic(.05)]
    (out / "coupling_work.json").write_text(
        json.dumps({"duration": 20, "grid": "20x20",
                    "horizontal_initial": .6, "results": results},
                   indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
