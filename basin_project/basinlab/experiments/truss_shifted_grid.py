"""Repeat the 2-DOF h=0.6 case on a spatially shifted IC grid."""
import json
import time
from pathlib import Path

import numpy as np

from ..metrics import mismatch
from ..truss import Truss, generalized_alpha, rk4
from .energy_loss_retest import interpolate_crossing


def main():
    started = time.perf_counter()
    out = Path(__file__).resolve().parents[2] / "data" / "truss_basins"
    out.mkdir(parents=True, exist_ok=True)
    model = Truss()
    n, h, end, energy_end = 80, .6, 60.0, 7.8
    dy, dv = 2.4/n, 2/n
    y = -1.2 + (np.arange(n)+.85)*dy
    vy = -1 + (np.arange(n)+.30)*dv
    Y, V = np.meshgrid(y, vy)
    q = np.stack((np.full_like(Y, .6), Y), axis=-1)
    v = np.stack((np.zeros_like(V), V), axis=-1)
    e0 = model.energy(q, v)
    select = (e0 > 1.5*model.barrier) & (e0 < 8*model.barrier)
    qs, vs, es = q[select], v[select], e0[select]
    ref = model.classify(*rk4(model, q, v, end, .005))
    coarse = model.classify(*rk4(model, q, v, end, .01))
    ref_e = rk4(model, qs, vs, energy_end, .005)
    delta_ref = es-model.energy(*ref_e)
    valid = np.isfinite(delta_ref) & (delta_ref > 1e-12)
    rhos = (.2, .25, .3, .35, .4, .45, .5, .6, .7, .8, .9, 1.0)
    rows = []
    for rho in rhos:
        candidate = model.classify(*generalized_alpha(model, q, v, end, h, rho))
        est = generalized_alpha(model, qs, vs, energy_end, h, rho)
        ratio = (es-model.energy(*est))[valid]/delta_ref[valid]
        rows.append({"rho": rho, "D_median": float(np.median(ratio)),
                     "full_mismatch": mismatch(ref, candidate),
                     "subset_mismatch": mismatch(ref[select], candidate[select])})
        print(f"rho={rho:.2f} D={rows[-1]['D_median']:.5f} "
              f"full={rows[-1]['full_mismatch']['fraction']:.5f} "
              f"subset={rows[-1]['subset_mismatch']['fraction']:.5f}", flush=True)
    crossing = interpolate_crossing(np.asarray(rhos),
                                     np.asarray([r["D_median"] for r in rows]))
    best = min(rows, key=lambda row: row["subset_mismatch"]["fraction"])
    result = {"configuration": {"n": n, "grid_offsets": [.85, .30],
                                 "horizontal_initial": .6, "h": h,
                                 "end": end, "D_end": energy_end,
                                 "D_population": "1.5 < initial E/barrier < 8"},
              "reference_h_comparison": mismatch(ref, coarse),
              "D_crossing_rho": crossing, "subset_best_rho": best["rho"],
              "rows": rows, "runtime_seconds": round(time.perf_counter()-started, 3)}
    np.savez_compressed(out / "2d_shifted_n80.npz", y=Y, vy=V, reference=ref)
    (out / "shifted_grid.json").write_text(json.dumps(result, indent=2, allow_nan=False)
                                            + "\n", encoding="utf-8")
    print(f"D*={crossing}, subset best={best['rho']}, "
          f"runtime {result['runtime_seconds']} s", flush=True)


if __name__ == "__main__":
    main()
