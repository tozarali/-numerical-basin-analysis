"""TEMPLATE — one script per figure: python -m basinlab.figures.<name>
Reads only data/<experiment>/ files (never recomputes), writes vector PDF
to figures/. Colour-blind-safe palette, sizes for a single journal column."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
COL1, COL2 = 3.35, 6.9          # inches; check the target journal's guide
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]


def main(experiment="template_smoke"):
    d = np.load(ROOT / "data" / experiment / "errors.npz")
    s = json.loads((ROOT / "data" / experiment / "summary.json").read_text())
    fig, ax = plt.subplots(figsize=(COL1, 2.4), constrained_layout=True)
    ax.plot(d["rhos"], 100*d["errors"].mean(axis=1), "o-", color=PALETTE[0])
    lo, hi = s["results"]["optimum"]["ci95"]
    ax.axvspan(lo, hi, color=PALETTE[1], alpha=.25, lw=0, label="95% CI of $\\rho^*$")
    ax.set_xlabel(r"$\rho_\infty$"); ax.set_ylabel("basin mismatch (%)")
    ax.legend(frameon=False, fontsize=7)
    out = ROOT / "figures" / f"{experiment}.pdf"
    fig.savefig(out); print(out)


if __name__ == "__main__":
    main()
