"""Plot the saved exploratory amplification curves."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "linear_exploratory"
    with (out / "curves.csv").open(newline="", encoding="utf-8") as file:
        rows = list(csv.DictReader(file))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    colors = {"generalized-alpha": "C0", "HHT": "C1", "WBZ": "C2", "Newmark": "C3"}
    for column, damped in enumerate((False, True)):
        for family, color in colors.items():
            for rho, linestyle in ((0.6, "-"), (1.0, "--")):
                group = [r for r in rows if r["family"] == family
                         and float(r["rho_target"]) == rho
                         and (float(r["physical_xi"]) > 0) == damped]
                omega = np.array([float(r["omega_h"]) for r in group])
                for row_index, key in enumerate(("xi_bar", "period_error")):
                    values = np.array([float(r[key]) for r in group])
                    axes[row_index, column].plot(
                        omega, values, color=color, linestyle=linestyle,
                        label=f"{family}, rho={rho:g}")
        axes[0, column].set_title("c=0" if not damped else "c=0.15, well frequency sqrt(2)")
        axes[1, column].set_xlabel("Omega h")
    axes[0, 0].set_ylabel("Total modal damping ratio")
    axes[1, 0].set_ylabel("Relative natural-period error")
    axes[0, 1].legend(fontsize=7)
    fig.suptitle("Exploratory linear oscillator curves")
    fig.tight_layout()
    fig.savefig(out / "curves.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
