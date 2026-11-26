"""Generate all manuscript figures from the result files in data/.

Usage (from the repository root):  python paper/make_figures.py
Environment variable BASIN_DATA overrides the data directory.
Figures are written to paper/figures/ as vector PDF."""
import collections
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("BASIN_DATA", HERE.parent / "data"))
OUT = HERE / "figures_new"
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": True,
    "axes.spines.right": True, "figure.dpi": 150, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "lines.linewidth": 1.2, "lines.markersize": 3.5,
})
W1, W2 = 3.5, 7.2           # single and double column width (inch)
COL = {"generalized-alpha": "#2a6fdb", "HHT": "#e0a100", "WBZ": "#1b9e77",
       "Newmark": "#d95f02", "discrete-gradient": "#c51b7d", "rho-Bathe": "#6a3d9a",
       "BDF2": "#555555", "trapezoid": "#333333"}
NAME = {"generalized-alpha": r"generalized-$\alpha$", "HHT": "HHT", "WBZ": "WBZ",
        "Newmark": "Newmark", "discrete-gradient": "discrete gradient",
        "rho-Bathe": r"$\rho_\infty$-Bathe"}
BASIN_CMAP = ListedColormap(["#1f5f99", "#f2f2f2", "#f0d33c"])   # -1, 0, +1


def J(p):
    return json.loads((DATA / p).read_text(encoding="utf-8"))


def tag(ax, s):
    ax.text(-0.16, 1.0, s, transform=ax.transAxes, fontweight="bold", va="bottom", ha="left",
            fontsize=9)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name)


# ----------------------------------------------------------------- regime map
RM = J("regime_map/summary.json")["results"]
CELLS = RM["cells"]
OPT = RM["family_optima"]


def regime_part(c, h):
    meta = json.loads((DATA / "regime_map/parts" / f"c{c}_h{h}.json").read_text())
    lab = np.load(DATA / "regime_map/parts" / f"c{c}_h{h}.npz")["labels"]
    ref = np.load(DATA / "regime_map/parts" / f"ref_c{c}.npz")["labels"]
    idx = {(f, r): i for i, (f, r) in enumerate(meta["cases"])}
    return meta, lab, idx, ref


def fig_basin_maps():
    c, h = 0.15, 0.4
    meta, lab, idx, ref = regime_part(c, h)
    opt = {o["family"]: o for o in OPT if o["c"] == c and o["h"] == h}
    rg = opt["generalized-alpha"]["rho_best"]
    panels = [("Reference (RK4)", ref, None),
              (r"trapezoidal ($\rho_\infty=1$)", lab[idx[("generalized-alpha", 1.0)]], ref),
              (rf"generalized-$\alpha$, $\rho_\infty={rg:.2f}$", lab[idx[("generalized-alpha", rg)]], ref),
              ("discrete gradient", lab[idx[("discrete-gradient", None)]], ref)]
    fig, axs = plt.subplots(1, 4, figsize=(W2, 2.0), sharey=True)
    for k, (ax, (title, L, R)) in enumerate(zip(axs, panels)):
        ax.imshow(L, origin="lower", extent=[-2, 2, -2, 2], cmap=BASIN_CMAP, vmin=-1, vmax=1,
                  interpolation="nearest")
        if R is not None:
            bad = np.ma.masked_where(L == R, np.ones_like(L, float))
            ax.imshow(bad, origin="lower", extent=[-2, 2, -2, 2], cmap=ListedColormap(["k"]),
                      alpha=0.55, interpolation="nearest")
            v = (L != 0) & (R != 0)
            title += f"\n{100*np.mean(L[v] != R[v]):.1f}% mismatch"
        else:
            title += "\n"
        ax.set_title(title)
        ax.set_xlabel(r"$x_0$")
        ax.set_aspect("equal")
        tag(ax, "abcd"[k])
    axs[0].set_ylabel(r"$\dot x_0$")
    save(fig, "fig01_basin_maps")


def fig_mismatch_curves():
    c = 0.15
    hs = (0.2, 0.4, 0.6)
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.2), sharey=True)
    for k, (ax, h) in enumerate(zip(axs, hs)):
        cs = [x for x in CELLS if x["c"] == c and x["h"] == h]
        for fam in ("generalized-alpha", "HHT", "WBZ", "Newmark"):
            cur = sorted((x["rho"], x["mismatch"]) for x in cs if x["family"] == fam)
            r, m = np.array(cur).T
            ax.plot(r, 100*m, "-", color=COL[fam], label=NAME[fam])
            i = int(np.argmin(m))
            ax.plot(r[i], 100*m[i], "*", color=COL[fam], ms=7, mec="k", mew=0.4)
        dg = [x for x in cs if x["family"] == "discrete-gradient"][0]["mismatch"]
        ax.axhline(100*dg, color=COL["discrete-gradient"], ls="--", label=NAME["discrete-gradient"])
        ax.axhline(30, color="0.6", lw=0.6, ls=":")
        ax.set_title(f"$c={c}$, $h={h}$")
        ax.set_xlabel(r"$\rho_\infty$")
        ax.set_ylim(0, 50)
        tag(ax, "abc"[k])
    axs[0].set_ylabel("basin mismatch (%)")
    hd, lb = axs[0].get_legend_handles_labels()
    fig.legend(hd, lb, frameon=False, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.08))
    save(fig, "fig02_mismatch_curves")


def informative():
    by = collections.defaultdict(list)
    for o in OPT:
        by[(o["c"], o["h"])].append(o)
    out = {}
    for k, v in by.items():
        ok = [o for o in v if not o["best_cell_saturated"]]
        if len(ok) == 4 and min(o["mismatch_best"] for o in ok) >= 0.01:
            out[k] = {o["family"]: o for o in ok}
    return out


def fig_optimum_map():
    inf = informative()
    cs = sorted({c for c, _ in (( o["c"], o["h"]) for o in OPT)})
    hs = sorted({o["h"] for o in OPT})
    fig, axs = plt.subplots(1, 4, figsize=(W2, 2.0), sharey=True)
    for k, (ax, fam) in enumerate(zip(axs, ("generalized-alpha", "HHT", "WBZ", "Newmark"))):
        A = np.full((len(cs), len(hs)), np.nan)
        for (c, h), d in inf.items():
            A[cs.index(c), hs.index(h)] = d[fam]["rho_best"]
        im = ax.imshow(A, origin="lower", cmap="viridis", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(hs)), [f"{h:g}" for h in hs], rotation=90)
        ax.set_yticks(range(len(cs)), [f"{c:g}" for c in cs])
        ax.set_title(NAME[fam])
        ax.set_xlabel("$h$")
        tag(ax, "abcd"[k])
    axs[0].set_ylabel("$c$")
    cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.02)
    cb.set_label(r"optimal $\rho_\infty$")
    save(fig, "fig03_optimum_map")


def fig_gain_maps():
    inf = informative()
    dg = {(x["c"], x["h"]): x["mismatch"] for x in CELLS if x["family"] == "discrete-gradient"}
    cs = sorted({o["c"] for o in OPT}); hs = sorted({o["h"] for o in OPT})
    A = np.full((len(cs), len(hs)), np.nan); B = A.copy()
    for (c, h), d in inf.items():
        best = min(o["mismatch_best"] for o in d.values())
        A[cs.index(c), hs.index(h)] = best/max(o["mismatch_rho1"] for o in d.values())
        B[cs.index(c), hs.index(h)] = dg[(c, h)]/best
    fig, axs = plt.subplots(1, 2, figsize=(W2*0.7, 2.2), sharey=True)
    for k, (ax, M, title, cmap, lim) in enumerate(zip(
            axs, (A, B), ("best tuned / trapezoidal", "discrete gradient / best tuned"),
            ("magma", "PuOr_r"), ((0, 1), (0, 2)))):
        im = ax.imshow(M, origin="lower", cmap=cmap, vmin=lim[0], vmax=lim[1], aspect="auto")
        for i in range(len(cs)):
            for j in range(len(hs)):
                if np.isfinite(M[i, j]):
                    ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                            color="w" if (k == 0 and M[i, j] < 0.5) else "k")
        ax.set_xticks(range(len(hs)), [f"{h:g}" for h in hs], rotation=90)
        ax.set_yticks(range(len(cs)), [f"{c:g}" for c in cs])
        ax.set_xlabel("$h$"); ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
        tag(ax, "ab"[k])
    axs[0].set_ylabel("$c$")
    save(fig, "fig04_gain_maps")


# ------------------------------------------------------------ confirmatory
def fig_confirmatory():
    cs = J("confirmatory_family_comparison/summary.json")
    r = cs["results"]
    grid_of = {}
    for cell in r["cells"]:
        grid_of.setdefault(cell["family"], []).append(cell["rho"])
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.5), gridspec_kw={"width_ratios": [1, 1, 1.5]})
    for fam, v in r["families"].items():
        rh = np.array(sorted(grid_of[fam]))
        m, D = np.array(v["curve"]), np.array(v["D_median"])
        n = min(len(rh), len(m))
        lab = NAME[fam] + ("" if v["reliable"] else " (unreliable)")
        axs[0].plot(rh[:n], 100*m[:n], color=COL[fam], label=lab)
        axs[0].plot(v["optimum"]["primary"], 100*np.interp(v["optimum"]["primary"], rh[:n], m[:n]),
                    "*", color=COL[fam], ms=7, mec="k", mew=0.4)
        axs[1].plot(rh[:n], D[:n], color=COL[fam])
        axs[1].plot(v["optimum"]["primary"], np.interp(v["optimum"]["primary"], rh[:n], D[:n]),
                    "*", color=COL[fam], ms=7, mec="k", mew=0.4)
        if v["D_crossing"] is not None:
            axs[1].plot(v["D_crossing"], 1.0, "v", color=COL[fam], ms=5, mec="k", mew=0.4)
    dd = r["H2_discrete_gradient"]["descriptive"]
    axs[0].axhline(100*dd["discrete_gradient_mismatch"], color=COL["discrete-gradient"], ls="--",
                   label=NAME["discrete-gradient"])
    axs[0].set_ylim(0, 40); axs[0].set_xlabel(r"$\rho_\infty$"); axs[0].set_ylabel("basin mismatch (%)")
    hd, lb = axs[0].get_legend_handles_labels()
    axs[1].axhline(1, color="0.5", lw=0.6)
    axs[1].set_xlabel(r"$\rho_\infty$"); axs[1].set_ylabel(r"median $D$")
    fm = {f: min(v["curve"]) for f, v in r["families"].items()}
    rows = [("discrete gradient", dd["discrete_gradient_mismatch"], COL["discrete-gradient"]),
            ("Newmark, best cell", fm["Newmark"], COL["Newmark"]),
            (r"best gen.-$\alpha$/HHT/WBZ", min(fm[x] for x in ("generalized-alpha", "HHT", "WBZ")),
             COL["generalized-alpha"]),
            (r"gen.-$\alpha$, $\rho_\infty=0.6$", dd["family_rho_0.6_mismatch"]["generalized-alpha"], "#7fa7e8"),
            (r"trapezoidal, $\rho_\infty=1$", dd["trapezoidal_rho_1_mismatch"], COL["trapezoid"])]
    y = np.arange(len(rows))[::-1]
    axs[2].barh(y, [100*v for _, v, _ in rows], color=[c for _, _, c in rows])
    for yi, (_, v, _) in zip(y, rows):
        axs[2].text(100*v + 0.3, yi, f"{100*v:.1f}", va="center", fontsize=6.5)
    axs[2].set_yticks(y, [a for a, _, _ in rows]); axs[2].set_xlabel("basin mismatch (%)")
    axs[2].set_xlim(0, 20)
    for k, ax in enumerate(axs):
        tag(ax, "abc"[k])
    fig.tight_layout(rect=(0, 0.09, 1, 1), w_pad=1.0)
    fig.legend(hd, lb, frameon=False, loc="lower left", ncol=5, fontsize=6.5,
               bbox_to_anchor=(0.02, 0.0))
    save(fig, "fig05_confirmatory")


# ------------------------------------------------- extended controls (Bathe, DG)
def _ext_points():
    """(model, label, best_ga, best_bathe, bathe_rho0, dg, trapezoid) per regime."""
    pts = []
    for s in J("asym_well/summary.json")["results"]["summary"]:
        pts.append(("SDOF" if s["eps"] == 0 else "SDOF tilted", s))
    for s in J("beam2mode_controls/summary.json")["results"]["summary"]:
        pts.append(("two-mode beam", s))
    for s in J("fe_dg_control/summary.json")["results"]["summary"]:
        pts.append(("FE beam", s))
    return pts


MK = {"SDOF": "o", "SDOF tilted": "s", "two-mode beam": "^", "FE beam": "D"}


def fig_controls_equal_step():
    pts = _ext_points()
    fig, axs = plt.subplots(1, 2, figsize=(W2*0.8, 2.7), sharey=True)
    for model in MK:
        P = [s for m, s in pts if m == model]
        x = np.array([s["best_ga"] for s in P])
        axs[0].scatter(100*x, 100*np.array([s["dg"] for s in P]), marker=MK[model], s=14,
                       facecolor="none", edgecolor=COL["discrete-gradient"], lw=0.8)
        axs[0].scatter(100*x, 100*np.array([s["best_bathe"] for s in P]), marker=MK[model], s=14,
                       color=COL["rho-Bathe"], lw=0)
        axs[1].scatter(100*np.array([s["trapezoid"] for s in P]), 100*np.array([s["bathe_rho0"] for s in P]),
                       marker=MK[model], s=14, color=COL["rho-Bathe"], lw=0)
        axs[1].scatter(100*np.array([s["trapezoid"] for s in P]), 100*np.array([s["dg"] for s in P]),
                       marker=MK[model], s=14, facecolor="none", edgecolor=COL["discrete-gradient"], lw=0.8)
    for ax, xl in zip(axs, (r"best tuned generalized-$\alpha$ (%)", r"trapezoidal rule, $\rho_\infty=1$ (%)")):
        ax.plot([0.3, 70], [0.3, 70], color="0.5", lw=0.6)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlim(0.3, 70); ax.set_ylim(0.3, 70)
        ax.set_xlabel(xl); ax.set_aspect("equal")
    axs[0].set_ylabel("mismatch of the control (%)")
    axs[0].set_title("vs. tuned generalized-$\\alpha$")
    axs[1].set_title("vs. trapezoidal rule")
    h = [plt.Line2D([], [], marker=MK[m], ls="", color="k", mfc="none", label=m) for m in MK]
    h += [plt.Line2D([], [], marker="o", ls="", color=COL["rho-Bathe"], label=r"$\rho_\infty$-Bathe (a: best $\rho_\infty$; b: $\rho_\infty=0$)"),
          plt.Line2D([], [], marker="o", ls="", mfc="none", color=COL["discrete-gradient"], label="discrete gradient")]
    fig.legend(handles=h, frameon=False, loc="lower center", ncol=3, fontsize=6, bbox_to_anchor=(0.5, -0.2))
    tag(axs[0], "a"); tag(axs[1], "b")
    save(fig, "fig06_controls_equal_step")


def fig_cost_matched():
    p = DATA / "cost_matched_controls/summary.json"
    if not p.exists():
        print("skip fig07 (no cost_matched_controls)"); return
    res = J("cost_matched_controls/summary.json")["results"]
    lab = {"sdof": "SDOF", "beam": "two-mode beam", "fe": "FE beam"}
    fig, axs = plt.subplots(1, 2, figsize=(W2*0.8, 2.7), sharey=True)
    for key, model in lab.items():
        if key not in res:
            continue
        S = res[key]["summary"]
        b = 100*np.array([s["best_bathe"] for s in S])
        b0 = 100*np.array([s["bathe_rho0"] for s in S])
        axs[0].scatter(100*np.array([s["best_ga_half"] for s in S]), b, marker=MK[model], s=14,
                       color=COL["rho-Bathe"], lw=0)
        axs[0].scatter(100*np.array([s["best_ga_half"] for s in S]), b0, marker=MK[model], s=14,
                       facecolor="none", edgecolor=COL["rho-Bathe"], lw=0.8)
        axs[1].scatter(100*np.array([s["dg_half"] for s in S]), b, marker=MK[model], s=14,
                       color=COL["rho-Bathe"], lw=0)
        axs[1].scatter(100*np.array([s["dg_half"] for s in S]), b0, marker=MK[model], s=14,
                       facecolor="none", edgecolor=COL["rho-Bathe"], lw=0.8)
    for ax, xl in zip(axs, (r"best generalized-$\alpha$ at $h/2$ (%)", r"discrete gradient at $h/2$ (%)")):
        ax.plot([0.1, 70], [0.1, 70], color="0.5", lw=0.6)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 70); ax.set_ylim(0.1, 70)
        ax.set_xlabel(xl); ax.set_aspect("equal")
    axs[0].set_ylabel(r"$\rho_\infty$-Bathe at step $h$ (%)")
    h = [plt.Line2D([], [], marker=MK[m], ls="", color="k", mfc="none", label=m)
         for m in ("SDOF", "two-mode beam", "FE beam")]
    h += [plt.Line2D([], [], marker="o", ls="", color=COL["rho-Bathe"], label=r"best $\rho_\infty$"),
          plt.Line2D([], [], marker="o", ls="", mfc="none", color=COL["rho-Bathe"], label=r"$\rho_\infty=0$")]
    fig.legend(handles=h, frameon=False, loc="lower center", ncol=5, fontsize=6, bbox_to_anchor=(0.5, -0.17))
    tag(axs[0], "a"); tag(axs[1], "b")
    save(fig, "fig07_cost_matched")


def fig_structures():
    fe = J("fe_dg_control/summary.json")["results"]["cells"]
    bm = J("beam2mode_controls/summary.json")["results"]["cells"]
    fig, axs = plt.subplots(1, 4, figsize=(W2, 2.1), sharey=False)
    for k, oh in enumerate((0.28, 0.57, 0.85)):
        ax = axs[k]
        for kind, fam in (("ga", "generalized-alpha"), ("bathe", "rho-Bathe")):
            cur = sorted((x["rho"], x["mismatch"]) for x in fe if x["omega_h"] == oh and x["kind"] == kind)
            r, m = np.array(cur).T
            ax.plot(r, 100*m, "-o", color=COL[fam], label=NAME[fam])
        dg = [x for x in fe if x["omega_h"] == oh and x["kind"] == "dg"][0]["mismatch"]
        ax.axhline(100*dg, color=COL["discrete-gradient"], ls="--", label=NAME["discrete-gradient"])
        ax.set_title("FE beam\n" + rf"$\omega_1 h={oh}$")
        ax.set_xlabel(r"$\rho_\infty$"); ax.set_ylim(0, 52)
        tag(ax, "abcd"[k])
    axs[0].set_ylabel("basin mismatch (%)")
    ax = axs[3]
    for kind, fam in (("ga", "generalized-alpha"), ("bathe", "rho-Bathe")):
        cur = sorted((x["rho"], x["mismatch"]) for x in bm if x["alpha"] == 1.0 and x["q20"] == 0.3
                     and x["h"] == 0.4 and x["kind"] == kind)
        r, m = np.array(cur).T
        ax.plot(r, 100*m, "-", color=COL[fam])
    dg = [x for x in bm if x["alpha"] == 1.0 and x["q20"] == 0.3 and x["h"] == 0.4 and x["kind"] == "dg"][0]
    ax.axhline(100*dg["mismatch"], color=COL["discrete-gradient"], ls="--")
    ax.set_title("two-mode beam\n" + r"$\alpha=1$, $q_{20}=0.3$, $h=0.4$")
    ax.set_xlabel(r"$\rho_\infty$"); ax.set_ylim(0, 52)
    tag(ax, "d")
    h = [plt.Line2D([], [], color=COL["generalized-alpha"], label=NAME["generalized-alpha"]),
         plt.Line2D([], [], color=COL["rho-Bathe"], label=NAME["rho-Bathe"]),
         plt.Line2D([], [], color=COL["discrete-gradient"], ls="--", label=NAME["discrete-gradient"])]
    fig.legend(handles=h, frameon=False, loc="lower center", ncol=3, fontsize=6, bbox_to_anchor=(0.5, -0.08))
    fig.tight_layout(w_pad=0.6)
    save(fig, "fig08_structures")


def fig_grid():
    g = J("grid_convergence/summary.json")["results"]
    fig, ax = plt.subplots(figsize=(W1, 2.2))
    for e in g["convergence"]:
        ns = sorted(int(n) for n in e["mismatch_by_n"])
        ax.plot(ns, [100*e["mismatch_by_n"][str(n)] for n in ns], "-o",
                label=rf"$h={e['h']:g}$, $\rho_\infty={e['rho']:g}$")
    ax.set_xscale("log", base=2); ax.set_xticks([100, 200, 400], ["100", "200", "400"])
    ax.set_xlabel("grid points per axis"); ax.set_ylabel("basin mismatch (%)")
    ax.legend(frameon=False, fontsize=5.5, loc="center left", bbox_to_anchor=(1.0, 0.5))
    save(fig, "figA1_grid_convergence")


if __name__ == "__main__":
    fig_basin_maps()
    fig_mismatch_curves()
    fig_optimum_map()
    fig_gain_maps()
    fig_confirmatory()
    fig_controls_equal_step()
    fig_cost_matched()
    fig_structures()
    fig_grid()
