"""Generate the LaTeX number macros and tables used in the manuscript from
the result files in data/. Every number quoted in the text is a macro defined
here, so nothing is typed by hand.

Usage (from the repository root):  python paper/make_tables.py
Environment variable BASIN_DATA overrides the data directory."""
import collections
import json
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("BASIN_DATA", HERE.parent / "data"))
OUT = HERE / "generated"
OUT.mkdir(exist_ok=True)


def J(p):
    return json.loads((DATA / p).read_text(encoding="utf-8"))


M = {}


def mac(name, val):
    M[name] = val


def pc(x, d=1):
    return f"{100*x:.{d}f}"


def tabular(spec, header, body):
    return ("\\begin{tabular}{" + spec + "}\\toprule\n" + header + " \\\\\\midrule\n"
            + body + "\n\\bottomrule\\end{tabular}\n")


# ------------------------------------------------------------- regime map
s = J("regime_map/summary.json")["results"]
cells, opt = s["cells"], s["family_optima"]
mac("nCells", len(cells))
mac("nSat", sum(c["saturated"] for c in cells))
mac("pSat", pc(sum(c["saturated"] for c in cells) / len(cells), 0))
dg = {(c["c"], c["h"]): c for c in cells if c["family"] == "discrete-gradient"}
bd = {(c["c"], c["h"]): c for c in cells if c["family"] == "BDF2"}
by = collections.defaultdict(list)
for o in opt:
    by[(o["c"], o["h"])].append(o)
inf, ratios, spread, tune_gain = [], [], [], []
for k, v in sorted(by.items()):
    ok = [o for o in v if not o["best_cell_saturated"]]
    if len(ok) < 4:
        continue
    b = min(o["mismatch_best"] for o in ok)
    if b < 0.01:
        continue
    inf.append(k)
    ratios.append(dg[k]["mismatch"] / b)
    rb = [o["rho_best"] for o in ok]
    spread.append(max(rb) - min(rb))
    tune_gain.append(b / max(o["mismatch_rho1"] for o in ok))
mac("nInf", len(inf))
mac("nDGwins", sum(r < 1 for r in ratios))
mac("medRatio", f"{np.median(ratios):.2f}")
mac("minRatio", f"{min(ratios):.2f}")
mac("maxRatio", f"{max(ratios):.2f}")
mac("medSpread", f"{np.median(spread):.2f}")
mac("maxSpread", f"{max(spread):.2f}")
mac("medTuneGain", f"{np.median(tune_gain):.2f}")
mac("nBDFworse", sum(bd[k]["mismatch"] > min(o["mismatch_best"] for o in by[k]) for k in inf))
exd = {o["family"]: o for o in by[(0.15, 0.4)]}
mac("exDG", pc(dg[(0.15, 0.4)]["mismatch"]))
for f, nm in [("generalized-alpha", "Gen"), ("HHT", "HHT"), ("WBZ", "WBZ"), ("Newmark", "New")]:
    mac("exRho" + nm, f"{exd[f]['rho_best']:.2f}")
    mac("exMis" + nm, pc(exd[f]["mismatch_best"]))
mac("exRhoOne", pc(exd["generalized-alpha"]["mismatch_rho1"]))
mac("rmNc", len({c["c"] for c in cells}))
mac("rmNh", len({c["h"] for c in cells}))

# ------------------------------------------------------- grid convergence
g = J("grid_convergence/summary.json")["results"]
mac("gridMaxChange", f"{max(e['last_change_pct_points'] for e in g['convergence']):.2f}")
mac("gridRefChanges", sum(r["step_check_changes"] + r["time_check_changes"] + r["unresolved"]
                          for r in g["reference_checks"]))
rows = []
for e in g["convergence"]:
    mb = e["mismatch_by_n"]
    rows.append(f"{e['h']:g} & {e['rho']:g} & {100*mb['100']:.2f} & {100*mb['200']:.2f} & "
                f"{100*mb['400']:.2f} & {e['last_change_pct_points']:.2f} \\\\")
(OUT / "tab_grid.tex").write_text(tabular(
    "rrrrrr", r"$h$ & $\rho_\infty$ & $n=100$ & $n=200$ & $n=400$ & change (pts)", "\n".join(rows)))

# -------------------------------------------------------- two-mode beam
b = J("beam2mode_basins/summary.json")["results"]["cells"]
mac("beamCells", len(b))
_bx = [x for x in b if x["alpha"] == 3.0 and x["q20"] == 0.3 and x["h"] == 0.4
       and x["family"] == "generalized-alpha"]
_bm = min(_bx, key=lambda x: x["mismatch"])
mac("beamExRho", f"{_bm['rho']:.1f}")
mac("beamExMis", pc(_bm["mismatch"]))
mac("beamExOne", pc([x for x in _bx if x["rho"] == 1.0][0]["mismatch"]))
_cp = J("beam2mode_basins/summary.json")["results"]["coupling_effect"]
mac("beamCoupOne", pc(_cp[0]["reference_label_change_q20_0_vs_0.3"], 0))
mac("beamCoupThree", pc(_cp[1]["reference_label_change_q20_0_vs_0.3"], 0))

# ---------------------------------------------------------------- FE beam
f = J("fe_arch_basins/summary.json")
mac("feN", f["configuration"]["n"] ** 2)
fc = f["results"]["cells"]
for oh in (0.28, 0.57, 0.85):
    cur = sorted((x["rho"], x["mismatch"]) for x in fc if x["omega_h"] == oh)
    r, mm = min(cur, key=lambda t: t[1])
    t = {0.28: "Lo", 0.57: "Mid", 0.85: "Hi"}[oh]
    mac(f"feBestRho{t}", f"{r:g}")
    mac(f"feBestMis{t}", pc(mm))
    mac(f"feRhoOne{t}", pc(dict(cur)[1.0]))
    mac(f"feRhoZero{t}", pc(dict(cur)[0.0]))
fv = J("fe_arch_validation/summary.json")["results"]
m16 = [m for m in fv["meshes"] if m["n_el"] == 16][0]
mac("feBuckErr", f"{100*abs(m16['buckling_factor_crossing']-1):.1f}")
mac("feAmpErr", f"{100*abs(m16['w_mid_up']/m16['w_mid_first_order']-1):.1f}")
mac("feOmegaOne", f"{m16['well_frequencies'][0]:.2f}")

# OpenSees cross-check
oc = J("fe_crosscheck_openseespy/summary.json")["results"]
if oc.get("status") == "done":
    rd = oc["relative_differences"]
    mac("osAmpRel", f"{rd['midspan_amplitude']:.0e}".replace("e-", r"\times10^{-").replace("e+", r"\times10^{") + "}")
    mac("osOmRel", f"{rd['omega1']:.0e}".replace("e-", r"\times10^{-").replace("e+", r"\times10^{") + "}")
    mac("osAxial", f"{-oc['opensees']['axial_load_factor_final']/oc['euler_load']:.3f}")
else:
    mac("osAmpRel", "n/a"); mac("osOmRel", "n/a"); mac("osAxial", "n/a")

# ------------------------------------------------------------ diagnostic
d = J("diagnostic_eval/summary.json")["results"]
for k, nm in (("two_rho", "Rho"), ("two_h", "H")):
    mac("diag" + nm + "Sp", f"{d[k]['spearman']:.2f}")
    mac("diag" + nm + "Hit", pc(d[k]["hit_rate"], 1))
    mac("diag" + nm + "FA", pc(d[k]["false_alarm_rate"], 0))
    mac("diag" + nm + "N", d[k]["n"])

# ----------------------------------------------------------- confirmatory
cs = J("confirmatory_family_comparison/summary.json")
r = cs["results"]
mac("confSha", cs["configuration"]["prereg_lock"]["sha256"][:16])
mac("confN", cs["configuration"]["n_target"])
h2 = r["H2_discrete_gradient"]
mac("hTwoDG", pc(h2["descriptive"]["discrete_gradient_mismatch"]))
mac("hTwoDiff", f"{100*h2['mismatch_minus_best_reliable_family']:.1f}")
mac("hTwoLo", f"{100*h2['ci95'][0]:.1f}")
mac("hTwoHi", f"{100*h2['ci95'][1]:.1f}")
mac("hTwoBDF", pc(h2["descriptive"]["BDF2_mismatch"]))
mac("hTwoTrap", pc(h2["descriptive"]["trapezoidal_rho_1_mismatch"]))
mac("hTwoGenSix", pc(h2["descriptive"]["family_rho_0.6_mismatch"]["generalized-alpha"]))
fm = {fam: min(v["curve"]) for fam, v in r["families"].items()}
mac("bestAlphaFam", pc(min(fm[x] for x in ("generalized-alpha", "HHT", "WBZ"))))
mac("bestNewmark", pc(fm["Newmark"]))
mac("hOneLo", f"{r['mean_error_difference_vs_baseline_ci95'][0]:.3f}")
mac("hOneHi", f"{r['mean_error_difference_vs_baseline_ci95'][1]:.3f}")
mac("hOneHits", len(r["hits_within_0.1"]))
mac("hOneRel", len(r["reliable_families"]))
names = {"generalized-alpha": r"generalized-$\alpha$", "HHT": "HHT", "WBZ": "WBZ", "Newmark": "Newmark"}
rows = []
for fam, v in r["families"].items():
    cross = "none" if v["D_crossing"] is None else f"{v['D_crossing']:.3f}"
    err = "--" if v["D_crossing"] is None else f"{v['prediction_error']:.3f}"
    rows.append(f"{names[fam]} & {v['optimum']['primary']:.3f} & {cross} & {err} & "
                f"{v['baseline_error']:.3f} & {'yes' if v['reliable'] else 'no'} \\\\")
(OUT / "tab_h1.tex").write_text(tabular(
    "lrrrrr", r"Family & optimum & $D=1$ crossing & $|$crossing$-$opt.$|$ & $|0.6-$opt.$|$ & reliable",
    "\n".join(rows)))
dd = h2["descriptive"]
rows = [("discrete gradient", dd["discrete_gradient_mismatch"]),
        (r"Newmark, best cell ($\rho_\infty=0.95$)", fm["Newmark"]),
        (r"best of generalized-$\alpha$/HHT/WBZ", min(fm[x] for x in ("generalized-alpha", "HHT", "WBZ"))),
        (r"generalized-$\alpha$, $\rho_\infty=0.6$", dd["family_rho_0.6_mismatch"]["generalized-alpha"]),
        (r"trapezoidal rule, $\rho_\infty=1$", dd["trapezoidal_rho_1_mismatch"]),
        ("BDF2 (unreliable)", dd["BDF2_mismatch"])]
(OUT / "tab_h2.tex").write_text(tabular(
    "lr", r"Integrator & mismatch (\%)", "\n".join(f"{a} & {100*v:.1f} \\\\" for a, v in rows)))

# ------------------------------------ extended controls at equal step size
ext = {"SDOF": [], "beam": [], "FE": []}
for s_ in J("asym_well/summary.json")["results"]["summary"]:
    ext["SDOF"].append(s_)
for s_ in J("beam2mode_controls/summary.json")["results"]["summary"]:
    ext["beam"].append(s_)
for s_ in J("fe_dg_control/summary.json")["results"]["summary"]:
    ext["FE"].append(s_)


def spreads(name, keys):
    c_ = J(name + "/summary.json")["results"]["cells"]
    gg = collections.defaultdict(list)
    for x in c_:
        if x["mismatch"] is not None:
            gg[(tuple(x[k] for k in keys), x["kind"])].append(x["mismatch"])
    return ({k: max(v) - min(v) for (k, kd), v in gg.items() if kd == "ga"},
            {k: max(v) - min(v) for (k, kd), v in gg.items() if kd == "bathe"})


sp = {"SDOF": spreads("asym_well", ["eps", "c", "h"]),
      "beam": spreads("beam2mode_controls", ["alpha", "q20", "h"]),
      "FE": spreads("fe_dg_control", ["omega_h"])}
rows = []
labels = {"SDOF": "double well ($\\varepsilon\\in\\{0,0.05,0.1\\}$)", "beam": "two-mode beam", "FE": "FE buckled beam"}
for k, S in ext.items():
    bg = np.array([x["best_ga"] for x in S]); dgv = np.array([x["dg"] for x in S])
    bb = np.array([x["best_bathe"] for x in S]); b0 = np.array([x["bathe_rho0"] for x in S])
    tr = np.array([x["trapezoid"] for x in S])
    sg, sb = sp[k]
    rows.append(f"{labels[k]} & {len(S)} & {int((dgv < bg).sum())} & {np.median(dgv/bg):.2f} & "
                f"{int((bb < bg).sum())} & {np.median(bb/bg):.2f} & {np.median(b0/tr):.2f} & "
                f"{100*np.median(list(sg.values())):.1f} & {100*np.median(list(sb.values())):.1f} \\\\")
    mac({"SDOF": "extSdof", "beam": "extBeam", "FE": "extFe"}[k] + "N", len(S))
allS = ext["SDOF"] + ext["beam"] + ext["FE"]
mac("extN", len(allS))
mac("extDGwins", sum(x["dg"] < x["best_ga"] for x in allS))
mac("extBathewins", sum(x["best_bathe"] < x["best_ga"] for x in allS))
mac("extBzerowins", sum(x["bathe_rho0"] < x["best_ga"] for x in allS))
mac("extBathevsDG", sum(x["best_bathe"] < x["dg"] for x in allS))
mac("extMedBzeroTrap", f"{np.median([x['bathe_rho0']/x['trapezoid'] for x in allS]):.2f}")
allsg = [v for k in sp for v in sp[k][0].values()]; allsb = [v for k in sp for v in sp[k][1].values()]
mac("extSpreadGA", pc(np.median(allsg))); mac("extSpreadBathe", pc(np.median(allsb)))
mac("extSpreadBatheMax", pc(max(allsb)))
(OUT / "tab_ext.tex").write_text(tabular(
    "lrrrrrrrr",
    r"Model & regimes & \multicolumn{2}{c}{DG vs.\ best GA} & \multicolumn{2}{c}{best Bathe vs.\ best GA} & "
    r"Bathe$_0$/TR & \multicolumn{2}{c}{spread over $\rho_\infty$ (pts)} \\" "\n"
    r" & & wins & ratio & wins & ratio & ratio & GA & Bathe",
    "\n".join(rows)))

# tilt effect
aw = {(x["eps"], x["c"], x["h"]): x for x in J("asym_well/summary.json")["results"]["summary"]}
dif = [abs(aw[(e, c, h)]["best_ga"] - aw[(0.0, c, h)]["best_ga"]) for e in (0.05, 0.1)
       for c in (0.05, 0.15, 0.3) for h in (0.2, 0.4, 0.6)]
drho = [abs(aw[(e, c, h)]["best_ga_rho"] - aw[(0.0, c, h)]["best_ga_rho"]) for e in (0.05, 0.1)
        for c in (0.05, 0.15, 0.3) for h in (0.2, 0.4, 0.6)]
mac("tiltMedDiff", pc(np.median(dif))); mac("tiltMaxDiff", pc(max(dif)))
mac("tiltMaxRho", f"{max(drho):.1f}")
fr = [x["fraction_right_well"] for x in J("asym_well/summary.json")["results"]["reference_checks"]
      if x["eps"] == 0.1]
mac("tiltFracLo", pc(min(fr), 0)); mac("tiltFracHi", pc(max(fr), 0))

# ------------------------------------------------- cost-matched comparison
p = DATA / "cost_matched_controls/summary.json"
if p.exists():
    cm = J("cost_matched_controls/summary.json")["results"]
    rows = []
    lab = {"sdof": "double well", "beam": "two-mode beam", "fe": "FE buckled beam"}
    allc = []
    for k in ("sdof", "beam", "fe"):
        if k not in cm:
            continue
        S = cm[k]["summary"]; allc += S
        bb = np.array([x["best_bathe"] for x in S]); bgh = np.array([x["best_ga_half"] for x in S])
        dgh = np.array([x["dg_half"] for x in S]); b0 = np.array([x["bathe_rho0"] for x in S])
        b1 = np.array([x["bathe_rho1"] for x in S]); trh = np.array([x["trapezoid_half"] for x in S])
        rows.append(f"{lab[k]} & {len(S)} & {np.median(bb/bgh):.2f} & {int((bb < bgh).sum())} & "
                    f"{np.median(b0/bgh):.2f} & {np.median(dgh/bb):.2f} & {int((dgh < bb).sum())} & "
                    f"{np.max(np.abs(b1-trh))*100:.2f} \\\\")
    (OUT / "tab_cost.tex").write_text(tabular(
        "lrrrrrrr",
        r"Model & regimes & Bathe/GA$_{h/2}$ & Bathe wins & Bathe$_0$/GA$_{h/2}$ & DG$_{h/2}$/Bathe & DG wins & "
        r"$\max|\mathrm{B}_1-\mathrm{TR}_{h/2}|$ (pts)",
        "\n".join(rows)))
    bb = np.array([x["best_bathe"] for x in allc]); bgh = np.array([x["best_ga_half"] for x in allc])
    dgh = np.array([x["dg_half"] for x in allc])
    mac("cmN", len(allc))
    mac("cmMedBatheGA", f"{np.median(bb/bgh):.2f}")
    mac("cmBatheWins", int((bb < bgh).sum()))
    mac("cmMedDGBathe", f"{np.median(dgh/bb):.2f}")
    mac("cmDGWins", int((dgh < bb).sum()))
    mac("cmMaxBoneTr", f"{100*max(abs(x['bathe_rho1']-x['trapezoid_half']) for x in allc):.2f}")
    sec = [(x["seconds_bathe"], x["seconds_ga_half"], x["seconds_dg_half"]) for x in allc]
    worst = []
    for k in ("sdof", "beam", "fe"):
        if k not in cm:
            continue
        gk = collections.defaultdict(list)
        keys = {"sdof": ("c", "H"), "beam": ("alpha", "q20", "H"), "fe": ("omega_H",)}[k]
        for x in cm[k]["cells"]:
            gk[tuple(x[q] for q in keys)].append(x)
        for v in gk.values():
            wb = max(x["mismatch"] for x in v if x["kind"] == "bathe")
            wg = max(x["mismatch"] for x in v if x["kind"] == "ga")
            worst.append((wb, wg))
    mac("cmMedWorst", f"{np.median([a_/b_ for a_, b_ in worst]):.2f}")
    mac("cmWorstBathe", pc(np.median([a_ for a_, _ in worst])))
    mac("cmWorstGA", pc(np.median([b_ for _, b_ in worst])))
    mac("cmTimeRatioGA", f"{np.median([a/b_ for a, b_, _ in sec]):.2f}")
    mac("cmTimeRatioDG", f"{np.median([c_/a for a, _, c_ in sec]):.2f}")

with open(OUT / "numbers.tex", "w", encoding="utf-8") as fh:
    for k, v in M.items():
        fh.write(f"\\newcommand{{\\{k}}}{{{v}}}\n")
print("wrote", OUT, len(M), "macros")
