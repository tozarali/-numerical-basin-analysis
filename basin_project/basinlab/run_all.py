"""Run every experiment of the project with one command.

    python -m basinlab.run_all                 # all exploratory experiments (resumable)
    python -m basinlab.run_all --quick         # smoke versions only (minutes)
    python -m basinlab.run_all --list          # registry with rough time estimates
    python -m basinlab.run_all --only regime_map,diagnostic_eval
    python -m basinlab.run_all --stage new     # legacy | new | all (default all)
    python -m basinlab.run_all --force         # rerun even if outputs exist
    python -m basinlab.run_all --confirmatory  # ALSO the locked confirmatory test
                                               # (refuses unless PREREG approved+committed)

Each experiment runs in its own process (`python -m <module>`), so one failure
or memory spike does not stop the others. Console output is mirrored to
data/_runs/<timestamp>/<name>.log; a REPORT.md and report.json summarise the
batch. An experiment is skipped when its outputs exist and its source file
fingerprint has not changed since the last successful run (state in
data/_runs/state.json); legacy outputs that predate the runner are kept
unless --force. Long experiments also resume internally from part files.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "data" / "_runs"


@dataclass
class Exp:
    name: str
    module: str
    stage: str                    # "legacy", "new", "confirmatory"
    outputs: list                 # relative to ROOT; quick runs use *_quick dirs
    minutes: float                # rough full-run estimate on a laptop
    deps: list = field(default_factory=list)
    quick_ok: bool = False        # supports BASIN_QUICK=1
    quick_outputs: list = field(default_factory=list)
    always_quick: bool = False    # smoke-only entry: run with BASIN_QUICK=1 even in full runs


E = "basinlab.experiments."
REGISTRY = [
    # --- legacy exploratory experiments (reproduce existing results) ------
    Exp("phase1", "basinlab.phase1", "legacy", ["data/phase1_results.json"], 2),
    Exp("phase1_refinement", "basinlab.phase1_refinement", "legacy",
        ["data/phase1_time_refinement.json"], 1, ["phase1"]),
    Exp("energy_loss_retest", E + "energy_loss_retest", "legacy",
        ["data/energy_loss_retest/summary.json"], 4),
    Exp("focused_rho_retest", E + "focused_rho_retest", "legacy",
        ["data/focused_rho_retest/summary.json"], 2),
    Exp("d_duration_sensitivity", E + "d_duration_sensitivity", "legacy",
        ["data/d_duration_sensitivity/summary.json"], 1),
    Exp("grid_refinement", E + "grid_refinement", "legacy", ["data/grid_refinement/summary.json"], 8),
    Exp("family_collapse_retest", E + "family_collapse_retest", "legacy",
        ["data/family_collapse_retest/summary.json"], 1),
    Exp("family_optimum_refinement", E + "family_optimum_refinement", "legacy",
        ["data/family_optimum_refinement/summary.json"], 1),
    Exp("linear_exploratory", E + "linear_exploratory", "legacy",
        ["data/linear_exploratory/summary.json"], 0.2, ["family_optimum_refinement"]),
    Exp("linear_exploratory_figure", "basinlab.figures.linear_exploratory", "legacy",
        ["data/linear_exploratory/curves.pdf"], 0.2, ["linear_exploratory"]),
    Exp("dsep_exploratory", E + "dsep_exploratory", "legacy", ["data/dsep_exploratory/summary.json"], 2),
    Exp("truss_basins", E + "truss_basins", "legacy", ["data/truss_basins/summary.json"], 5),
    Exp("truss_analysis", E + "truss_analysis", "legacy",
        ["data/truss_basins/same_population_analysis.json"], 1, ["truss_basins"]),
    Exp("truss_coupling_work", E + "truss_coupling_work", "legacy",
        ["data/truss_basins/coupling_work.json"], 1, ["truss_basins"]),
    Exp("truss_step_refinement", E + "truss_step_refinement", "legacy",
        ["data/truss_basins/step_refinement.json"], 2, ["truss_basins"]),
    Exp("truss_shifted_grid", E + "truss_shifted_grid", "legacy",
        ["data/truss_basins/shifted_grid.json"], 3, ["truss_basins"]),
    # --- new exploratory experiments (WP1, WP2, WP4, WP5, WP6) --------------
    Exp("fe_arch_validation", E + "fe_arch_validation", "new",
        ["data/fe_arch_validation/summary.json"], 2, quick_ok=True,
        quick_outputs=["data/fe_arch_validation_quick/summary.json"]),
    Exp("beam2mode_basins", E + "beam2mode_basins", "new",
        ["data/beam2mode_basins/summary.json"], 15, quick_ok=True,
        quick_outputs=["data/beam2mode_basins_quick/summary.json"]),
    Exp("regime_map", E + "regime_map", "new", ["data/regime_map/summary.json"], 70,
        quick_ok=True, quick_outputs=["data/regime_map_quick/summary.json"]),
    Exp("diagnostic_eval", E + "diagnostic_eval", "new", ["data/diagnostic_eval/summary.json"], 0.5,
        ["regime_map"], quick_ok=True, quick_outputs=["data/diagnostic_eval_quick/summary.json"]),
    Exp("fe_arch_basins", E + "fe_arch_basins", "new", ["data/fe_arch_basins/summary.json"], 25,
        ["fe_arch_validation"], quick_ok=True,
        quick_outputs=["data/fe_arch_basins_quick/summary.json"]),
    Exp("grid_convergence", E + "grid_convergence", "new", ["data/grid_convergence/summary.json"], 70,
        quick_ok=True, quick_outputs=["data/grid_convergence_quick/summary.json"]),
    Exp("confirmatory_smoke", E + "confirmatory_family_comparison", "new",
        ["data/confirmatory_family_comparison_smoke/summary.json"], 0.5, quick_ok=True,
        quick_outputs=["data/confirmatory_family_comparison_smoke/summary.json"], always_quick=True),
    # --- locked confirmatory test --------------------------------------------
    Exp("confirmatory_family_comparison", E + "confirmatory_family_comparison", "confirmatory",
        ["data/confirmatory_family_comparison/summary.json"], 40),
]
BY_NAME = {e.name: e for e in REGISTRY}


def source_file(exp):
    return ROOT / (exp.module.replace(".", "/") + ".py")


def fingerprint(exp, quick):
    h = hashlib.sha256(source_file(exp).read_bytes())
    h.update(b"quick" if quick else b"full")
    return h.hexdigest()[:16]


def outputs_exist(exp, quick):
    outs = exp.quick_outputs if quick else exp.outputs
    return bool(outs) and all((ROOT / o).exists() for o in outs)


def summary_completed(exp, quick):
    """True unless a summary.json says the run failed."""
    for o in (exp.quick_outputs if quick else exp.outputs):
        p = ROOT / o
        if p.name == "summary.json" and p.exists():
            try:
                if str(json.loads(p.read_text(encoding="utf-8")).get("outcome", "completed")).startswith("failed"):
                    return False
            except (ValueError, OSError):
                return False
    return True


def select(args):
    stages = {"legacy": ["legacy"], "new": ["new"], "all": ["legacy", "new"]}[args.stage]
    if args.confirmatory:
        stages = stages + ["confirmatory"]
    chosen = [e for e in REGISTRY if e.stage in stages]
    if args.quick:
        chosen = [e for e in chosen if e.quick_ok]
    if args.only:
        names = [n.strip() for n in args.only.split(",") if n.strip()]
        unknown = [n for n in names if n not in BY_NAME]
        if unknown:
            raise SystemExit(f"unknown experiment(s): {unknown}; see --list")
        chosen = [BY_NAME[n] for n in names]
    if args.skip:
        skip = {n.strip() for n in args.skip.split(",")}
        chosen = [e for e in chosen if e.name not in skip]
    return chosen


def load_state():
    p = RUNS / "state.json"
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def save_state(state):
    RUNS.mkdir(parents=True, exist_ok=True)
    tmp = RUNS / "state.json.tmp"
    tmp.write_text(json.dumps(state, indent=1), encoding="utf-8")
    tmp.replace(RUNS / "state.json")


def peak_memory_watch(pid, box, stop):
    try:
        import psutil
    except ImportError:
        box["peak_mb"] = None
        return
    try:
        proc = psutil.Process(pid)
        peak = 0
        while not stop.is_set():
            try:
                peak = max(peak, proc.memory_info().rss)
            except psutil.Error:
                break
            time.sleep(0.5)
        box["peak_mb"] = round(peak/2**20, 1)
    except psutil.Error:
        box["peak_mb"] = None


def run_one(exp, quick, log_dir):
    env = os.environ.copy()
    env["BASIN_QUICK"] = "1" if (quick or exp.always_quick) else "0"
    env["PYTHONUNBUFFERED"] = "1"
    log = log_dir / f"{exp.name}.log"
    started = time.perf_counter()
    box, stop = {}, threading.Event()
    with log.open("w", encoding="utf-8") as fh:
        fh.write(f"$ {sys.executable} -m {exp.module}  (BASIN_QUICK={env['BASIN_QUICK']})\n")
        proc = subprocess.Popen([sys.executable, "-m", exp.module], cwd=ROOT, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, encoding="utf-8", errors="replace")
        watcher = threading.Thread(target=peak_memory_watch, args=(proc.pid, box, stop), daemon=True)
        watcher.start()
        tail = []
        for line in proc.stdout:
            print(f"  [{exp.name}] {line}", end="", flush=True)
            fh.write(line)
            tail = (tail + [line.rstrip()])[-15:]
        code = proc.wait()
        stop.set()
        watcher.join(timeout=2)
    return {"returncode": code, "seconds": round(time.perf_counter() - started, 1),
            "peak_mb": box.get("peak_mb"), "log": str(log.relative_to(ROOT)), "tail": tail}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true", help="smoke sizes, writes *_quick outputs")
    ap.add_argument("--stage", choices=["legacy", "new", "all"], default="all")
    ap.add_argument("--only", help="comma-separated experiment names")
    ap.add_argument("--skip", help="comma-separated experiment names to leave out")
    ap.add_argument("--force", action="store_true", help="rerun even when outputs are up to date")
    ap.add_argument("--confirmatory", action="store_true",
                    help="also run the locked confirmatory test (needs approved, committed PREREG)")
    ap.add_argument("--stop-on-error", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    if args.list:
        print(f"{'name':32s} {'stage':13s} {'~min':>6s}  quick  deps")
        for e in REGISTRY:
            print(f"{e.name:32s} {e.stage:13s} {e.minutes:6.1f}  {'yes' if e.quick_ok else ' - ':5s}  "
                  f"{','.join(e.deps)}")
        print("Minutes are rough laptop estimates for full runs, not guarantees.")
        return 0

    chosen = select(args)
    state = load_state()
    plan = []
    for e in chosen:
        fp = fingerprint(e, args.quick)
        key = e.name + ("@quick" if args.quick else "")
        prev = state.get(key, {})
        if args.force:
            reason = None
        elif prev.get("ok") and prev.get("fingerprint") == fp and outputs_exist(e, args.quick):
            reason = "up to date"
        elif not prev and outputs_exist(e, args.quick) and summary_completed(e, args.quick):
            reason = "existing outputs (predate runner; --force to rerun)"
        else:
            reason = None
        plan.append((e, key, fp, reason))
    todo = [p for p in plan if p[3] is None]
    est = sum(e.minutes for e, *_ in todo) if not args.quick else None
    print(f"Plan: {len(todo)} to run, {len(plan) - len(todo)} skipped"
          + (f"; rough estimate {est:.0f} min" if est else ""), flush=True)
    for e, key, fp, reason in plan:
        print(f"  {'SKIP' if reason else 'RUN '} {e.name}" + (f"  ({reason})" if reason else ""))
    if args.dry_run:
        return 0

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + ("_quick" if args.quick else "")
    log_dir = RUNS / stamp
    log_dir.mkdir(parents=True, exist_ok=True)
    results, failed, ran = [], set(), set()
    for e, key, fp, reason in plan:
        row = {"name": e.name, "module": e.module, "stage": e.stage}
        if reason and any(d in ran for d in e.deps):
            reason = None                      # a dependency was just recomputed
        if reason:
            row.update(status="skipped", reason=reason)
        elif any(d in failed for d in e.deps):
            row.update(status="blocked", reason=f"dependency failed: {[d for d in e.deps if d in failed]}")
            failed.add(e.name)
        else:
            print(f"\n=== {e.name} ({e.module}) ===", flush=True)
            res = run_one(e, args.quick, log_dir)
            ok = res["returncode"] == 0
            ran.add(e.name)
            row.update(status="ok" if ok else "FAILED", **res)
            state[key] = {"ok": ok, "fingerprint": fp, "finished": datetime.now(timezone.utc).isoformat(),
                          "seconds": res["seconds"], "log": res["log"]}
            save_state(state)
            if not ok:
                failed.add(e.name)
                print(f"!!! {e.name} failed (exit {res['returncode']}); see {res['log']}", flush=True)
                if args.stop_on_error:
                    results.append(row)
                    break
        results.append(row)

    report = {"started_utc": stamp, "quick": args.quick, "python": sys.version.split()[0],
              "conda_prefix": os.environ.get("CONDA_PREFIX"), "results": results}
    (log_dir / "report.json").write_text(json.dumps(report, indent=1), encoding="utf-8")
    lines = [f"# Run report {stamp}", "", f"Python {report['python']}, conda: {report['conda_prefix']}",
             "", "| experiment | status | seconds | peak MB | note |", "|---|---|---:|---:|---|"]
    for r in results:
        note = r.get("reason") or ("" if r["status"] == "ok" else (r.get("tail") or [""])[-1][:120])
        lines.append(f"| {r['name']} | {r['status']} | {r.get('seconds', '')} | "
                     f"{r.get('peak_mb', '') or ''} | {note} |")
    (log_dir / "REPORT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (RUNS / "LATEST.txt").write_text(str(log_dir.relative_to(ROOT)) + "\n", encoding="utf-8")
    n_fail = sum(r["status"] in ("FAILED", "blocked") for r in results)
    print(f"\nDone. {sum(r['status'] == 'ok' for r in results)} ok, "
          f"{sum(r['status'] == 'skipped' for r in results)} skipped, {n_fail} failed/blocked.")
    print(f"Report: {(log_dir / 'REPORT.md').relative_to(ROOT)}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
