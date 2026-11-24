# Research status

Updated: 2026-11-24. A Git repository was initialized at the project root
with user approval; no commit or committed preregistration exists, and all
numerical results remain exploratory. Git for Windows is installed at
`C:\Program Files\Git\cmd\git.exe`; commit identity is not yet configured.

| Gate | Status | Evidence / outstanding work |
|---|---|---|
| 0 Literature | Open | User handles literature; no search performed here. |
| 1 Infrastructure | Open | `basin_project/docs/phase1.md`; complete 200-to-400 convergence, box sensitivity and signed-boundary checks remain. |
| 2 Regime map | Partial | Family and damping sweeps exist; no preregistered regime map or uncertainty for the minima. |
| 3 D mechanism | Open, problematic | `basin_project/docs/linear_exploratory.md`, `basin_project/docs/dsep_exploratory.md`; true near-barrier full cycle absent at c=0.15 for tested phases. |
| 4 Structural / FE | Open | `basin_project/docs/truss_results.md`; 1- and 2-DOF variants are one structure, not separate independent structures; FE validation outstanding. |

## Most important weaknesses

1. D depends on time window and initial-condition population; the 2-DOF
   truss provides a counterexample to a stable, universal optimizer.
2. The actual proposed near-barrier full-cycle D_sep is undefined in the
   tested damped cases; replacing it with a fixed window does not solve the
   window-dependence objection. The phase-error regression is descriptive.
3. Independent confirmation is missing: one c/h family comparison with
   correlated rho points, grid-dependent/shallow optima, no preregistered
   test, separate structure or verified FE case.

The observation that family optima differ is consistent with the earlier
assessment, but neither general rho non-transferability nor D transferability
is established. Equal-looking minimum errors do not establish a common
phase-error floor.

## This session

- Baseline 28/28 and final 35/35 tests passed.
- Generated exploratory linear curves and pooled fit; checked actual
  reference cycle existence and computed two clearly named finite-time
  alternatives. See the two documents above and their data directories.
- `basin_project/protocol/PREREG.md` is an unapproved, uncommitted draft.
  Do not launch confirmatory runs from it.
- Preparation smoke (24x24 grid, 284 selected ICs, three rho values per
  family) and 37/37 tests passed in the Miniconda `basin-project` environment
  (Python 3.12.16). Linear extrapolation to 32768 ICs estimates 45.1 minutes
  compute only, excluding bootstrap, I/O and memory effects; this is not a
  benchmark of the eventual long-run code.
- `basin_project/run_all.cmd` / `run_all.ps1` locate Miniconda and create a
  missing `basin-project` environment; they execute tests and this smoke only.
  They are not yet long-run or confirmatory runners.

## Session 2026-11-24 (afternoon, Claude chat)

- One-shot runner `basinlab.run_all` + `run_all.cmd/.ps1/.sh` + `environment.yml`:
  tests first, then every experiment in its own process, resumable, with
  logs and `data/_runs/<stamp>/REPORT.md`. The launcher no longer stops at
  preparation; the confirmatory test stays locked behind an approved,
  committed PREREG.
- New exploratory experiments: WP1 `grid_convergence`, WP2 `regime_map`
  (incl. BDF2 and discrete-gradient controls), WP4 `beam2mode_basins`,
  WP5 `fe_arch_validation` + `fe_arch_basins`, WP6 `diagnostic_eval`.
  All passed quick (smoke) runs only; no full-size result exists yet.
- FE model changed: a stress-free shallow clamped arch (rise 0.1, EA/EI=1e4)
  was found NOT bistable when unloaded, so WP5 now uses a buckled
  clamped-clamped beam (prescribed end shortening 1.6 Delta_cr). Verified
  against closed forms: Euler load (16 el: +1.7%, 32 el: <1%) and
  first-order post-buckling amplitude (~1% low, as expected).
- Fixes: discrete-gradient Newton failed on the full grid (now exact
  difference quotient); signed boundary shift sign logic was wrong (rewritten,
  tested); per-IC masked Newton for MDOF gen-alpha (rho=1 needed >30 iterations
  in FE). 52 tests pass (4 Windows-only skipped here).
- Still open: rho-infinity Bathe (needs original paper), asymmetric well
  (epsilon), FE literature benchmark/openseespy cross-check, PREREG approval.

## Next

NPZ outputs remain local and are ignored by Git. Finalize the primary
comparison and event-based metric before committing a protocol. The user
intends to run long experiments locally: prepare a 24x24 smoke test, time
estimate, memory profiling, resumable long-run code, and reproducible command
before any long run, and obtain approval as required.
