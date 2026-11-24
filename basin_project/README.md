# Basin project

**One command (Windows + conda):** `run_all.cmd` (or `run_all.cmd --quick`
first). Linux/macOS: `./run_all.sh`. Inside an activated environment:
`python -m basinlab.run_all --help`. See `../CALISTIRMA.md` (Turkish) for
options, resume behaviour, the run report and the locked confirmatory test.
Environment: `environment.yml` (conda-forge, Python 3.12).

New exploratory experiments (2026-11-24): `fe_arch_validation`,
`fe_arch_basins` (FE buckled clamped beam; corotational elements with
analytic tangent), `beam2mode_basins`, `regime_map`, `diagnostic_eval`,
`grid_convergence`, and the locked `confirmatory_family_comparison`.

This is the initial Phase 0–1 research infrastructure. Original pilot scripts
remain in `../pilots/` and are not part of this package. The reusable
generalized-alpha, HHT, WBZ and dissipative Newmark parameterizations are
implemented in `basinlab`; the 2-DOF pilot and energy-loss experiments have not
all been migrated. New, independently implemented truss and energy-loss
experiments are available; no results here establish the proposed D criterion.

Run from this directory with Python 3.12 and installed NumPy and SciPy:

```
python -m unittest discover -s tests -v
python -m basinlab.phase1
python -m basinlab.phase1_refinement
python -m basinlab.experiments.energy_loss_retest
python -m basinlab.experiments.focused_rho_retest
python -m basinlab.experiments.d_duration_sensitivity
python -m basinlab.experiments.grid_refinement
python -m basinlab.experiments.family_collapse_retest
python -m basinlab.experiments.family_optimum_refinement
python -m basinlab.experiments.truss_basins
python -m basinlab.experiments.truss_analysis
python -m basinlab.experiments.truss_coupling_work
python -m basinlab.experiments.truss_step_refinement
python -m basinlab.experiments.truss_shifted_grid
python -m basinlab.experiments.linear_exploratory
python -m basinlab.figures.linear_exploratory
python -m basinlab.experiments.dsep_exploratory
```

`data/phase1_results.json` records configuration, environment, convergence,
metrics and runtime; `data/phase1_n*.npz` contain grid coordinates and labels.
The refinement command additionally writes `data/phase1_time_refinement.json`.
The energy-loss retest writes 210 matched-time points to
`data/energy_loss_retest/`; it can take a few minutes on a laptop.
The focused retest refines rho by 0.025 where D=1 and basin minima differed,
and compares both the full grid and the same energy-selected IC subset.
The duration experiment tests the sensitivity of the D crossing to the
finite horizon used to measure energy loss.
The denser-grid experiment reruns the key matched-time cases at 200x200.
The family comparison uses identical initial conditions and samples each
integrator family over the range represented in its configuration.
The family-optimum command refines rho by 0.025 around the candidate minima.
The structural truss sweep takes several minutes; its follow-up analyses
read or complement the outputs in `data/truss_basins/`.
The linear and near-separatrix experiments are exploratory; see
`docs/linear_exploratory.md` and `docs/dsep_exploratory.md`. The plotting
command requires matplotlib. The legacy `dtheory.d_sep` measures a fixed-time
proxy, not necessarily a completed damped cycle.
Rerunning overwrites these named outputs. The experiment is a small validation
case, not the full planned parameter study.

The physical equation is x'' + c x' - x + x^3 = 0, with
E = v^2/2 + (x^2-1)^2/4. For c > 0, a state with E < 1/4 is
trapped in its current physical well. Other states are marked 0 (unresolved),
never counted as a third well. Numerical classifications use this same
endpoint criterion; it does not guarantee that a discrete scheme cannot
subsequently cross the barrier, so the end-time comparison is essential.

Entropy excludes complete boxes containing unresolved samples and reports
their number. Its box width is fixed in phase space as grid density changes;
finite sample errors and box placement remain. The boundary metric uses
midpoints between adjacent different labels, with Euclidean distance in the
unscaled (x,v) plane; it reports unsigned symmetric mean and Hausdorff
distances. It is not a signed displacement and is unavailable for unresolved
grids. These metrics depend on coordinates, window and resolution.

See `docs/literature.md` for the unfinished literature protocol and
`docs/phase1.md` and `docs/truss_results.md` for interpretation and gate status. `figures/`,
`manuscript/`, and `zenodo/` are reserved for later deliverables, not
publication-ready outputs.

## Miniconda preparation launcher

On another Windows computer, copy the complete `basin_project` folder, install
Miniconda and ensure there is internet access for the first environment setup.
Double-click `run_all.cmd`, or in a PowerShell terminal run:

```powershell
powershell -ExecutionPolicy Bypass -File .\run_all.ps1
```

Run the command from the copied `basin_project` directory. The launcher finds
Conda on PATH, via `CONDA_EXE` or at standard Miniconda/Anaconda locations.
It creates `basin-project` from conda-forge if missing; an existing environment
is checked but not modified. No activation is needed. If a custom installation
is not found, open a Miniconda terminal before running it. Creation may need
substantial download time and disk space.

**Preparation only:** it runs the test suite, then an exploratory 24x24 timing
smoke. It does not implement or start the long confirmatory scan. Each smoke
run creates a new directory under `data/confirmatory_preparation/`; existing
results are not overwritten. The reported time is only a linear compute-time
extrapolation from a small grid, omitting bootstrap, I/O and memory effects.
Environment versions and smoke configuration are recorded in `summary.json`.
