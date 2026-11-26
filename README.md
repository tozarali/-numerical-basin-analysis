# Numerical dissipation and basins of attraction in bistable structural dynamics

Code, preregistration protocol and results for the paper

> A. Tozar, *Does the time integrator choose the basin? Numerical dissipation,
> sub-stepping and energy consistency in bistable structural dynamics* (submitted).

The study measures how implicit time integrators change **which equilibrium a
bistable structure reaches** after a transient. For every initial condition the
final equilibrium computed with a candidate integrator is compared with a
converged RK4 reference; the fraction of disagreeing initial conditions is the
*basin mismatch*.

Integrators: generalized-α, HHT, WBZ and dissipative Newmark (parameterized by
the high-frequency spectral radius ρ∞), the ρ∞-Bathe composite method, an
energy-consistent (Gonzalez) discrete-gradient scheme, BDF2 and the trapezoidal
rule. Models: damped double-well oscillator (symmetric and tilted), a two-mode
buckled-beam reduction, and a finite-element buckled clamped beam (corotational
Euler–Bernoulli elements), cross-checked against OpenSeesPy.

## Repository layout

| path | content |
|---|---|
| `basinlab/` | models, integrators, metrics, experiment runner |
| `basinlab/experiments/` | one module per experiment (see `docs/EXPERIMENTS.md`) |
| `protocol/PREREG.md` | preregistration of the confirmatory test (byte-identical to the file used; SHA-256 `bab82e7255431ddb…` when read as UTF-8 text) |
| `data/` | result files: `summary.json` (configuration, environment, results) and per-cell part files |
| `paper/` | `make_tables.py` and `make_figures.py`, which regenerate every number, table and figure of the paper from `data/` |
| `tests/` | unit tests (integrator order/spectral radius, ρ∞-Bathe identities, energy conservation of the discrete gradient, FE model checks) |

## Quick start

```bash
conda env create -f environment.yml      # or: pip install -r requirements.txt
conda activate basin-project
python -m unittest discover -s tests     # ~30 s
python paper/make_tables.py              # numbers and tables from data/
python paper/make_figures.py             # figures from data/
```

Re-running experiments:

```bash
python -m basinlab.run_all --list        # registry with rough run times
python -m basinlab.run_all --quick       # smoke versions (minutes)
python -m basinlab.run_all --only cost_matched_controls --force
python -m basinlab.run_all --confirmatory  # also the locked confirmatory test
```

On Windows, `run_all.cmd` creates the conda environment if needed, runs the
tests and then the runner. The FE cross-check needs `pip install openseespy`
(it reports `skipped` otherwise). Full runs of all experiments take several
hours on a laptop; experiments resume from their part files.

## Exploratory and confirmatory results

Only `confirmatory_family_comparison` is confirmatory: it was specified in
`protocol/PREREG.md` before it was run, and the runner refuses to execute it
unless that file is approved and unmodified in Git. All other experiments are
exploratory. Earlier exploratory experiments that are not reported in the paper
(`phase1`, `energy_loss_retest`, `truss_*`, `dsep_exploratory`, …) are kept for
completeness.

## License

Code: MIT (see `LICENSE`). Result files in `data/` and figures: CC BY 4.0.
If you use this material, please cite the paper (see `CITATION.cff`).
