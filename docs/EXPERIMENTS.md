# Experiments

All experiments write `data/<name>/summary.json` with the configuration, the
software environment, results and run time. Long experiments store one part
file per cell in `data/<name>/parts/` and resume from them. `BASIN_QUICK=1`
runs a small smoke version into `data/<name>_quick/`.

## Reported in the paper

| module | status | content |
|---|---|---|
| `regime_map` | exploratory | double well, 5 damping values × 8 step sizes × 4 families (ρ∞ step 0.05) + BDF2 + discrete gradient, 100×100 initial conditions |
| `grid_convergence` | exploratory | basin mismatch on 100², 200², 400² grids (c = 0.15) |
| `confirmatory_family_comparison` | **confirmatory** | preregistered test of H1 (D = 1 rule) and H2 (discrete gradient) at c = 0.15, h = 0.4 on 32 768 Sobol initial conditions |
| `diagnostic_eval` | exploratory | two-run checks (ρ∞ pair, step pair) without a reference |
| `beam2mode_basins` | exploratory | two-mode beam, three families |
| `fe_arch_validation` | exploratory | FE buckled beam: Euler load, post-buckling amplitude, equilibria, frequencies |
| `fe_crosscheck_openseespy` | exploratory | same beam in OpenSeesPy (statics and frequencies) |
| `fe_arch_basins` | exploratory | FE beam, generalized-α over ρ∞ at three step sizes |
| `asym_well` | exploratory | tilted double well; generalized-α, ρ∞-Bathe, trapezoid, discrete gradient at equal step |
| `beam2mode_controls` | exploratory | two-mode beam with the same controls |
| `fe_dg_control` | exploratory | FE beam with the same controls |
| `cost_matched_controls` | exploratory | equal-cost comparison: ρ∞-Bathe at h against generalized-α and discrete gradient at h/2 |
| `confirmatory_smoke` | smoke | 24×24 smoke version of the confirmatory run (no inferential use) |

## Earlier exploratory work (not reported in the paper)

`phase1`, `phase1_refinement`, `energy_loss_retest`, `focused_rho_retest`,
`d_duration_sensitivity`, `grid_refinement`, `family_collapse_retest`,
`family_optimum_refinement`, `linear_exploratory`, `dsep_exploratory`,
`truss_basins` and its analyses. They informed the design of the
preregistered test (in particular the known weaknesses of the energy-loss
ratio D) and are kept for completeness.

## Definitions

* **Basin mismatch**: fraction of initial conditions with resolved reference and
  candidate labels whose labels differ at the common end time. Unresolved
  states (energy above the barrier at the end time) are excluded and counted,
  never treated as a third basin.
* **Saturated cell**: mismatch ≥ 30 % or candidate basin entropy ≥ 0.9 ln 2.
* **Reliable cell**: no unresolved labels, no label change between T and 2T,
  mismatch < 30 %.
* **Equal cost**: equal number of implicit solves per unit time (ρ∞-Bathe solves
  two systems per step).
