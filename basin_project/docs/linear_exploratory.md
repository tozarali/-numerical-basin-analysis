# Linear analysis and two-component fit

Status: exploratory, 2026-11-24. Not a preregistered or independent test.

The existing 3-state linear amplification map uses the same alpha-method
convention as `integrators.step`. The linearized well has omega=sqrt(2),
xi=0.15/(2 sqrt(2)); for h=0.4, Omega h=0.5657. Total modal damping
includes physical damping. Relative period error is relative to the natural
period, not a phase lag accumulated along nonlinear trajectories.
`data/linear_exploratory/curves.csv` and `curves.pdf` cover four families,
rho targets 0.6 and 1, undamped and physically damped, Omega h=0.02..3.
The map uses the unique positive-imaginary eigenvalue when present; complex
root identity is not independently tracked through high-frequency
coalescences. Critical/overdamped continuous systems and zero frequency
are explicitly marked nonoscillatory, even if a numerical parasitic pair
appears. A full closed-form generalized-alpha comparison remains outstanding.

The code tests the matrix against an independently rearranged direct update
for all families, the trapezoidal complex eigenvalue against the bilinear
closed form, and nonoscillatory cases. `rho_at_1000` is measured separately:
the existing parameter-mapping return for generalized-alpha is nominal rho,
not the actual finite-frequency spectral radius.

## Descriptive regression

The existing `family_optimum_refinement/rho_refinement.csv` supplies 45
rows (11 generalized-alpha, 8 HHT, 17 WBZ, 9 Newmark). They share the same
100x100 IC grid; the selected band 0.6<E0<2 has 5040 grid points. The
source stores aggregate mismatch and D values, not IC-level errors or the
per-row valid-D count. The sampled rho ranges were selected for refinement
around previously observed optima, and observations along rho are dependent.
The single c=0.15, h=0.4 regime and one IC population limit transferability.

OLS, pooled in-sample, predicts subset mismatch:

| Predictors | R² |
|---|---:|
| abs(D_median-1) | 0.9449 |
| abs(relative linear period error) | 0.0979 |
| both | 0.9669 |
| both + interaction | 0.9847 |

Partial Spearman of mismatch and abs(D-1), controlling period error, is
0.815; mismatch and period error, controlling abs(D-1), is 0.598. These
are descriptive residual-rank correlations, not significance tests.
The interaction fit is poorly conditioned (raw design condition number
about 11180; scaling also affects this number). No independent validation,
confidence intervals or causal decomposition are available. The improved
in-sample fit does **not** show that phase error causes the residual basin
error, or that D transfers across regimes or structures.

Reproduce via `python -m basinlab.experiments.linear_exploratory` and
`python -m basinlab.figures.linear_exploratory`; outputs are overwritten.
