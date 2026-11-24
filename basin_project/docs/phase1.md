# Phase 1 — validation protocol

Gate: **not yet passed**. The baseline experiment reports observed numbers in
`../data/phase1_results.json`; the script does not automatically pronounce a
gate passed. Its 16×16 independent reference comparison and 24–96 point-per-axis
grids are preliminary, not the planned 200–400 grid study. The grid and
physical-box sizes must be refined further, and additional damping regimes and
families examined. Predeclare quantitative tolerances for the intended
scientific study rather than reverse-engineering them from these results.

The reference checks compare DOP853 at two tolerance settings against
independently implemented RK4 and SciPy RK45 on identical 16×16 initial
conditions; RK4 step and time-horizon comparisons are also made on every
study grid. Resolved comparison counts and exclusions are always reported.
The end-time labels on the generalized-alpha scheme are not proofs of
asymptotic numerical capture.

For c=0.15, h=0.4 and rho=0.4, mismatch at n=24, 48, 96 is 8.33%,
12.07%, 11.00%, respectively; corresponding mean unsigned boundary
distances are 0.0377, 0.0518, 0.0458. Grid convergence is not yet
demonstrated. The h=0.01 versus 0.005 RK4 reference labels and T=80 versus
120 labels agree at all three grids, with no exclusions. On the 16×16
independent sample, DOP853 tolerance, RK45 and RK4 comparisons also show
zero label disagreements, not a proof near every basin boundary. All 16
entropy boxes at this coarse box width contain both labels; other widths
are needed. On the 48×48 grid, candidate mismatch decreases from 12.07%
to 3.91%, 1.30% and 0.26% at h=0.4, 0.2, 0.1 and 0.05, respectively.
Exact values and runtime are stored with the data.

Sb and Sbb here are empirical within-box Shannon entropies with natural logs.
Fixed phase-space box width allows sampling refinement to be separated from
changing the measurement scale. A separate box-width sensitivity study is
needed. With two outcomes Sbb cannot exceed ln(2); failing the ln(2) test
cannot establish smoothness. Symmetric mean distance and Hausdorff are
different unsigned measures; the latter is a maximum, not a mean or a signed
distance. Unresolved cells prevent boundary-distance reporting.

Pilot caveat: `IDEA.md` section 5 originally printed the opposite signs for
the restoring force; it has been corrected to match the pilots and potential.
Pilot fixed-final-time signs were not certified capture, so the preliminary
percentages require re-evaluation. Pilot 4 also compares mismatched times
for h=0.3 and 0.6 because of rounding T/h, but not for h=0.4 or 0.5;
its effect on the D optimum is not yet quantified.

## Matched-time D retest (exploratory)

`basinlab.experiments.energy_loss_retest` reruns c=0.05, 0.15; h=0.2–0.6;
rho step 0.05, with T_D=8 and T_basin=60. A final shortened step now makes
the requested end times exact. D uses the median ΔE_num/ΔE_ref over
0.6<E0<2.0; the initial basin mismatch is reported separately. A focused
rho-by-0.025 sweep compares basin mismatch over that same energy-selected
initial-condition subset, avoiding the earlier mismatch between D and basin
populations.

For c=0.15 and h=0.4, 0.5, 0.6, the T_D=8 D=1 crossings are rho=0.451,
0.473, 0.482; the same-subset basin minima are rho=0.475, 0.425, 0.425.
The discrepancies are 0.024, 0.048, 0.057. At c=0.05 and h=0.5, 0.6,
the crossings are 0.653, 0.659 while the minima are 0.750, 0.775; minimum
mismatch remains 0.429 and 0.344. Thus the original ±0.1 criterion is
promising in this moderate-damping sample, but it does not hold uniformly.
At c=0.05 the basin response is weakly predictive from D near the crossing:
the 200x200 test at c=0.05, h=0.5 has a 0.653 crossing and a 0.750
minimum (mismatch 0.431).

The optimum also depends on which initial-condition population defines
mismatch: at c=0.15, h=0.4 and 0.6, the 100x100 full-grid minima are near
rho=0.35, while the 0.6<E0<2.0 subset minima are 0.475 and 0.425. The
energy subset is the pilot's existing D population, but follow-up should
freeze both populations and the primary mismatch measure before inspecting
new results. Correlations across adjacent rho points share the same initial
conditions and reference; their rank coefficients are descriptive, not
independent-sample significance tests.

The dense 200x200 check for c=0.15, h=0.6 moves the same-subset minimum
from rho=0.425 at n=100 to 0.450 at n=200; D=1 is rho=0.482, a
0.031 discrepancy at n=200. Reference labels at T=60 and T=120 agree on
these dense grids. This is a useful resolution check, not a full grid
convergence study; only two parameter cases and a narrow rho window were
checked, without uncertainty intervals for the location of the shallow
minimum.

D also depends on the measurement horizon. For c=0.15 and h=0.4, 0.5,
0.6, D=1 crosses at rho=0.427, 0.435, 0.419 for T_D=4, versus 0.452,
0.473, 0.482 for T_D=8. At T_D=2 the crossings become 0.563, 0.620,
0.696. Any proposed rule must specify how its energy-loss window is chosen;
the existing finite-horizon median is not yet a horizon-free criterion.

At c=0.15, h=0.4 and T_D=4, a family-matched sweep gives (D-crossing,
same-subset basin minimum): generalized-alpha (0.426, 0.475), HHT
(no crossing in its tested range, best at rho=0.560 with D=0.986),
WBZ (0.673, 0.750), and Newmark (0.957, 0.950). The family-specific
rho refinement is stored separately. These single-case family results
support further testing, not generality: rho samples are correlated through
their shared ICs, the HHT range is bounded, and one pooled rank coefficient
cannot establish an independent collapse. The next decisive experiment
should freeze the time window, IC selection and comparison metrics, then
repeat on a new IC grid with uncertainty on both D crossing and basin-minimum
location. Phase 3 and the structural/FE validation remain open; the proposed
D rule is therefore not a passed gate.

Phase 3's analytical/averaged D and predictive test remain future work; if the
prediction gate fails, retain the measurement/regime-paper alternative.
