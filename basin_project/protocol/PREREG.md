# Preregistration: independent confirmatory test at one regime (c=0.15, h=0.4)

Status: APPROVED
Version: 1.0 (prepared 2026-11-24)

This file becomes binding when the line above reads `Status: APPROVED` and
the file is committed unmodified in Git (`python protocol/approve_prereg.py`
does both after an explicit confirmation by the author). Any later change to
a threshold, population, grid or rule requires a new version and must be
reported as a deviation.

## 0. What is already known (exploratory, not confirmatory)

Exploratory runs (regime map, grid convergence, two-mode beam, FE beam,
diagnostic evaluation; all on populations other than the one below) suggested:

- the best rho-infinity differs between integrator families and across
  regimes;
- an energy-consistent discrete-gradient integrator gave lower basin
  mismatch than the best-tuned family member in most informative SDOF cells;
- the D-based selection rule has known weaknesses (window and population
  dependence, light-damping breakdown).

This test is prospective on a new Sobol population. It confirms or rejects
two hypotheses at one regime. It does not establish general transferability,
a causal mechanism, or behaviour in other models; those remain exploratory.

## 1. Model, population, references

- Model: x'' + 0.15 x' - x + x^3 = 0 (barrier Eb = 0.25), h = 0.4, T = 60.
- Population: scrambled Sobol on [-2,2]^2, seed 20261124, generate 2^17 points,
  keep the first 32768 with 0.6 < E0 < 2.0. If fewer exist, stop; do not
  change the procedure silently.
- Reference labels: RK4 with h = 0.005 to T = 60. Checks on the same ICs:
  h = 0.0025, and T = 120. The reference is accepted only if both checks
  change zero labels and no label is unresolved. Unresolved is never a
  third basin.
- Common physical end time for all integrators (exact, not rounded).

## 2. Integrators

- Families (nominal rho grids, step 0.025): generalized-alpha 0..1;
  HHT 0.51 and 0.525..1; WBZ 0.05..1; dissipative Newmark 0.6..1.
  rho = 0.6 lies in each grid. Nominal and measured spectral radius are
  recorded.
- Controls: energy-consistent discrete gradient (Gonzalez-type, exact
  difference quotient), BDF2, trapezoidal rule (generalized-alpha rho = 1).
- Not included in this test: rho-infinity-Bathe, asymmetric well, FE/beam
  models (planned as separate exploratory or later registered studies).

## 3. Outcomes and reliability

- Primary outcome: fraction of resolved candidate/reference labels that
  disagree at T = 60 (basin mismatch).
- A cell (family, rho) or control is reliable only if: zero unresolved
  labels, zero label changes between candidate T = 60 and T = 120, and
  mismatch < 30 %. Failed checks are reported, never filtered silently.
- Observed optimum of a family: local quadratic vertex over 3 grid points on
  each side of the discrete minimum; if the fit is flagged (edge, not
  convex, vertex outside window), the discrete argmin is used. Duplicate
  minima: smallest rho. Both values are reported.
- A family is reliable if the reference passes and every cell within
  +-0.1 rho of its discrete optimum is reliable.
- D: median over valid ICs of (E0 - E_num(4)) / (E0 - E_ref(4)); exclude
  only invalid or reference loss <= 1e-12 and report counts. D crossing:
  unique linear-interpolated crossing of D = 1 on the grid; missing or
  multiple crossings are prediction failures; never extrapolate.
- Bootstrap: 2000 paired resamples of IC indices (identical across all
  integrators and rho), seed 20261125, recomputing every statistic per
  replicate.

## 4. Hypothesis H1 (the D rule; limited gate)

H1: D = 1 crossings predict the observed family optima better than the fixed
common rho = 0.6.

Pass only if (a) among reliable families at least 80 % have a unique D
crossing within 0.1 of the observed optimum, and (b) the paired bootstrap 95 %
interval of the mean over reliable families of
(|D crossing - optimum| - |0.6 - optimum|) lies strictly below zero.
In bootstrap replicates a missing or multiple crossing counts as the largest
possible error on that family's grid. If fewer than three families are
reliable the result is inconclusive, not passed.

Expectation recorded in advance: H1 may well fail. A failed gate is reported
as a negative result and thresholds are not changed.

## 5. Hypothesis H2 (energy-consistent integration)

H2: at this regime the discrete-gradient integrator has lower basin mismatch
than the best-tuned member of the four families.

- Statistic: mismatch(discrete gradient) minus the minimum, over all reliable
  family cells, of the family mismatch curve. The minimum is recomputed in every
  bootstrap replicate. This favours the tuned families (it picks the best of
  about 150 noisy cells with knowledge of the reference), so the test is
  conservative against H2.
- Pass only if the paired bootstrap 95 % interval of the statistic lies
  strictly below zero and the discrete-gradient control is reliable (Section 3).
  Otherwise failed or inconclusive (unreliable control).
- Descriptive, not gated: mismatch of discrete gradient, BDF2, trapezoidal,
  and each family at rho = 0.6.

## 6. Reporting rules

- Report both hypotheses whatever the outcome, with the full per-cell tables,
  reliability flags, reference checks, seeds, file hash and Git commit.
- Everything not named in Sections 4 and 5 is exploratory, including other
  regimes (c, h), other models, phase-error analyses, and half-period or
  first-turn proxies.
- The limited gates here cannot by themselves support multi-regime
  claims. A multi-regime claim needs its own preregistration.
- The runner stores this file's SHA-256 and the Git commit in the output.
  Runtime estimate and a 24x24 smoke test (`confirmatory_smoke`) precede
  the full run.
