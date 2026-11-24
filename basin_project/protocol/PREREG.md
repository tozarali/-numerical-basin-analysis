# Draft preregistration: independent family comparison

**DRAFT — NOT APPROVED, NOT COMMITTED, NOT AUTHORIZATION TO RUN.**
Prepared 2026-11-24 after seeing exploratory results. This proposal is for
a prospective new population, not retroactive confirmation of existing data.
Critical design choices below require resolution before final approval and
commit; subsequent threshold changes require a new version.

## Broader `opus.md` programme

The focused design below is a proposed first confirmatory comparison, not
a substitute for the full WP2/WP3 programme. The broader design must cover
c={0.05, 0.1, 0.15, 0.25, 0.4}, eight log-spaced h values from 0.05
to 1, epsilon={0, 0.1}, and at least 21 admissible rho points per family.
It must also include BDF2, rho-infinity Bathe, trapezoidal and a discrete
gradient control, together with mismatch, entropy, signed boundary shift,
reference checks and saturation flags. The asymmetric force/energy and its
barrier, Bathe, some metrics and the full runner are not yet validated;
they require explicit definitions and tests before preregistration can be
finalized. WP3's multi-regime gate and baseline must be specified on those
reliable cells; the limited gate below cannot pass it.

The user intends to run long computations locally. Before launching them,
prepare a reproducible script, 24x24 smoke test, elapsed-time and memory
estimate, checkpoints, and a user-approved committed final protocol.

## Question and limits

At c=0.15 and h=0.4, on an independent Sobol IC population, do finite-time
D=1 crossings predict observed subset basin optima across four integrator
families better than a fixed common rho=0.6? This is one regime, not proof
of general transferability, a causal phase mechanism, or a full regime map.
The proposed real-cycle D_sep cannot be a primary measure: it is undefined
in the tested c=0.15 near-barrier reference trajectories.

## Proposed design

- Model: x''+0.15 x'-x+x³=0; Eb=0.25.
- Population: scrambled Sobol on [-2,2]², seed 20261124, generate 2^17
  points; take the first 32768 with 0.6<E0<2. If insufficient, stop rather
  than silently change the procedure. Record coordinates and selection.
- Families: generalized-alpha rho=0..1 step .025; HHT rho=.51 plus
  .525..1 step .025; WBZ rho=.05..1 step .025; dissipative Newmark
  rho=.6..1 step .025. Record nominal and measured mapping. The latter
  three ranges are restricted and rho=0.6 lies in each.
- Reference labels: RK4 h=.005, T=60. Separately compare reference labels
  with h=.0025 and T=120 on the same ICs; compare candidate T=60/120 for
  all rho values. Count unresolved labels; do not treat them as a third well.
- D: median of (E0-E_num(4))/(E0-E_ref(4)) on this same IC population,
  exact common physical end time. Record numerator/denominator per IC;
  exclude only invalid or reference-loss <=1e-12, and report counts.
- Primary outcome: fraction of resolved candidate/reference labels that
  disagree at T=60. A cell is reliable only if reference refinements and
  candidate T checks show zero label changes, no unresolved labels remain,
  and mismatch is <30%. A failed check is reported, not silently filtered.
- Estimate the observed optimum using discrete argmin and a local quadratic
  bootstrap fit; report both, fit failures and boundary optima. Select a
  unique D=1 crossing by linear interpolation on the prespecified rho grid;
  missing or multiple crossings are prediction failures, never extrapolate.
- Joint paired bootstrap resamples IC indices identically across all
  families and rho values, recalculating mismatch, D medians, crossings and
  optima per replicate. Proposal: 2000 replicates, seed 20261125.

## Proposed gate

Among reliable family cases, at least 80% must have a unique D crossing within
0.1 of the observed optimum; missing/multiple crossings count as failures.
Additionally, the paired bootstrap 95% interval for mean absolute prediction
error (D crossing minus fixed rho=0.6 baseline) must lie strictly below zero.
If fewer than three families are reliable, the gate is inconclusive, not
passed. This four-family, single-regime gate is a limited prospective test;
it cannot by itself pass the broader WP3 multi-regime gate in `opus.md`.

## Unresolved before approval / commit

1. Specify the exact quadratic fit window, duplicate minima policy, and
   handling of bootstrap replicates with missing/multiple D crossings.
2. Decide whether zero changes in all label reference checks is feasible;
   if not, specify a quantitative reliability rule *before* confirmation.
3. Validate that all needed per-IC D and label arrays, bootstrap and
   reference checks are implemented and tested; estimate memory/disk.
4. Run and record a 24x24 smoke test and runtime extrapolation before any
   >10-minute run; obtain the user's approval for the long run.
5. Decide repository and large-data strategy, obtain user approval, finalize
   this protocol and commit it before starting confirmatory computation.

Half-conservative-period and first-turn proxies, phase regression, full-grid
analyses and additional h/c regimes are exploratory until separately
preregistered. A failed gate does not justify changing its thresholds.

## Proposed resolutions implemented in code (2026-11-24, still DRAFT)

`basinlab/experiments/confirmatory_family_comparison.py` implements the
design above. Where the draft was open it uses the defaults below; copy,
change or reject them here before approval. The runner refuses the full
run until this file has a `Status: APPROVED` line, no draft marker, and is
committed unmodified in Git.

1. Quadratic fit window: 3 grid points each side of the discrete minimum
   (+-0.075 on the 0.025 grid); duplicate minima: smallest rho.
2. Primary observed optimum: quadratic vertex; discrete argmin if the fit
   is flagged (edge, not convex, vertex outside window). Both reported.
3. Family reliability: reference checks pass and every cell within +-0.1
   of the discrete optimum is reliable (zero unresolved, zero T=60/120
   candidate changes, mismatch < 30%). A strict all-cells flag is reported too.
4. Bootstrap replicates with a missing or multiple D crossing: that
   family's prediction error is set to the worst error possible on its rho
   grid (conservative, counts against the D rule).
5. Gate statistic: mean over reliable families of (|D crossing - optimum|
   - |0.6 - optimum|); 2000 paired replicates, seed 20261125.
Smoke check (512 ICs, coarse grid) ran end to end; it is not evidence.
