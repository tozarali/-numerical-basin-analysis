# Exploratory structural test: nonlinear von Mises truss

Status: exploratory, not a Phase 4 pass. This is one idealised structure
studied with one and two active coordinates, not two independent structural
benchmarks and not an FE verification.

## Model and experiment

Two identical pin-jointed elastic axial bars connect fixed supports at
(±a,0) to a unit-mass apex (u,y). The stress-free length is L0=sqrt(a²+b²),
and V = k/2 sum_i (Li-L0)² with a=1, b=0.7, k=EA/L0=1. Each active
coordinate has viscous damping 0.15. Wells are (0,±0.7). The two model
configurations are: (i) u fixed at zero (1 DOF), (ii) u free with u0=0.6
and u-dot0=0 (2 DOF). The latter is a large horizontal perturbation,
chosen as an exploratory nonlinear coupling test.

The potential barrier on y=0 is E_b=k(L0-a)²=0.04868887685 for this
geometry. Reference labels are assigned only if physical energy is below
E_b; other states are unresolved, not counted as a well. For the discrete
scheme this endpoint energy test alone does not certify future numerical
capture; T=60 versus T=120 candidate checks were run at rho=0, 0.4 and 1
for every h and grid.

The (y0, y-dot0) domain is [-1.2,1.2]×[-1,1], sampled at cell centres on
48×48 and 96×96 grids. Generalized-alpha was tested at h=0.2,0.4,0.6 and
rho=0:0.1:1, for 132 total runs. D is the median ratio of physical energy
losses for points with 1.5<E0/E_b<8; it is measured at T=8 for h=0.2,0.4
and at the **same time** T=7.8 for h=0.6. D is finite-time, not an
analytically derived per-cycle quantity. Full-grid mismatch and D-subset
mismatch must not be interchanged.

## Reference and solver checks

RK4 h=0.005 versus 0.01 and reference T=60 versus 120 yielded no label
disagreements or exclusions on all four grids/configurations. On 24
deterministically selected ICs per model at n=48, DOP853 and RK4 labels
also agree, and tightening DOP853 tolerances preserves all labels. This
does not certify every IC near the basin boundary. Candidate T=60 versus
120 labels agree at rho=0,0.4,1 in all cases; this is not a time check for
the other eight rho values. No trajectories were unresolved in the sweep.

At rho=0.4 on the 96×96 grid, mismatch for h=0.2,0.1,0.05 is
0.673%, 0.152%, 0.087% (1 DOF) and 0.781%, 0.239%, 0.043% (2 DOF).
Tests verify the force and Hessian by finite differences, stable equilibria,
the barrier along y=0, the symmetric restriction, and short-time
RK4/generalized-alpha convergence against DOP853.

## Observed basin mismatch, full grid, n=96

| Model | h | rho=0 | best tested rho (step 0.1) | best | rho=1 |
|---|---:|---:|---:|---:|---:|
| 1 DOF | 0.4 | 5.23% | 0.3 | 1.41% | 1.89% |
| 1 DOF | 0.6 | 10.76% | 0.3 or 0.4 | 3.21% | 4.25% |
| 2 DOF | 0.4 | 7.14% | 0.3 | 2.47% | 2.54% |
| 2 DOF | 0.6 | 12.48% | 0.4 | 4.51% | 5.45% |

Thus rho affects the outcome in both configurations, but improvement over
rho=1 is small in these examples, especially for 2 DOF at h=0.4. No
continuous-parameter optimum or practical benefit is established by this
coarse rho grid.

## Coupling and D limitations

In the sampled 2-DOF trajectories (20×20 ICs and 41 times from 0 to 20),
the RMS change in vertical restoring force due to the actual u(t), relative
to the RMS force evaluated at u=0 and the same y(t), is 0.230. The median
integral of absolute inter-coordinate work divided by vertical physical
damping loss is 0.832; the median magnitude of *net* work divided by
vertical damping loss is 0.206. The quadrature comparison at 0.1 and 0.05
sampling intervals gives respectively 0.832 and 0.832 for the first
ratio; the maximum energy-balance residual falls from 0.00105 to 0.000261.
These are sample-specific nonlinear coupling diagnostics, not a universal
definition of strong coupling.

Changing from the constrained u0=0 case to the 2-DOF u0=0.6 case changes
23.31% of reference labels on the n=96 grid. Since this also changes initial
energy and excites an extra coordinate, it is not an isolated causal
measurement of coupling.

On the *same* D-selected IC population, for the 2-DOF n=96 case at h=0.6,
the D=1 crossing is rho≈0.339, whereas the lowest mismatch in the tested
grid is at rho=0.9 and 1 (3.70%; rho=0.3 gives 4.44%). The minimum is
shallow and its location changes with the grid: n=48 chooses rho=0.6.
On an 80×80 shifted grid the crossing is 0.336, the minimum is again at
rho=0.9 and 1 (3.58%), and rho=0.35 gives 4.38%. RK4 h=0.005 and 0.01
agree in all 6400 reference labels on that grid. The shifted-grid check was
chosen after seeing the earlier outcome, so is a reproducibility check, not
a prospectively registered validation. Therefore neither a precise optimizer
nor the original universal ±0.1 prediction is established. The full-grid
optimum differs from the subset optimum; this is a population-dependence
result, not a contradiction.

## Gate and next experiment

The effect is present in this idealised structural system and coupling is
measurable, but the original Phase 4 gate requires a separate structural
model and FE case, and the proposed D rule fails to give a stable,
population-independent optimizer here. It remains **open/not passed**.
Next, repeat the 2-DOF experiment on larger or independently sampled
grids, quantify uncertainty in the flat minima, and test a different
geometry/structural model before making a general claim. A measurement/regime
paper remains the fallback if the predictive criterion fails.

Reproducible outputs are in `../data/truss_basins/`: `summary.json`,
`rho_sweep.csv`, reference/candidate `.npz` labels,
`same_population_analysis.json`, `coupling_work.json`, and
`step_refinement.json`, `shifted_grid.json`. Commands are listed in the
project README.
