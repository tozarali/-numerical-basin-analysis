# Near-separatrix loss: exploratory event check

Status: exploratory, 2026-11-24. The proposed full-cycle D_sep is undefined
in the tested c=0.15 cases, not a confirmed prediction criterion.

## Definition and reference

For each delta in {0.02, 0.05, 0.1}, start at x=0, v>0 with
E0=0.25(1+delta), on the same initial phase. A full cycle means the next
positive-direction x=0 crossing after the initial crossing. DOP853 records
the first downward energy crossing E=0.25, the first v=0 turning point and
all positive-direction x=0 crossings through T=60. Its energy balance is
checked against the integrated c*v² dissipation, and reference tolerances
are tightened from 1e-11/1e-13 to 1e-12/1e-14.

| delta | barrier capture time | first turn time | full cycle |
|---:|---:|---:|---|
| 0.02 | 1.6827 | 4.3060 | none |
| 0.05 | 1.6945 | 3.8025 | none |
| 0.10 | 1.7158 | 3.4162 | none |

After physical energy falls below the barrier, monotonic energy loss for
c>0 prevents any later x=0 crossing. The reference returns no full cycle
in all three cases. Final energy-balance residuals are at most 1.1e-12;
the largest change in capture time under tightened tolerances is 2.4e-11.
These checks do not establish numerical event accuracy for other models,
initial phases or damping values. A test at c=0.001 also confirms that the
event code can detect a cycle when one actually occurs.

## Two explicitly finite-time proxies

For h=0.4, 0.5, 0.6 and four families, 21 rho points each, compare energy
loss at (a) the **reference first-turn time** and (b) half the undamped
orbit period. Both use the same end time for reference and numerical
solution and shorten the final numerical step. Neither is a full-cycle
measure, a numerical-event-matched measure, or window independent. The
historical `dtheory.d_sep` function is proxy (b), despite its name.
All 756 cases and every interpolated crossing are saved in
`data/dsep_exploratory/ratios.csv` and `summary.json`.

The first-turn proxy has no D=1 crossing in any of the 36 family/h/delta
cases over the tested rho ranges. The half-conservative-period proxy has a
crossing in 35/36 cases; crossings are interpolated on the sampled grid.
For h=0.4, delta=0.05, its crossings versus existing subset basin optima are:

| family | fixed-time proxy crossing | observed basin optimum |
|---|---:|---:|
| generalized-alpha | 0.535 | 0.475 |
| HHT | 0.703 | 0.560 |
| WBZ | 0.824 | 0.750 |
| Newmark | 0.985 | 0.950 |

The half-period crossing changes with delta (e.g. HHT 0.703 at 0.05 to
0.837 at 0.1), and is not the proposed full-cycle D_sep. Existing optima
are from a different initial-energy population and a single h=0.4 regime;
they have no bootstrap uncertainty. No comparable four-family basin optima
at h=0.5 or h=0.6 were produced here. Therefore these comparisons cannot
establish predictive accuracy or superiority over the existing T_D median.

Reproduce with `python -m basinlab.experiments.dsep_exploratory`.
The command overwrites its output directory's named files.
