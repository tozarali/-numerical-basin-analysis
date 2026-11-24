# Technical decisions

2026-11-24 — Literature search remains with the user; no web search.

2026-11-24 — Reuse the supplied skeleton modules rather than replace them.
At zero frequency or critical/overdamped physical damping do not interpret a
numerically complex parasitic pair as a physical oscillation.

2026-11-24 — For near-barrier experiments define a cycle by the next positive
x=0 reference crossing; distinguish physical capture from an unresolved
60-unit horizon. The reference first-turn-time and conservative half-period
losses are finite-time proxies, not true-cycle D_sep or validated replacements.

2026-11-24 — Fit the existing finite-time median and local linear period
error only as in-sample descriptive covariates; record dependent observations,
single-regime selection, missing uncertainty, and raw matrix conditioning.

2026-11-24 — User approved Git initialization at `paper2_idea` and plans
to run long experiments locally. Git executable was unavailable; initialized
the repository using Dulwich installed in session storage. No commit was
created. The user chose to keep NPZ results local, outside Git; existing
files are retained. Do not start a confirmatory scan
until the preregistration is complete, approved and committed.

2026-11-24 — The user chose Miniconda. Keep an isolated `basin-project`
environment instead of modifying base. The first 24x24 timing run used the
existing system Python; the launcher was subsequently tested in the conda
environment. Both timing estimates are exploratory extrapolations only. The
`run_all.ps1` launcher deliberately runs preparation only until a full
approved, resumable confirmatory runner exists.

2026-11-24 — User asked for one-shot conda execution of all experiments.
`run_all` runs exploratory experiments only by default; the confirmatory
comparison needs `--confirmatory` AND an approved, committed PREREG (code
lock). Legacy outputs that predate the runner are skipped unless --force.

2026-11-24 — WP5 FE case switched from stress-free shallow arch to a
buckled clamped-clamped beam because the unloaded arch was not bistable
(midspan reaction never returned to zero). Classification barrier = flat
compressed state, verified as index-1 by a (non-exhaustive) saddle search.

2026-11-24 — Asymmetric well and rho-infinity Bathe deliberately left out
of regime_map: not implemented/validated; Bathe coefficients must come from
the original paper, not memory.
