# Phase 0 — literature search protocol

Status: **open; not a systematic review and no novelty claim established**.
The exploratory web searches on 2026-11-24 did not substitute for Scopus,
Web of Science, arXiv and forward/backward citation searches. Record exact
search date, platform, complete query, filters, hit count, export, duplicates,
screening decisions and reasons for exclusion before evaluating the gate.
No database hit counts are claimed here.

## Search strings

Adapt syntax to each database; search titles, abstracts and keywords without
year restriction, then document any language/document-type filters:

1. `("numerical damping" OR "algorithmic dissipation" OR "numerical dissipation" OR "time integration") AND ("basin of attraction" OR "basin boundary" OR multistab* OR bistab* OR "equilibrium selection")`
2. `(Newmark OR "generalized-alpha" OR "generalised-alpha" OR HHT OR "Hilber-Hughes-Taylor" OR "Wilson-theta" OR "Bathe method") AND (basin* OR multistab* OR bistab* OR "snap-through")`
3. `("basin entropy" OR "basin stability" OR "uncertainty exponent") AND ("numerical integrat*" OR discretization OR discretisation OR "time step")`
4. `(Duffing OR "double well" OR "shallow arch" OR "von Mises truss") AND ("numerical damping" OR "time integration" OR "numerical dissipation") AND (basin* OR attractor* OR "snap-through")`

Include numerical studies that compare basin or equilibrium outcomes across
integrators, time steps or numerical dissipation, even when their terminology
differs. Retain structural-dynamics integrator papers and basin-method papers
as contextual records. Exclude unrelated static optimization basins and
physical-damping-only studies from the direct-overlap category, but log
borderline cases. Screen title/abstract, then full text; snowball cited and
citing work for each near match. Keep bibliographic exports and a screening
table with DOI/URL, method, model, outcome, relationship to proposed claim and
decision.

## Verified starting points, not a complete search

* Daza et al., *Basin entropy: a new tool to analyze uncertainty in dynamical
  systems* (2016): https://arxiv.org/abs/1605.02342 . Defines Sb and Sbb;
  Sbb > ln(2) is sufficient, not necessary, for a fractal boundary, and cannot
  be satisfied by an exact two-outcome entropy estimate.
* Menck et al., *How basin stability complements the linear-stability
  paradigm* (2013): https://www.nature.com/articles/nphys2516 .
* *Comparison of numerically dissipative schemes for structural dynamics:
  Generalized-alpha versus energy-decaying methods* (2020):
  https://www.sciencedirect.com/science/article/pii/S0263823120309496 .
* *For direct time integrations: A comparison of the Newmark and
  rho-infinity-Bathe schemes* (2019):
  https://www.sciencedirect.com/science/article/abs/pii/S0045794919304304 .

Gate: **not passed**. Full database coverage and screening are pending; a
direct overlap would require a positioning/differentiation analysis rather
than an automatic assertion of novelty.
