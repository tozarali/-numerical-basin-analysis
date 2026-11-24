# PREREG — <experiment name>   (commit BEFORE running; then do not edit)

- Date / commit of code used:
- Question & hypothesis (H#):
- Models and parameters: family list, rho grid (fine step near expected optimum), h values, c values, epsilon
- IC population: generator (sobol_box seed=…, N=…), box, energy band (E_b, k*E_b], primary vs secondary population
- Reference: method, step/tolerance, T_basin; check (T vs 2T, h vs h/2) on the same ICs
- Primary measure: mismatch on primary population. Secondary: Sb/Sbb (box sizes …), net transfer, signed shift
- Saturation rule: cells with mismatch >= 30 % OR Sb >= 0.9 ln2 are flagged and excluded from correlations
- Optimum estimator: stats.bootstrap_optimum(half_width=…, n_boot=…); report discrete argmin as well
- D definition: (T_D median | D_sep delta=… | D_th level …); crossing by dtheory.crossing
- Success criterion (gate): e.g. |rho*_pred - rho*_obs| <= 0.1 in >= 80 % of reliable cells AND better than nominal-rho baseline
- What would count as failure, and what we will write then:
- Expected runtime / memory:
