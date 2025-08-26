# Validation.md
_via codex_

The (initial) evidence supports that TSCSMethods.jl is estimating correctly.

## Conclusion

ATTs, the package recovers effects with low error (mae ≈ 0.023, mxe ≈ 0.040), consistent with the naive event-time DiD oracle.
- Root cause of earlier failures: Low SNR, insufficient treated support per f, and windowing constraints — not structural bias in the estimator.
- Assumptions satisfied: Single-day treatment indicator, proper L/F window alignment, adequate support, and no confounding. Matching is optional here and does not drive correctness.

## What this does not yet prove

- Confounding: Performance when treatment depends on covariates or pre-trends.
- Alternative patterns: Dense/duration treatment, different AR dynamics,
missingness.
- Robustness: Sensitivity to L/F choices and parameter grids.

## Fast next validations

- Confounded DGP: Make treatment depend on x1/x2 or pre-trend; verify raw bias and show matching/balancing (calipers/refinement) reduce it.
- Sensitivity: Sweep L/F and treated density; confirm stability of estimates.
- Placebo/permutation (real data): Randomize event times within unit; Type I error near 5%.
- Coverage: Bootstrap CI coverage under a null DGP.

If you’d like, I can wire the passing synthetic test into runtests.jl
(short version for PRs, full for nightly) and add a confounded DGP test to demonstrate bias reduction with matching.
