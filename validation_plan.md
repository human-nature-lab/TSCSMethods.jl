# Validation Plan for TSCSMethods.jl

## Objectives
- Validate estimator correctness under noiseless and noisy synthetic scenarios.
- Emphasize standard inference calibration: 95% coverage and 5% Type I.
- Detect regressions via CI gates on calibration; report point-error metrics without gating by default.

## Immediate Steps (short-term)
- Noiseless synthetic (exact recovery):
  - DGP: randomized event-time; σy=0, φ=0; constant per‑f deltas; adequate L/F.
  - Expectation: `ATT_hat[f] == delta[f]` (machine precision). Fix windowing/matching if not.
- Randomized high‑SNR correctness (unit test):
  - Keep current test (per‑f ATT recovery) with modest noise and healthy support.
  - Reduce bootstrap iterations for speed (e.g., `model.iterations = 100`).
- Confounded DGP test:
  - Treatment probability increases with `x1`, `x2`, and pre‑trend.
  - Show naive bias; verify ≥30% MAE reduction with matching/calipers and MAE ≤ 0.05.

## Medium Steps
- Coverage (null DGP): 95% CI empirical coverage near nominal; gate overall coverage ∈ [0.93, 0.97].
- Placebo/permutation (real data): Type I near 5% with ≥200 permutations; gate overall Type I ∈ [0.03, 0.07].
- Sensitivity sweeps: Grid over `L`, `F`, treated share, φ; confirm stability and adequate support per f.
- Cross‑method baselines: TWFE DiD and event‑study OLS as sanity checks on simulations.

## Automation & CI
- Scripts: `scripts/run_validations.jl` (driver), `scripts/sim_bias_coverage.jl`, `scripts/seed_sweep.jl`, `scripts/placebo_permutation.jl` (artifacts in `test/validation/`).
- CI jobs:
  - PR: unit tests + seed sweep (report-only); upload artifacts.
  - Nightly: coverage (gate) + placebo (gate); optionally enable calibrated seed‑sweep gates.
- Gates (standard practice):
  - Coverage: overall ∈ [0.93, 0.97].
  - Placebo Type I: overall ∈ [0.03, 0.07].
  - Point-error (MAE/MxE): report-only by default; gated only with calibrated thresholds.

## Reporting
- Persist metrics (JSON/CSV) and small plots as CI artifacts.
- Summarize “last good” calibration results in docs; keep seed-sweep reports for trend monitoring.

## Timeline
- Week 1: Noiseless exact-recovery + high‑SNR correctness + confounded test.
- Week 2: Coverage and placebo gates + sensitivity sweeps.
- Week 3: CI wiring (PR vs nightly), documentation polish, optional calibrated gates for seed‑sweep.
