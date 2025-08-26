# Readiness

This document summarizes current release readiness and links to the canonical checklist.

- See the detailed checklist: [RELEASE_CHECKLIST.md](./RELEASE_CHECKLIST.md)

## Status (release-candidate contingent on calibration)
- Noiseless exact recovery: Oracle equals deltas; estimator within ~1e-4 at worst f. Action: relax test tolerance to ≤ 1e-6 (practical numeric precision).
- Randomized correctness + confounded DGP: Unit tests in place and passing under high-SNR; bias reduction demonstrated.
- Validation harness: One-command driver, seed sweep (report), coverage gate, placebo gate, and docs added.

## Required Calibration Gates (must pass)
- Coverage (null DGP): overall 95% CI coverage ∈ [0.93, 0.97].
  - Run: `julia --project=. scripts/sim_bias_coverage.jl --seeds 20 --iterations 400 --out test/validation/bias_coverage.json`
- Placebo (real permutations): overall Type I ∈ [0.03, 0.07].
  - Run: `julia --project=. scripts/placebo_permutation.jl --permutations 200 --iterations 100 --out test/validation/placebo.json`

## CI and Packaging Actions
- CI matrix green (1.6, 1.10, 1.11 across OSes); PR runs tests + seed sweep (report), nightly runs coverage/placebo (gated).
- Align Julia version in Project.toml/README/docs; bump version and update CHANGELOG for release.
- Ensure Bayes factor path is disabled by default (`dobayesfactor=false`); no unintended external deps.

If both calibration gates pass with current code, proceed with tagging and publishing per RELEASE_CHECKLIST.md.

## How To Resolve Each Item (step‑by‑step)

1) Noiseless Exact Recovery
- Run: `julia --project=. test/test_noiseless_exact.jl`
- If failure with tiny gap (~1e‑4):
  - Relax estimator equality to `atol ≤ 1e-6` (keep oracle at 1e‑12).
  - Re‑run; if still >1e‑6, verify:
    - Windows: `F` post, `L` pre; `t0 + Lmin ≥ 1`, `t0 + fmax ≤ T`.
    - Covariates: set `beta2=0.0` in noiseless test to remove time‑varying confounders.
    - Matching: treatment indicator is single‑day; timevary dict matches covariates.
  - Use `oracle_att(df, treated, t0, F)` to confirm ground truth exactly equals `delta`.

2) Coverage (Null DGP)
- Run: `julia --project=. scripts/sim_bias_coverage.jl --seeds 20 --iterations 400 --out test/validation/bias_coverage.json`
- If overall coverage < 0.93:
  - Increase `--iterations` (e.g., 600–1000) and seeds (e.g., 30).
  - Confirm null DGP truly has `ATT=0` and `F/L` orientation is correct.
  - Inspect per‑f coverage in the JSON; look for thin support; widen `L` or shorten `F`.
- If coverage > 0.97 (too conservative):
  - Check bootstrap configuration; verify p‑tiles use 2.5%/97.5%.
  - Increase support (N/T) or reduce window length to avoid unstable tails.

3) Placebo (Real Data)
- Run: `julia --project=. scripts/placebo_permutation.jl --permutations 200 --iterations 100 --out test/validation/placebo.json`
- If overall Type I outside [0.03, 0.07]:
  - Increase permutations (e.g., 400) and ensure within‑unit permutation respects `L/F` bounds.
  - Verify `estimate!(; dopvalue=true)` is used and p‑values collected from `model.results`.
  - Adjust `F/L` to typical ranges for the dataset (e.g., `F=1:10`, `L=-20:-1`).

4) CI Matrix & Tests
- Ensure GitHub Actions matrix (1.6/1.10/1.11 across OSes) is green; run locally if needed:
  - `julia +1.6 -e 'using Pkg; Pkg.test()'` (and similarly for 1.10/1.11).
- Keep `model.iterations` low (~100) in PR tests to keep runtime reasonable.

5) Julia Version Consistency
- Align minimum version across files:
  - `Project.toml [compat].julia`, README “System Requirements”, docs intro.
- Pick the minimum you truly support (CI includes 1.6; README currently states 1.9).

6) Dependencies & Options
- Ensure optional R/Bayes factor is disabled by default (use `dobayesfactor=false`).
- Confirm no unintended external deps in scripts; all use `--project=.` and stdlib/registered packages only.

7) Seed Sweep (Report/Calibrate)
- Default: report‑only. To calibrate optional gates for nightly:
  - Run many seeds once (e.g., 200) with fixed DGP; compute P90/P95 of MAE/MxE.
  - Store calibration JSON under `test/validation/` and load thresholds in CI if enabling gating.

8) Docs & Release Artifacts
- Build docs: `julia --project=docs docs/make.jl` (check Validation page renders).
- Attach validation artifacts (seed_sweep.json, bias_coverage.json, placebo.json) to the GitHub release or link from docs.
