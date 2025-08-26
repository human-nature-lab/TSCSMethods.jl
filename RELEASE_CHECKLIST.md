# Release Checklist — TSCSMethods.jl

Use this checklist to validate statistical correctness, finalize packaging, and publish a tagged release. Items marked Gate must pass before release.

## Statistical Validation
- Noiseless exact recovery (sanity): Gate
  - Run: `julia --project=. test/test_noiseless_exact.jl`
  - Expect: oracle ATT equals `delta` exactly; estimator ATT within atol ≤ 1e-6 of `delta`.
- Coverage (null DGP): Gate
  - Run: `julia --project=. scripts/sim_bias_coverage.jl --seeds 20 --iterations 400 --out test/validation/bias_coverage.json`
  - Expect: overall 95% CI coverage ∈ [0.93, 0.97].
- Placebo/Permutation (real data): Gate
  - Run: `julia --project=. scripts/placebo_permutation.jl --permutations 200 --iterations 100 --out test/validation/placebo.json`
  - Expect: overall Type I error ∈ [0.03, 0.07].
- Seed sweep (randomized DGP): Report-only
  - Run: `julia --project=. scripts/seed_sweep.jl --out test/validation/seed_sweep.json`
  - Action: Attach JSON to release or link from docs; do not gate by default.

## Tests & CI
- Unit tests: Gate
  - Run: `julia --project=. -e 'using Pkg; Pkg.test()'`
  - Ensure new tests pass (noiseless, randomized correctness, confounded DGP, coverage unit).
- CI matrix: Gate
  - Green on GitHub Actions for Julia 1.6, 1.10, 1.11 across OSes.
- Validation workflow:
  - PR: runs tests + seed sweep (report), uploads JSON artifacts.
  - Nightly: runs coverage and placebo (gated), uploads artifacts.

## Packaging & Metadata
- Versioning: Gate
  - Bump `Project.toml` version (semantic versioning); add `CHANGELOG.md` entry.
- Julia version consistency: Gate
  - Align minimum Julia version across `Project.toml`, README, and docs (CI currently tests 1.6/1.10/1.11).
- Dependencies: Gate
  - Ensure optional R/Bayes factor path is disabled by default (`dobayesfactor=false`); no unintended external deps.

## Quality & Docs
- Formatting: Gate
  - Run formatter check (GitHub Actions Format Check); fix any diffs.
- Documentation: Gate
  - Docs build passes: `julia --project=docs docs/make.jl`.
  - Validation docs present: `docs/src/validation.md`; README links included.
  - Add brief release notes summarizing validation results and any API changes.

## Tag & Publish
- Create a git tag (e.g., `vX.Y.Z`) and push; verify CI and docs deploy green.
- Attach or link validation artifacts (seed_sweep.json, bias_coverage.json, placebo.json) in the release.
- Announce version, validation status, and highlights in release notes.

Notes
- If any Gate fails, fix or adjust configuration and re-run before tagging.
- Retain validation artifacts for traceability; they demonstrate calibration at release time.
