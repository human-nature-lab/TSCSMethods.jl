# Validation Tests Overview

This suite validates statistical correctness, calibration, and robustness of TSCSMethods.jl. It follows standard practice: gate on inference calibration (coverage, Type I), and report point‑estimate error (bias/MAE/MSE) without universal raw‑error cutoffs.

## Noiseless Exact Recovery
- Files: `test/test_noiseless_exact.jl`, `test/simulate_tscs.jl`.
- Goal: Prove alignment and logic correctness by recovering injected per‑f deltas exactly.
- Method: Generate synthetic data with σy=0 and φ=0 (no noise), then run the pipeline; ATT must match `delta` up to machine precision. Also verifies oracle ATT from the simulator equals `delta`.

## Randomized Correctness (Synthetic)
- Files: `test/test_synthetic_known_effects.jl`, `test/simulate_tscs.jl`.
- Goal: Recover known per‑f ATT under a randomized, event‑time DGP with adequate support and SNR.
- Method: Build a model with `F` (post) and `L` (pre); run `match!` → `estimate!(; dobayesfactor=false)`; compare ATT to ground truth `delta`.
- Notes: High‑SNR configuration asserts small MAE/MxE; demonstrates estimator behavior with realistic noise.

## Confounded DGP — Bias Reduction
- File: `test/test_confounded_dgp.jl`.
- Goal: Show matching/calipers reduce bias when treatment depends on covariates/pre‑trends.
- Method: Generate confounded data; compare naive event‑time DiD vs TSCSMethods (optionally `autobalance`).
- Pass: Mean absolute bias reduced by ≥30% vs naive and absolute MAE ≤ 0.05.

## Coverage (Null DGP)
- Files: `test/test_coverage_type1.jl`, `scripts/sim_bias_coverage.jl`.
- Goal: Empirical 95% CI coverage near nominal under ATT=0.
- Method: Generate null DGP, run estimation with bootstrap; compute coverage using 2.5%/97.5% percentiles.
- Gate (script): overall coverage ∈ [0.93, 0.97].

## Placebo/Permutation (Real Data)
- File: `scripts/placebo_permutation.jl`.
- Goal: Type I error near 5% on permuted event times in `example_data()`.
- Method: Permute `gub` event day within unit (respecting L/F bounds); estimate ATT and collect bootstrap p‑values.
- Gate: overall Type I ∈ [0.03, 0.07].

## Seed Sweep (Reporting)
- File: `scripts/seed_sweep.jl`.
- Goal: Track MAE/MxE across seeds for a fixed randomized DGP.
- Output: JSON summary (`test/validation/seed_sweep.json`) with per‑seed errors and aggregate stats.
- Gates: Off by default (report‑only). Optional via `--gate 1` or `TSCS_GATES=1` with `--mae-thresh/--mxe-thresh`.

## One‑Command Runner
- File: `scripts/run_validations.jl`.
- Runs seed sweep (report), coverage (gate), and placebo (gate) in sequence; prints a compact summary and aggregates exit codes.

## How to Run
- Unit tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.
- All validations: `julia --project=. scripts/run_validations.jl`.
- Coverage gate: `julia --project=. scripts/sim_bias_coverage.jl --seeds 20 --iterations 400 --out test/validation/bias_coverage.json`.
- Placebo gate: `julia --project=. scripts/placebo_permutation.jl --permutations 200 --iterations 100 --out test/validation/placebo.json`.

## CI Guidance
- PR: run tests + seed sweep (report), upload artifacts.
- Nightly: run coverage and placebo (gated), optionally add seed‑sweep gates calibrated to your environment.
