# Validation Implementation Spec

Purpose: Provide precise specifications for automation (scripts, tests, CI) to validate TSCSMethods.jl. This document is for implementation by a coding agent.

## 1) Data Generating Processes (DGPs)

Common parameters
- Inputs: `N::Int`, `T::Int`, `treated_share::Float64`, `F::UnitRange{Int}`, `L::UnitRange{Int}`, `phi::Float64` (AR(1) outcome), `rho::Float64` (AR covariate), `sigma_y::Float64`, `beta1::Float64`, `beta2::Float64`, `seed::Int`.
- Output DataFrame columns: `:id, :t, :gub, :y, :x1, :x2`.
- Event-time treatment: one-day indicator `gub=1` at `t0[id]`; 0 otherwise. Ensure full window bounds: `t0 + minimum(L) ≥ 1` and `t0 + maximum(F) ≤ T`.

DGP A — Randomized (Unconfounded)
- Treatment: Uniform random subset of units with share `treated_share`; event times uniform in valid range.
- Outcome baseline: `y0 = unit_FE + time_FE + beta1*x1 + beta2*x2 + phi*y_{t-1} + noise`.
- Known effects: apply `delta[f_index]` at `t = t0 + f` for each `f ∈ F`.

DGP B — Confounded (Selection on observables and pre-trend)
- Treatment assignment: logit probability increases with `x1`, `x2`, and a pre-trend slope `pretrend[id]`. Draw treated ids by probability; draw `t0` similarly with weak dependence on `x1`.
- Pre-trend: add unit-specific slope to baseline outcome prior to treatment: `y0 += pretrend[id]*(t - t0[id])` for `t < t0[id]`.
- Apply known effects as in DGP A.

DGP C — Null (ATT = 0)
- Same as DGP A but `delta .= 0.0`.

Implementation file: `test/simulate_tscs.jl`
- Function signature:
  - `simulate_tscs(; N, T, treated_share, F, L, delta, phi, rho, sigma_y, beta1, beta2, seed, confounded::Bool=false, conf_strength::Float64=0.7, pretrend::Bool=false, pretrend_strength::Float64=0.02)` -> `(DataFrame, Vector{Int}, Dict{Int,Int})`
- Provide thin wrappers: `simulate_randomized(...)`, `simulate_confounded(...)`, `simulate_null(...)` using the same engine.

## 2) Tests

A) Randomized DGP (already present)
- File: `test/test_synthetic_known_effects.jl` (keep as-is with current parameters).
- Ensure runtime control by setting `model.iterations = 100` before `estimate!`.

B) Confounded DGP — Bias reduction with matching
- File: `test/test_confounded_dgp.jl`
- Steps:
  1) Generate confounded data (DGP B) with small true `delta` (e.g., peak 0.05) and moderate confounding.
  2) Compute naive event-time DiD ATT (oracle baseline under randomization; here it will be biased).
  3) Run TSCS pipeline:
     - `timevary = Dict(:x1=>false, :x2=>true)`
     - `makemodel(...); match!(...);` optionally `autobalance(...; threshold=0.1, refinementnum=3, calmin=0.1, step=0.05)`
     - `estimate!(...; dobayesfactor=false)`.
  4) Metrics: per-f absolute bias vs truth for naive and matched.
  5) Assertions (PR-safe thresholds):
     - Mean absolute bias reduced by ≥30% vs naive.
     - Post-matching mean absolute bias ≤ 0.05.
- Keep runtime ≤ ~60s: N≈600, T≈200, treated_share≈0.10, `iterations=100`.

C) Coverage and Type I (optional to start)
- File: `test/test_coverage_type1.jl`
- DGP C with ATT=0, run `estimate!` with `iterations=400` and compute empirical coverage for 95% CIs across 10 seeds. Assert coverage ∈ [0.90, 0.99] for PR (tighten in nightly).

Wire in runtests
- Add includes under existing `@testset "Correctness tests"`.

## 3) Scripts

A) Seed Sweep — Randomized DGP
- File: `scripts/seed_sweep.jl`
- CLI: `--seeds 10 --N 600 --T 200 --treated-share 0.10 --F 1:8 --L -30:-1 --phi 0.0 --rho 0.5 --sigma-y 0.2 --beta1 0.3 --beta2 0.5 --out test/validation/seed_sweep.json`
- Behavior: for each seed, generate DGP A, run TSCS with `iterations=100`, compute MAE, MxE, and per-f errors. Emit JSON:
```
{
  "config": {...},
  "results": [
     {"seed": 123, "mae": 0.023, "mxe": 0.041, "per_f": [ ... ]},
     ...
  ],
  "summary": {"mae_mean": ..., "mae_sd": ..., "mxe_max": ...}
}
```
- Exit nonzero if `mae_mean > 0.03` or `mxe_max > 0.06` (PR gates).

B) Bias + Coverage — Randomized and Null
- File: `scripts/sim_bias_coverage.jl`
- Runs: (i) DGP A with known `delta`, reports bias/RMSE; (ii) DGP C (null), reports 95% coverage across seeds.
- JSON output: `test/validation/bias_coverage.json` with sections `bias` and `coverage`.

C) Placebo/Permutation — Real Data
- File: `scripts/placebo_permutation.jl`
- Load `example_data()`. For K permutations (default 200), randomly permute `gub` event times within unit, run TSCS (iterations=100), record per-f ATTs; compute two-sided Type I error against 0. Output `test/validation/placebo.json`. Exit nonzero if Type I not in [0.03, 0.07].

Implementation details
- All scripts must run with `julia --project=.`.
- Set `model.iterations` before `estimate!` to control runtime.
- Always call `estimate!(; dobayesfactor=false)`.
- Ensure deterministic runs with the provided seeds.

## 4) CI Integration

PR CI (quick)
- Run package tests (existing).
- Run `scripts/seed_sweep.jl` with 6–10 seeds and low iterations; upload JSON artifacts.

Nightly CI (full)
- Seed sweep with more seeds (e.g., 30) and `iterations=400`.
- Run `scripts/sim_bias_coverage.jl`.
- Optionally run `scripts/placebo_permutation.jl`.

Gates
- Randomized DGP: `mae_mean ≤ 0.03`, `mxe_max ≤ 0.06`.
- Coverage: 95% in [0.93, 0.97].
- Placebo Type I: in [0.03, 0.07].

## 5) Examples

- Run seed sweep:
```
julia --project=. scripts/seed_sweep.jl --seeds 10 --N 600 --T 200 --treated-share 0.10 \
  --F 1:8 --L -30:-1 --phi 0.0 --rho 0.5 --sigma-y 0.2 --beta1 0.3 --beta2 0.5 \
  --out test/validation/seed_sweep.json
```

- Run bias + coverage:
```
julia --project=. scripts/sim_bias_coverage.jl --seeds 20 --iterations 400 \
  --out test/validation/bias_coverage.json
```

## 6) Guardrails
- Do not add R dependencies; always call `estimate!(; dobayesfactor=false)`.
- Respect `.JuliaFormatter.toml` (2-space indent, 92 margin).
- Keep scripts idempotent; validate inputs; fail with clear messages.
- Keep tests deterministic and under PR runtime budgets (≤ 10 min overall).
