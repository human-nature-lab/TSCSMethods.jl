# Validation Infrastructure

This directory contains the validation infrastructure for TSCSMethods.jl, implementing automated testing of statistical correctness across different data generating processes (DGPs).

## Structure

### Data Generating Processes (`test/simulate_tscs.jl`)
- **DGP A (Randomized)**: Unconfounded treatment assignment for baseline validation
- **DGP B (Confounded)**: Selection on observables and pre-trends for bias reduction testing  
- **DGP C (Null)**: ATT = 0 for coverage and Type I error testing

### Validation Scripts (`scripts/`)
- `seed_sweep.jl`: Randomized DGP validation across multiple seeds
- `sim_bias_coverage.jl`: Bias and coverage analysis across DGPs
- `placebo_permutation.jl`: Placebo tests using real data with permuted treatment times

### CI Integration (`.github/workflows/`)
- **PR CI**: Quick validation with seed sweep (8-10 seeds, 100 iterations)
- **Nightly CI**: Comprehensive testing (20-30 seeds, 400 iterations, all scripts)

## Quality Gates

### PR Gates (Seed Sweep)
- Mean Absolute Error (MAE) ≤ 0.03
- Maximum Error (MxE) ≤ 0.06

### Nightly Gates
- **Bias Reduction**: ≥30% reduction vs naive estimator
- **Coverage**: 95% CI empirical coverage ∈ [0.93, 0.97]  
- **Type I Error**: Placebo test Type I rate ∈ [0.03, 0.07]

## Usage

### Manual Testing
```bash
# Seed sweep with default parameters
julia --project=. scripts/seed_sweep.jl --seeds 10 --out test/validation/manual_sweep.json

# Bias and coverage analysis  
julia --project=. scripts/sim_bias_coverage.jl --seeds 15 --iterations 200 --out test/validation/manual_bias_coverage.json

# Placebo permutation tests
julia --project=. scripts/placebo_permutation.jl --permutations 100 --iterations 50 --out test/validation/manual_placebo.json
```

### CI Artifacts
- JSON files with detailed results are uploaded as GitHub Actions artifacts
- Validation summaries provide high-level quality metrics
- Failed gates cause CI to exit with non-zero status

## Output Format

All scripts generate JSON files with this general structure:
```json
{
  "config": { /* input parameters */ },
  "results": [ /* detailed per-seed/permutation results */ ],
  "summary": { /* aggregate metrics for quality gates */ }
}
```

This infrastructure ensures statistical correctness is maintained across package changes and validates the methodology against known benchmarks.