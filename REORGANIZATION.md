# TSCSMethods.jl Reorganization Plan

## Overview

This document outlines the reorganization plan for TSCSMethods.jl to improve maintainability and developer experience while preserving all existing functionality and dependencies.

## Current Structure Issues

- **Flat Source Directory**: 38+ files in `src/` with no logical grouping
- **Mixed Test Organization**: Unit tests and integration tests intermixed
- **Scattered Validation**: Scripts in `scripts/` directory separate from validation results
- **Unclear Dependencies**: Dependency flow not obvious from file structure

## Proposed New Structure

### Source Code Organization (`src/`)

```
src/
├── TSCSMethods.jl           # Main module (updated funlist paths)
├── core/                    # Core types and construction
│   ├── types.jl
│   └── construction.jl
├── matching/                # Matching subsystem
│   ├── matching_setup.jl
│   ├── retrieve_matches.jl
│   ├── retrieve_matches_missing.jl
│   ├── retrieve_matches_utilities.jl
│   ├── distancing.jl
│   ├── distancing_utilities.jl
│   ├── match.jl
│   ├── ranking.jl
│   └── caliper.jl
├── balancing/               # Balancing subsystem
│   ├── balancing.jl
│   ├── balancing_utilities.jl
│   ├── meanbalancing.jl
│   ├── overallbalancing.jl
│   ├── fullbalancing.jl
│   ├── mean_fullbalancing.jl
│   └── autobalancing.jl
├── estimation/              # Estimation subsystem
│   ├── estimation.jl
│   ├── estimation_setup.jl
│   ├── estimation_utilities.jl
│   ├── estimation_observationweights.jl
│   ├── estimation_stratified.jl
│   ├── overall.jl
│   ├── bayesfactor.jl
│   ├── bootstrapping.jl
│   └── resampling.jl
├── advanced/                # Advanced features
│   ├── stratification.jl
│   ├── refine.jl
│   └── groupindices.jl
└── utilities/               # Utility functions
    ├── model_utilities.jl
    ├── storage.jl
    ├── information.jl
    ├── imputation.jl
    ├── inspection.jl
    ├── filterunits!.jl
    └── example_data.jl
```

### Test Structure (`test/`)

```
test/
├── runtests.jl              # Main test runner (updated includes)
├── unit/                    # Unit tests by subsystem
│   ├── test_construction.jl
│   ├── test_matching.jl
│   ├── test_distance.jl
│   ├── test_balancing.jl
│   ├── test_estimation.jl
│   └── test_utilities.jl
├── integration/             # Integration and workflow tests
│   ├── test_basic_workflow.jl
│   ├── test_workflow.jl
│   └── test_simple.jl
├── correctness/             # Statistical correctness tests
│   ├── test_noiseless_exact.jl
│   ├── test_synthetic_known_effects.jl
│   ├── test_confounded_dgp.jl
│   ├── test_coverage_type1.jl
│   └── test_statistical_correctness.jl
├── validation/              # Validation framework (existing)
│   └── [existing validation files...]
├── benchmark/               # Performance tests
│   ├── benchmark_correctness.jl
│   └── baseline_benchmark.json
└── support/                 # Test utilities
    ├── simulate_tscs.jl
    ├── test_input_validation.jl
    └── test_results.md
```

### Other Directories

```
validation/                  # Move from scripts/ to top-level
├── compatibility_test.jl
├── placebo_permutation.jl
├── run_validations.jl
├── seed_sweep.jl
└── sim_bias_coverage.jl

examples/                    # Rename from vignette/
├── enhanced_tutorial.md
├── example_data.csv
├── example_data.jls
├── simpledata.jld2
└── tutorial.ipynb

docs/                        # Keep existing structure
├── [existing structure unchanged]
```

## Updated funlist Array

The critical `funlist` array in `TSCSMethods.jl` will be updated to reflect new paths while maintaining **exact same dependency order**:

```julia
funlist = [
    "core/types.jl",                           # 1: Base types
    "core/construction.jl",                    # 2: Model construction
    "matching/matching_setup.jl",              # 3: Matching setup
    "advanced/groupindices.jl",                # 4: Group utilities
    "matching/retrieve_matches_utilities.jl",   # 5: Match retrieval utilities
    "matching/retrieve_matches.jl",            # 6: Match retrieval
    "matching/retrieve_matches_missing.jl",     # 7: Missing data handling
    "matching/distancing_utilities.jl",        # 8: Distance utilities
    "matching/distancing.jl",                  # 9: Distance calculations
    "matching/match.jl",                       # 10: Core matching
    "matching/ranking.jl",                     # 11: Match ranking
    "matching/caliper.jl",                     # 12: Caliper matching
    "balancing/meanbalancing.jl",              # 13: Mean balancing
    "balancing/balancing_utilities.jl",        # 14: Balancing utilities
    "balancing/overallbalancing.jl",           # 15: Overall balancing
    "balancing/balancing.jl",                  # 16: Main balancing
    "estimation/estimation_setup.jl",          # 17: Estimation setup
    "estimation/estimation_utilities.jl",      # 18: Estimation utilities
    "estimation/estimation_observationweights.jl", # 19: Observation weights
    "estimation/estimation.jl",                # 20: Main estimation
    "estimation/estimation_stratified.jl",     # 21: Stratified estimation
    "estimation/overall.jl",                   # 22: Overall results
    "estimation/bayesfactor.jl",               # 23: Bayesian factors
    "estimation/resampling.jl",                # 24: Resampling
    "estimation/bootstrapping.jl",             # 25: Bootstrap methods
    "advanced/stratification.jl",              # 26: Stratification
    "advanced/refine.jl",                      # 27: Refinement
    "balancing/autobalancing.jl",              # 28: Automatic balancing
    "utilities/model_utilities.jl",            # 29: Model utilities
    "utilities/storage.jl",                    # 30: Storage utilities
    "utilities/information.jl",                # 31: Information utilities
    "utilities/imputation.jl",                 # 32: Imputation
    "utilities/inspection.jl",                 # 33: Inspection utilities
    "utilities/filterunits!.jl",               # 34: Unit filtering
    "utilities/example_data.jl"                # 35: Example data
];
```

## Implementation Strategy

### Phase 1: Preparation
- [ ] Create new directory structure (empty folders)
- [ ] Verify all tests pass in current state
- [ ] Create backup branch: `git checkout -b backup-before-reorganization`

### Phase 2: File Movement
- [ ] Use `git mv` commands to preserve file history:
  ```bash
  # Core files
  mkdir -p src/core
  git mv src/types.jl src/core/
  git mv src/construction.jl src/core/
  
  # Matching files
  mkdir -p src/matching
  git mv src/matching_setup.jl src/matching/
  git mv src/retrieve_matches*.jl src/matching/
  git mv src/distancing*.jl src/matching/
  git mv src/match.jl src/matching/
  git mv src/ranking.jl src/matching/
  git mv src/caliper.jl src/matching/
  
  # Continue for all subsystems...
  ```

### Phase 3: Update Main Module
- [ ] Update `funlist` array in `src/TSCSMethods.jl` with new paths
- [ ] Ensure exact same include order is maintained

### Phase 4: Update Test Structure
- [ ] Create new test subdirectories
- [ ] Move test files to appropriate locations using `git mv`
- [ ] Update `test/runtests.jl` with new include paths:
  ```julia
  @testset "Backend tests" begin
      include("unit/test_simple.jl")
      include("unit/test_construction.jl")
      include("unit/test_matching.jl")
      include("unit/test_distance.jl")
      include("unit/test_balancing.jl")
      include("unit/test_estimation.jl")
      include("unit/test_utilities.jl")
  end

  @testset "Correctness tests" begin
      include("correctness/test_noiseless_exact.jl")
      include("correctness/test_synthetic_known_effects.jl")
      include("correctness/test_confounded_dgp.jl")
      include("correctness/test_coverage_type1.jl")
  end
  ```

### Phase 5: Move Supporting Directories
- [ ] `git mv scripts/ validation/`
- [ ] `git mv vignette/ examples/`
- [ ] Update any references to old paths in documentation

### Phase 6: Validation
- [ ] Run full test suite: `julia --project=. -e "using Pkg; Pkg.test()"`
- [ ] Run validation scripts in new location
- [ ] Verify documentation builds: `julia --project=docs/ docs/make.jl`
- [ ] Check all examples work correctly

### Phase 7: Documentation Updates
- [ ] Update `CLAUDE.md` with new structure information
- [ ] Update any hardcoded paths in documentation files
- [ ] Update development commands if needed

### Phase 8: Final Verification
- [ ] Clean test run on fresh Julia environment
- [ ] Verify package can be imported and basic workflow works
- [ ] Run comprehensive validation suite

## Benefits of This Reorganization

### Developer Experience
- **Logical Grouping**: Related functionality is co-located
- **Clear Navigation**: Easier to find relevant files
- **Scalable Structure**: New features can be added to appropriate subsystems
- **Better Testing**: Test structure mirrors source organization

### Maintainability
- **Dependency Clarity**: Flow from core → matching → balancing → estimation is obvious
- **Modular Organization**: Each subsystem is self-contained
- **Consistent Patterns**: Similar naming and organization across subsystems

### Backward Compatibility
- **No API Changes**: All exported functions remain identical
- **Preserved Dependencies**: Include order maintained exactly
- **Git History**: File history preserved through `git mv`

## Risk Mitigation

### Dependency Preservation
- Maintain exact same `funlist` order to preserve include dependencies
- Use `git mv` to preserve file history
- Test after each major phase

### Rollback Plan
- Backup branch created before any changes
- Each phase can be individually reverted if issues arise
- Full test suite run after each phase

### Validation
- Comprehensive test coverage ensures functionality is preserved
- Multiple validation checkpoints throughout process
- Fresh environment testing to catch path issues

## Post-Reorganization

After successful reorganization:
- Update this document with "COMPLETED" status
- Archive this document or move to `docs/development/`
- Update contributor documentation with new structure
- Consider updating package version to reflect structural improvement

## Notes

- This reorganization maintains 100% backward compatibility
- All existing code using TSCSMethods.jl will continue to work unchanged
- Internal organization improvement only - no functional changes
- Careful preservation of the critical `funlist` dependency order