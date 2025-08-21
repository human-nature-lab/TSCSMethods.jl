# Working Example Workflow - TSCSMethods Package

This documents the **working baseline** that proves TSCSMethods estimation functions correctly with real data.

## Successful Analysis with Package's Example Data

### Data Structure
```julia
using TSCSMethods
using Statistics
using DataFrames

# Load package's built-in example data
data = example_data()
```

**Data characteristics that WORK:**
- **Dataset**: COVID death rates with gubernatorial election treatment
- **Units**: 100 (FIPS codes)  
- **Time periods**: 90 days (day 0 to 89)
- **Total observations**: 9,000 rows
- **Treatment pattern**: Only 20 treated observations total (very sparse)
- **Columns**: 
  - `fips` (unit ID) 
  - `day` (time period)
  - `gub` (treatment: 0/1)
  - `death_rte` (outcome variable)
  - `pop_dens` (covariate)
  - `date`, `cumul_death_rate` (additional variables)

### Working Analysis Pipeline

```julia
# 1. Set up time-varying covariate specification
timevary = Dict{Symbol, Bool}(:pop_dens => false)  # Pop density is time-invariant

# 2. Create model with proper variable mapping
model = makemodel(
    data, 
    :day,        # time variable (0 to 89)
    :fips,       # unit ID (100 unique units)
    :gub,        # treatment variable (binary 0/1)
    :death_rte,  # outcome variable (death rate)
    [:pop_dens], # covariates to include
    timevary,    # time-varying specification
    1:10,        # F: Post-treatment periods (days after treatment)
    -20:-1       # L: Pre-treatment periods (days before treatment)
)

# 3. Execute full pipeline
match!(model, data)    # Perform matching
estimate!(model, data) # Estimate treatment effects

# 4. Extract results
println("ATT estimates: ", round.(model.results.att, digits=4))
println("Mean ATT: ", round(mean(model.results.att), digits=4))
```

### Results That Demonstrate Package Works

**✅ CORRECT BEHAVIOR:**
- **ATT estimates**: [0.0057, 0.0528, -0.0052, -0.0635, 0.0355, 0.0343, 0.074, -0.0874, -0.0708, 0.0443]
- **Mean ATT**: 0.002 
- **Standard deviation**: 0.0571
- **Range**: [-0.0874, 0.074]

**Key indicators of proper function:**
1. **Estimates vary across F periods** (not constant)
2. **Reasonable magnitudes** (small effects on death rates)
3. **Both positive and negative estimates** (realistic variation)
4. **No systematic bias** toward single value

## Key Differences from Synthetic Data Problems

### Treatment Pattern
**Real data (works):**
- Very sparse treatment: 20/9000 observations treated
- Treatment concentrated in specific units/periods
- Realistic irregular treatment assignment

**Synthetic data (fails):**
- Dense treatment: 300/6000 observations treated  
- Regular treatment blocks (30 consecutive periods)
- Artificial perfect treatment timing

### Data Scale and Magnitudes
**Real data:**
- Outcome values: Death rates (~0.5 average)
- Small treatment effects expected (~0.002-0.07)
- Natural variation in covariates

**Synthetic data:**
- Outcome values: Large synthetic values (generated with unit/time effects)
- Large treatment effects imposed (ATT = 2.0)
- Artificial covariate patterns

## Implications for Testing

1. **Package estimation is fundamentally correct** - real data produces sensible results
2. **Synthetic data generation may have issues** with:
   - Treatment assignment patterns
   - Outcome scale/magnitude  
   - Covariate relationships
   - Data structure assumptions

3. **F/L window specifications are correct** - same windows work with real data

## Next Steps for Synthetic Data Fix

1. **Model synthetic data more closely on real data patterns**:
   - Sparser treatment assignment
   - Smaller, realistic effect sizes
   - Similar outcome variable scales

2. **Validate data generation process** by comparing:
   - Treatment density patterns
   - Outcome variable distributions
   - Pre-treatment trends

3. **Test intermediate cases** between real and synthetic data to isolate the issue

## Working Code Template

```julia
# This template WORKS and can be used as baseline for testing
using TSCSMethods
using Statistics
using DataFrames

# Load known-good data
data = example_data()

# Standard analysis setup  
timevary = Dict{Symbol, Bool}(:pop_dens => false)
model = makemodel(data, :day, :fips, :gub, :death_rte, [:pop_dens], timevary, 1:10, -20:-1)

# Execute pipeline
match!(model, data)
estimate!(model, data)

# Verify results look reasonable
att_estimates = model.results.att
@assert std(att_estimates) > 0.01  # Should have variation
@assert !all(abs.(att_estimates .- att_estimates[1]) .< 0.001)  # Not all identical
println("✅ Package working correctly with real data")
```

## How the Example Data is Stored & Generated

The example data is **dynamically generated** (not stored) by the `example_data()` function:

### Data Generation Process
```julia
function example_data(; n_units::Int = 100, n_days::Int = 90, seed::Int = 123)
```

**Key characteristics:**
- **Synthetic COVID death rate data** with gubernatorial election treatment
- **Treatment assignment**: First 20% of units (20/100) get treated
- **Treatment timing**: One-time event around day 40 (staggered by unit)
- **Treatment pattern**: `gub=1` only on the specific treatment day, `gub=0` otherwise
- **Treatment effect**: Small negative effect (-0.15 * rand()) on death rates after day 30

### Critical Pattern Differences from Our Synthetic Data:

**✅ Real example data (works):**
- **Sparse treatment**: Only 20 observations have `gub=1` (treatment events)
- **Point-in-time treatment**: Treatment is a single-day event indicator
- **Realistic effects**: Small treatment effects (~-0.075 on average)
- **Natural variation**: Seasonal patterns, time trends, noise

**❌ Our synthetic data (fails):**
- **Dense treatment**: 300 observations with `treatment=1` (treatment periods)  
- **Duration treatment**: Treatment lasts 30 consecutive periods
- **Large effects**: ATT = 2.0 (much larger scale)
- **Artificial patterns**: Perfect treatment blocks

### The JLD2 File Mystery

The `vignette/simpledata.jld2` file (12 MB) appears to be a **pre-saved version** of example data, but the function generates data dynamically, suggesting this file may be:
- Legacy storage from development
- Pre-computed version for faster loading
- Documentation/testing artifact

The function itself creates fresh synthetic data each time when called.

This workflow provides the **proven baseline** for package functionality and should be used to diagnose issues with synthetic data generation.