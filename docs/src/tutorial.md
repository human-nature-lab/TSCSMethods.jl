# Tutorial

This tutorial walks through a complete analysis using TSCSMethods.jl, explaining each step in detail.

## Complete Workflow Overview

The following diagram shows the complete user workflow from data to results:

![User Workflow](assets/images/user_workflow.svg)

Each step is explained in detail below.

## Understanding the Data Structure

TSCSMethods expects **staggered treatment design** data where:

- **Units** (e.g., counties, countries) are treated at specific times
- **Treatment** is binary (0/1) and occurs on specific dates, not continuously
- **Time periods** are relative to treatment: negative for pre-treatment, positive for post-treatment

### Example Data Format

```julia
using TSCSMethods, DataFrames

# Generate example data
dat = example_data(n_units=20, n_days=50, seed=123)
first(dat, 10)
```

Key columns:
- `fips`: Unit identifier
- `day`: Time period (0-based)  
- `gub`: Treatment indicator (1 only on treatment day)
- `death_rte`: Outcome variable
- `pop_dens`: Covariate

## Step 1: Model Creation

The `makemodel` function sets up your causal inference model:

```julia
model = makemodel(
    dat,              # Your data
    :day,             # Time variable
    :fips,            # Unit identifier  
    :gub,             # Treatment variable
    :death_rte,       # Outcome variable
    [:pop_dens],      # Covariates for matching
    Dict(:pop_dens => false),  # Which covariates are time-varying
    5:10,             # F: Post-treatment periods to estimate
    -15:-10           # L: Pre-treatment periods for matching
)
```

### Key Parameters Explained

- **F periods** (`5:10`): How many periods after treatment to estimate effects
- **L periods** (`-15:-10`): Which pre-treatment periods to use for matching
  - **Must be negative** for pre-treatment
  - Used to find similar control units
- **Time-varying covariates**: Set `true` if covariate changes over time

## Step 2: Matching

Find control units similar to treated units in pre-treatment periods:

```julia
match!(model, dat)

# Check how many matches were found
println("Found matches for $(length(model.matches)) treated observations")
```

The matching algorithm:
1. Identifies treated units and their treatment times
2. Finds control units with similar covariate values in pre-treatment periods
3. Creates matched sets for each treated observation

## Step 3: Balancing

Assess how well matching achieved covariate balance:

```julia
balance!(model, dat)

# Check balance results
model.meanbalances
```

Good balance means treated and control groups have similar covariate distributions in pre-treatment periods.

## Step 4: Estimation

Estimate average treatment effects with bootstrap inference:

```julia
# Run estimation (without Bayesian factors for simplicity)
estimate!(model, dat; dobayesfactor=false)

# View results
model.results
```

Results include:
- `f`: Time periods relative to treatment
- `att`: Average treatment effect estimates
- `q025`, `q975`: 95% confidence intervals
- `treated`: Number of treated units
- `matches`: Number of control units matched

## Interpreting Results

```julia
# Look at results
println("Treatment Effects by Time Period:")
select(model.results, :f, :att, :q025, :q975)
```

- **Positive ATT**: Treatment increased the outcome
- **Negative ATT**: Treatment decreased the outcome  
- **Confidence intervals**: Statistical uncertainty around estimates
- **Multiple time periods**: See how effects evolve over time

## Complete Example

```julia
using TSCSMethods, DataFrames

# 1. Load/generate data
dat = example_data(n_units=30, n_days=60, seed=42)

# 2. Create model  
model = makemodel(dat, :day, :fips, :gub, :death_rte, 
                 [:pop_dens], Dict(:pop_dens => false),
                 3:8, -20:-10)

# 3. Run complete workflow
match!(model, dat)
balance!(model, dat) 
estimate!(model, dat; dobayesfactor=false)

# 4. Examine results
println("Number of treated observations: ", model.treatednum)
println("Average treatment effect in period 3: ", model.results[1, :att])
```