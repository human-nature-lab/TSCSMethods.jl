# TSCSMethods.jl: Complete Tutorial

This tutorial provides a comprehensive walkthrough of causal inference with time-series cross-sectional data using TSCSMethods.jl.

## Background: What is TSCS Causal Inference?

Time-series cross-sectional (TSCS) data combines:
- **Cross-sectional variation**: Different units (counties, countries, etc.)
- **Time-series variation**: Observations over multiple time periods
- **Treatment variation**: Units receive treatment at different times (**staggered design**)

### Example Research Questions
- Do gubernatorial elections affect COVID-19 outcomes?
- How do policy changes impact economic indicators?
- What is the effect of protests on subsequent political behavior?

## Tutorial Dataset

We'll use a simulated dataset based on the COVID-19 analysis from Feltham et al. (2023):

```julia
using TSCSMethods, DataFrames, Statistics

# Generate realistic example data
dat = example_data(n_units=50, n_days=80, seed=42)

# Examine the data structure
first(dat, 10)
```

### Key Variables
- `fips`: County identifier (units)
- `day`: Time period (0-based from start date)
- `gub`: Treatment indicator (gubernatorial election = 1)
- `death_rte`: Outcome variable (COVID death rate)
- `pop_dens`: Population density (time-invariant covariate)
- `cumul_death_rate`: Cumulative death rate (time-varying covariate)

### Understanding Treatment Structure

```julia
# Check treatment pattern - should be event-based, not continuous
treatment_summary = combine(groupby(dat, :fips), 
    :gub => sum => :total_treatments,
    :gub => (x -> findall(==(1), x)) => :treatment_days
)

# Show units with treatment
filter(:total_treatments => >(0), treatment_summary)
```

**Key insight**: Each treated unit has exactly one treatment event (gub=1) on a specific day, not continuous treatment.

## Step 1: Understanding Time Periods

TSCSMethods uses **relative time** periods:

```julia
# Example: If a unit is treated on day 30
treatment_day = 30

# L periods (pre-treatment): -10:-1 means days 20-29
L_periods = -10:-1
actual_L_days = treatment_day .+ L_periods  # [20, 21, ..., 29]

# F periods (post-treatment): 1:5 means days 31-35  
F_periods = 1:5
actual_F_days = treatment_day .+ F_periods  # [31, 32, 33, 34, 35]

println("Pre-treatment matching window: days ", actual_L_days)
println("Post-treatment estimation window: days ", actual_F_days)
```

## Step 2: Model Creation

```julia
# Define matching covariates
matching_covariates = [:pop_dens, :cumul_death_rate]

# Specify which covariates vary over time
timevary = Dict(
    :pop_dens => false,        # Population density is roughly constant
    :cumul_death_rate => true  # Cumulative death rate changes over time
)

# Create the model
model = makemodel(
    dat,                    # Data
    :day,                   # Time variable  
    :fips,                  # Unit identifier
    :gub,                   # Treatment variable
    :death_rte,             # Outcome variable
    matching_covariates,    # Covariates for matching
    timevary,              # Time-varying specification
    3:8,                   # F: Estimate effects 3-8 periods after treatment
    -15:-5;                # L: Match on periods 15-5 before treatment
    title = "COVID_Elections"
)

println("Created model with $(model.treatednum) treated observations")
```

## Step 3: Matching

Find control units similar to treated units in pre-treatment periods:

```julia
# Perform matching
match!(model, dat)

# Check matching results
println("Matching completed:")
println("- Treated observations: $(length(model.observations))")
println("- Total potential controls: $(length(model.ids))")
println("- Matches found: $(length(model.matches))")

# Examine a single match
if length(model.matches) > 0
    first_match = model.matches[1]
    println("First treated unit has $(sum(first_match.mus)) total matches")
end
```

### How Matching Works

1. **For each treated unit at treatment time τ**:
   - Look at covariate values in pre-treatment periods (L)
   - Find control units with similar covariate values in same periods
   - Calculate distances and select closest matches

2. **Distance calculation**:
   - Standardized by variance among treated units
   - Weighted by time period importance

## Step 4: Balance Assessment  

Check if matching achieved good covariate balance:

```julia
# Calculate balance statistics
balance!(model, dat)

# Examine balance results
println("Balance assessment completed")
println("Mean balances structure: ", size(model.meanbalances))

# Check overall balance
if haskey(model, :grandbalances)
    println("Grand balances available for covariates: ", keys(model.grandbalances))
end
```

### Interpreting Balance
- **Small balance statistics**: Good match between treated and control groups
- **Large balance statistics**: Poor match, may need refinement or calipers

## Step 5: Treatment Effect Estimation

Estimate average treatment effects with bootstrap confidence intervals:

```julia
# Perform estimation (without Bayesian factors for simplicity)
estimate!(model, dat; dobayesfactor=false)

# Examine results
println("\\nTreatment Effect Results:")
println("Periods estimated: ", model.results.f)
println("Sample ATT estimates:")
select(model.results, :f, :att, :q025, :q975, :treated)
```

### Understanding Results

```julia
# Detailed interpretation
results = model.results
for row in eachrow(results)
    period = row.f
    att = round(row.att, digits=4)
    ci_lower = round(row.q025, digits=4)  
    ci_upper = round(row.q975, digits=4)
    
    significance = (ci_lower > 0 && ci_upper > 0) || (ci_lower < 0 && ci_upper < 0)
    sig_marker = significance ? " **" : ""
    
    println("Period $period: ATT = $att [$ci_lower, $ci_upper]$sig_marker")
end
```

## Step 6: Advanced Features

### Refinement: Improve Match Quality

```julia
# Refine to best 10 matches per treated unit
refined_model = refine(
    model, dat;
    refinementnum = 10,
    dobalance = true,
    doestimate = true
)

println("Refined model results:")
select(refined_model.results, :f, :att, :q025, :q975)
```

### Calipers: Restrict Match Quality

```julia
# Apply calipers to ensure minimum match quality
initial_calipers = Dict(:cumul_death_rate => 0.5)

caliper_model = caliper(model, dat, initial_calipers)
balance!(caliper_model, dat)
estimate!(caliper_model, dat; dobayesfactor=false)

println("Caliper model with stricter matching:")
select(caliper_model.results, :f, :att, :q025, :q975)
```

### Model Comparison

```julia
# Compare different model specifications
models = [model, refined_model, caliper_model]
model_names = ["Base", "Refined", "Caliper"]

println("Model Comparison (Period 5 ATT):")
for (i, (mod, name)) in enumerate(zip(models, model_names))
    if nrow(mod.results) >= 3  # Ensure period 5 exists
        period_5_att = mod.results[3, :att]  # Period 5 is 3rd row (3:8 range)
        println("$name: $(round(period_5_att, digits=4))")
    end
end
```

## Step 7: Diagnostic Checks

### Treatment Pattern Validation

```julia
# Verify treatment is event-based (not continuous)
treatment_check = combine(groupby(dat, :fips), 
    :gub => sum => :total_treatments)

continuous_treatment = filter(:total_treatments => >(1), treatment_check)
if nrow(continuous_treatment) > 0
    @warn "Found units with continuous treatment - this may cause issues"
end
```

### Time Period Coverage

```julia
# Check sufficient data coverage
time_range = extrema(dat.day)
println("Data covers days $(time_range[1]) to $(time_range[2])")

# Verify L and F periods are feasible
L_range = model.L
F_range = model.F
println("L periods: $L_range (need data $(L_range[1]) periods before treatment)")  
println("F periods: $F_range (need data $(F_range[end]) periods after treatment)")
```

### Sample Size Check

```julia
# Check adequate sample sizes
println("Sample size diagnostics:")
println("- Treated units: $(model.treatednum)")
println("- Control pool: $(length(model.ids))")
println("- Treated/Control ratio: $(round(model.treatednum/length(model.ids), digits=2))")

if model.treatednum < 10
    @warn "Small number of treated units may lead to imprecise estimates"
end
```

## Complete Analysis Example

Putting it all together:

```julia
# Complete workflow function
function tscs_analysis(data; title="Analysis")
    
    # 1. Model specification
    model = makemodel(data, :day, :fips, :gub, :death_rte,
                     [:pop_dens, :cumul_death_rate], 
                     Dict(:pop_dens => false, :cumul_death_rate => true),
                     3:8, -15:-5; title=title)
    
    # 2. Matching and balancing  
    match!(model, data)
    balance!(model, data)
    
    # 3. Estimation
    estimate!(model, data; dobayesfactor=false)
    
    # 4. Summary
    println("\\n=== $title Results ===")
    println("Treated units: $(model.treatednum)")
    println("Time periods: $(model.F)")
    println("\\nTreatment Effects:")
    
    for row in eachrow(model.results[1:3, :])  # First 3 periods
        period = row.f
        att = round(row.att, digits=4)
        ci = "[$(round(row.q025, digits=3)), $(round(row.q975, digits=3))]"
        println("  Period $period: $att $ci")
    end
    
    return model
end

# Run complete analysis
final_model = tscs_analysis(dat; title="Elections and COVID")
```

## Next Steps

- **Explore**: Try different L and F period specifications
- **Validate**: Check robustness with different matching strategies  
- **Extend**: Use multiple outcomes or stratification
- **Document**: Save results and create visualizations

This tutorial covered the complete TSCSMethods.jl workflow for causal inference with time-series cross-sectional data!