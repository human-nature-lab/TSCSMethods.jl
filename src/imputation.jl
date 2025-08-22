# imputation.jl

"""
    ImputationResults

Container for counterfactual imputation results.

# Fields
- `results`: DataFrame with counterfactual trajectories and treatment effects
- `matched_pretreatment_avg`: Average pre-treatment outcome for matched controls
- `treated_pretreatment_avg`: Average pre-treatment outcome for treated units
- `baseline_difference`: Difference in pre-treatment averages (treated - matched)
"""
struct ImputationResults
    results::DataFrame
    matched_pretreatment_avg::Float64
    treated_pretreatment_avg::Float64
    baseline_difference::Float64
    
    function ImputationResults(results, matched_avg, treated_avg)
        new(results, matched_avg, treated_avg, treated_avg - matched_avg)
    end
end

function Base.show(io::IO, ir::ImputationResults)
    println(io, "ImputationResults:")
    println(io, "  Results: $(nrow(ir.results)) time periods")
    println(io, "  Pre-treatment averages:")
    println(io, "    Treated units: $(round(ir.treated_pretreatment_avg, digits=3))")
    println(io, "    Matched controls: $(round(ir.matched_pretreatment_avg, digits=3))")
    print(io, "    Baseline difference: $(round(ir.baseline_difference, digits=3))")
end

"""
Explanation of imputation.jl

This file implements counterfactual outcome imputation for causal
inference analysis. Here's what it does:

Core Purpose

The impute_results function creates counterfactual outcomes Y(0) - what
would have happened to treated units if they had not been treated - by
using the matched control units.

Key Functions

make_timeunit_lookup(dat, variable) (lines 3-9):
- Creates a lookup dictionary mapping (time, unit_id) tuples to outcome
values
- Problem: Hardcoded column names r.running and r.fips - should be
parameterized

impute_results(m, matches, dat, tvar; stratum = 1) (lines 14-107):
Main imputation function that:

1. Data Preparation (lines 16-44):
- Extracts model results and outcome variable
- Handles stratified vs non-stratified models
- Creates lookup dictionary for outcome values
- Joins matches with treatment information
2. Counterfactual Calculation (lines 45-53):
- For each treated observation, calculates:
    - Actual outcome: odict[(τ, treatedunit)]
    - Counterfactual: mean([odict[(τ, mu)] for mu in mus]) (average of
matched controls)
3. Pre-treatment Baseline (lines 55-70):
- Calculates pre-treatment outcome averages for both treated and
matched units
- Used for assessing balance and baseline equivalence
4. Results Augmentation (lines 71-87):
- Combines ATT estimates with counterfactual predictions
- Creates confidence intervals for counterfactual outcomes
- counter_post = observed - att (the counterfactual trajectory)
5. Percentage Change Calculation (lines 91-99):
- Converts ATT estimates to percentage changes
- Useful for interpretation (e.g., "treatment increased outcomes by
X%")
6. Day-to-day Variability (lines 101-105):
- Calculates percent change in matched control outcomes between
consecutive days
- Helps assess stability of the counterfactual
"""

"""
    make_timeunit_lookup(dat, variable)

Create lookup dictionary mapping (time, unit) tuples to variable values.
Returns Dict{Tuple, Union{Float64, Missing}} for fast value retrieval.
"""
function make_timeunit_lookup(dat, variable)
    odict = Dict{Tuple, Union{Float64, Missing}}();
    for r in eachrow(dat)
        odict[(r.running, r.fips)] = r[variable]
    end
    return odict
end

"""
    impute_results(m, matches, dat, tvar; stratum = 1)

Generate counterfactual outcomes Y(0) using matched control units.

# Arguments
- `m`: Fitted TSCSMethods model
- `matches`: DataFrame of matched units from matching procedure  
- `dat`: Original panel data
- `tvar`: Time variable name
- `stratum`: Stratum number for stratified models (default: 1)

# Returns
- `ImputationResults`: Structured container with counterfactual analysis

The results DataFrame includes:
- `counterfactual_trajectory`: What treated outcomes would have been without treatment
- `counterfactual_values`: Average of matched control outcomes at each time
- `counterfactual_lower`/`counterfactual_upper`: Confidence intervals for counterfactuals
- `pct_change`/`pct_change_lower`/`pct_change_upper`: Treatment effects as percentage changes
- `daily_counterfactual_change`: Day-to-day variation in counterfactual outcomes

Creates counterfactual predictions by averaging matched control outcomes
at each time period, enabling visualization of treated vs counterfactual trajectories.
"""
function impute_results(m, matches, dat, tvar; stratum = 1)

    model_results = m.results
    outcome_var = m.outcome
    
    if m.stratifier != Symbol("")
        @subset!(model_results, :stratum .== 1)
    end

    outcome_lookup = make_timeunit_lookup(dat, m.outcome)

    treatment_info = if m.stratifier != Symbol("")
        DataFrame(
            :timetreated => [obs[1] for obs in m.observations],
            :treatedunit => [obs[2] for obs in m.observations],
            :stratum => m.strata
        )
    else
        DataFrame(
            :timetreated => [obs[1] for obs in m.observations],
            :treatedunit => [obs[2] for obs in m.observations],
        )
    end
    
    matches_with_treatment = leftjoin(matches, treatment_info, on = [:timetreated, :treatedunit])

    if m.stratifier != Symbol("")
        # Filter to requested stratum
        @subset!(matches_with_treatment, :stratum .== stratum)
    end

    matches_with_treatment = @transform(matches_with_treatment, $(tvar) = :timetreated + :f)
    matches_with_treatment[!, outcome_var] .= 0.0
    matches_with_treatment[!, :counterfactual_values] .= 0.0
    
    for row in eachrow(matches_with_treatment)
        time_period = row[:running]
        row[outcome_var] = outcome_lookup[(time_period, row[:treatedunit])]
        matched_units = row[:matchunits]
        row[:counterfactual_values] = mean([outcome_lookup[(time_period, unit)] for unit in matched_units])
    end

    treated_obs = unique(matches_with_treatment[!, [:timetreated, :treatedunit]])
    matched_obs = unique(flatten(matches_with_treatment, :matchunits)[!, [:timetreated, :matchunits]])
    matched_obs[!, :pretreatment_time] = matched_obs.timetreated .- 1
    treated_obs[!, :pretreatment_time] = treated_obs.timetreated .- 1

    matched_obs[!, outcome_var] .= 0.0
    treated_obs[!, outcome_var] .= 0.0

    for unit_obs in eachrow(matched_obs)
        unit_obs[outcome_var] = outcome_lookup[(unit_obs[:pretreatment_time], unit_obs[:matchunits])]
    end

    for unit_obs in eachrow(treated_obs)
        unit_obs[outcome_var] = outcome_lookup[(unit_obs[:pretreatment_time], unit_obs[:treatedunit])]
    end

    time_aggregated = @chain matches_with_treatment begin
        groupby(:f)
        combine(
            outcome_var => mean => outcome_var,
            :counterfactual_values => mean => :counterfactual_values
        )
    end

    augmented_results = leftjoin(model_results, time_aggregated, on = :f)
    augmented_results[!, :counterfactual_trajectory] = augmented_results[!, outcome_var] .- augmented_results.att

    # Calculate confidence intervals for counterfactual trajectory
    conf_lower, conf_upper = sort([
        augmented_results[!, outcome_var] .- augmented_results[!, Symbol("2.5%")], 
        augmented_results[!, outcome_var] .- augmented_results[!, Symbol("97.5%")]
    ])

    augmented_results[!, :counterfactual_lower] = conf_lower
    augmented_results[!, :counterfactual_upper] = conf_upper

    matched_pretreatment_avg = mean(matched_obs[!, m.outcome])
    treated_pretreatment_avg = mean(treated_obs[!, m.outcome])

    # Calculate percentage changes
    augmented_results[!, :pct_change] .= 0.0
    augmented_results[!, :pct_change_lower] .= 0.0
    augmented_results[!, :pct_change_upper] .= 0.0

    for row in eachrow(augmented_results)
        row[:pct_change] = change_pct(row[outcome_var], row[:att])
        row[:pct_change_lower] = change_pct(row[outcome_var], row[Symbol("2.5%")])
        row[:pct_change_upper] = change_pct(row[outcome_var], row[Symbol("97.5%")])
    end

    # Day-to-day variability in counterfactual
    daily_variation = 100 .* diff(augmented_results.counterfactual_values) .* 
                     inv.(augmented_results.counterfactual_values[2:end])
    augmented_results[!, :daily_counterfactual_change] = Vector{Union{Missing, Float64}}(undef, nrow(augmented_results))
    augmented_results[2:end, :daily_counterfactual_change] = daily_variation

    return ImputationResults(augmented_results, matched_pretreatment_avg, treated_pretreatment_avg)
end
