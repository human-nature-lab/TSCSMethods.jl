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
    _validate_imputation_inputs(m, matches, dat, tvar, unit_var, stratum)

Validate inputs for imputation function. Throws ArgumentError if validation fails.
"""
function _validate_imputation_inputs(m, matches, dat, tvar, unit_var, stratum)
    # Input validation
    isnothing(m) && throw(ArgumentError("Model cannot be nothing"))
    isnothing(matches) && throw(ArgumentError("Matches DataFrame cannot be nothing"))
    isnothing(dat) && throw(ArgumentError("Data DataFrame cannot be nothing"))
    
    # Check required columns exist
    tvar ∉ names(dat) && throw(ArgumentError("Time variable '$tvar' not found in data"))
    unit_var ∉ names(dat) && throw(ArgumentError("Unit variable '$unit_var' not found in data"))
    
    # Check model has required fields
    !hasproperty(m, :results) && throw(ArgumentError("Model must have results field"))
    !hasproperty(m, :outcome) && throw(ArgumentError("Model must have outcome field"))
    !hasproperty(m, :observations) && throw(ArgumentError("Model must have observations field"))
    
    # Check outcome variable exists in data
    m.outcome ∉ names(dat) && throw(ArgumentError("Outcome variable '$(m.outcome)' not found in data"))
    
    # Check matches has required columns
    required_match_cols = [:timetreated, :treatedunit, :matchunits, :f]
    for col in required_match_cols
        col ∉ names(matches) && throw(ArgumentError("Required column '$col' not found in matches"))
    end
    
    # Check stratum is valid for stratified models
    if m.stratifier != Symbol("") 
        stratum < 1 && throw(ArgumentError("Stratum must be >= 1"))
        !hasproperty(m, :strata) && throw(ArgumentError("Stratified model must have strata field"))
    end
    
    return nothing
end

"""
    _calculate_counterfactuals!(matches_with_treatment, outcome_lookup, outcome_var, tvar)

Calculate observed outcomes and counterfactual values for each matched observation.
Modifies matches_with_treatment in place.
"""
function _calculate_counterfactuals!(matches_with_treatment, outcome_lookup, outcome_var, tvar)
    for row in eachrow(matches_with_treatment)
        time_period = row[tvar]
        treated_unit = row[:treatedunit]
        matched_units = row[:matchunits]
        
        # Check for required data
        if !haskey(outcome_lookup, (time_period, treated_unit))
            throw(ArgumentError("Missing outcome data for treated unit $treated_unit at time $time_period"))
        end
        
        row[outcome_var] = outcome_lookup[(time_period, treated_unit)]
        
        # Calculate counterfactual from matched units
        matched_outcomes = Float64[]
        for unit in matched_units
            if haskey(outcome_lookup, (time_period, unit))
                push!(matched_outcomes, outcome_lookup[(time_period, unit)])
            else
                @warn "Missing outcome data for matched unit $unit at time $time_period, skipping"
            end
        end
        
        if isempty(matched_outcomes)
            throw(ArgumentError("No valid matched outcomes found for time $time_period"))
        end
        
        row[:counterfactual_values] = mean(matched_outcomes)
    end
    return nothing
end

"""
    _calculate_pretreatment_averages(matches_with_treatment, outcome_lookup, outcome_var)

Calculate pre-treatment outcome averages for treated and matched units.
Returns (matched_avg, treated_avg).
"""
function _calculate_pretreatment_averages(matches_with_treatment, outcome_lookup, outcome_var)
    # Extract unique treated and matched units for pre-treatment baseline calculation
    treated_obs = unique(matches_with_treatment[!, [:timetreated, :treatedunit]])
    matched_obs = unique(flatten(matches_with_treatment, :matchunits)[!, [:timetreated, :matchunits]])
    
    # Calculate pre-treatment time period (one period before treatment)
    # This is used to assess baseline balance between treated and matched units
    matched_obs[!, :pretreatment_time] = matched_obs.timetreated .- 1
    treated_obs[!, :pretreatment_time] = treated_obs.timetreated .- 1

    matched_obs[!, outcome_var] .= 0.0
    treated_obs[!, outcome_var] .= 0.0

    # Get pre-treatment outcomes for matched units with error handling
    for unit_obs in eachrow(matched_obs)
        key = (unit_obs[:pretreatment_time], unit_obs[:matchunits])
        if haskey(outcome_lookup, key)
            unit_obs[outcome_var] = outcome_lookup[key]
        else
            @warn "Missing pre-treatment data for matched unit $(unit_obs[:matchunits]) at time $(unit_obs[:pretreatment_time])"
            unit_obs[outcome_var] = missing
        end
    end

    # Get pre-treatment outcomes for treated units with error handling
    for unit_obs in eachrow(treated_obs)
        key = (unit_obs[:pretreatment_time], unit_obs[:treatedunit])
        if haskey(outcome_lookup, key)
            unit_obs[outcome_var] = outcome_lookup[key]
        else
            @warn "Missing pre-treatment data for treated unit $(unit_obs[:treatedunit]) at time $(unit_obs[:pretreatment_time])"
            unit_obs[outcome_var] = missing
        end
    end

    # Calculate pre-treatment averages with missing data handling
    matched_pretreatment_avg = mean(skipmissing(matched_obs[!, outcome_var]))
    treated_pretreatment_avg = mean(skipmissing(treated_obs[!, outcome_var]))
    
    # Validate final results
    if isnan(matched_pretreatment_avg)
        throw(ArgumentError("Unable to calculate matched pre-treatment average - insufficient data"))
    end
    if isnan(treated_pretreatment_avg)
        throw(ArgumentError("Unable to calculate treated pre-treatment average - insufficient data"))
    end
    
    return (matched_pretreatment_avg, treated_pretreatment_avg)
end

"""
    make_timeunit_lookup(dat, variable, time_col, unit_var)

Create lookup dictionary mapping (time, unit) tuples to variable values.
Returns Dict{Tuple, Union{Float64, Missing}} for fast value retrieval.

# Arguments
- `dat`: DataFrame with panel data
- `variable`: Column name for the variable to look up
- `time_col`: Column name for time variable
- `unit_col`: Column name for unit ID variable
"""
function make_timeunit_lookup(dat, variable, time_col, unit_col)
    odict = Dict{Tuple, Union{Float64, Missing}}()
    for r in eachrow(dat)
        odict[(r[time_col], r[unit_col])] = r[variable]
    end
    return odict
end

"""
    impute_results(m, matches, dat, tvar, unit_var; stratum = 1)

Generate counterfactual outcomes Y(0) using matched control units.

# Arguments
- `m`: Fitted TSCSMethods model
- `matches`: DataFrame of matched units from matching procedure  
- `dat`: Original panel data
- `tvar`: Time variable name in data
- `unit_var`: Unit ID variable name in data
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
function impute_results(m, matches, dat, tvar, unit_var; stratum = 1)
    
    # Validate all inputs
    _validate_imputation_inputs(m, matches, dat, tvar, unit_var, stratum)

    model_results = m.results
    outcome_var = m.outcome
    
    if m.stratifier != Symbol("")
        @subset!(model_results, :stratum .== 1)
    end

    # Create outcome lookup with error handling
    try
        outcome_lookup = make_timeunit_lookup(dat, m.outcome, tvar, unit_var)
    catch e
        throw(ArgumentError("Failed to create outcome lookup: $(e.msg)"))
    end
    
    # Check for missing critical data
    nrow(model_results) == 0 && throw(ArgumentError("Model results are empty"))
    nrow(matches) == 0 && throw(ArgumentError("Matches DataFrame is empty"))

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
    
    # Calculate outcomes and counterfactuals
    _calculate_counterfactuals!(matches_with_treatment, outcome_lookup, outcome_var, tvar)

    # Calculate pre-treatment baseline averages
    matched_pretreatment_avg, treated_pretreatment_avg = _calculate_pretreatment_averages(
        matches_with_treatment, outcome_lookup, outcome_var
    )

    time_aggregated = @chain matches_with_treatment begin
        groupby(:f)
        combine(
            outcome_var => mean => outcome_var,
            :counterfactual_values => mean => :counterfactual_values
        )
    end

    augmented_results = leftjoin(model_results, time_aggregated, on = :f)
    # Calculate counterfactual trajectory: what treated outcomes would have been without treatment
    # counterfactual = observed - treatment_effect
    augmented_results[!, :counterfactual_trajectory] = augmented_results[!, outcome_var] .- augmented_results.att

    # Calculate confidence intervals for counterfactual trajectory
    # Since counterfactual = observed - ATT, confidence bounds are:
    # counterfactual_bounds = observed - ATT_bounds
    # Use element-wise min/max to ensure lower < upper for each time period
    counterfactual_bound1 = augmented_results[!, outcome_var] .- augmented_results[!, Symbol("2.5%")]
    counterfactual_bound2 = augmented_results[!, outcome_var] .- augmented_results[!, Symbol("97.5%")]
    
    augmented_results[!, :counterfactual_lower] = min.(counterfactual_bound1, counterfactual_bound2)
    augmented_results[!, :counterfactual_upper] = max.(counterfactual_bound1, counterfactual_bound2)

    # Calculate percentage changes
    augmented_results[!, :pct_change] .= 0.0
    augmented_results[!, :pct_change_lower] .= 0.0
    augmented_results[!, :pct_change_upper] .= 0.0

    for row in eachrow(augmented_results)
        row[:pct_change] = change_pct(row[outcome_var], row[:att])
        row[:pct_change_lower] = change_pct(row[outcome_var], row[Symbol("2.5%")])
        row[:pct_change_upper] = change_pct(row[outcome_var], row[Symbol("97.5%")])
    end

    # Calculate day-to-day percentage change in counterfactual outcomes
    # Formula: 100 * (value[t] - value[t-1]) / value[t-1]
    # This helps assess stability of the counterfactual trajectory
    daily_variation = 100 .* diff(augmented_results.counterfactual_values) .* 
                     inv.(augmented_results.counterfactual_values[2:end])
    
    # Initialize column with missing values, then fill from second row onward
    augmented_results[!, :daily_counterfactual_change] = Vector{Union{Missing, Float64}}(undef, nrow(augmented_results))
    augmented_results[2:end, :daily_counterfactual_change] = daily_variation

    return ImputationResults(augmented_results, matched_pretreatment_avg, treated_pretreatment_avg)
end
