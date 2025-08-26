# mean_fullbalancing.jl

"""
    outcome_to_lag_indices(outcome_period::Int, matching_window_length::Int, min_outcome_period::Int)

Convert outcome period to corresponding lag period indices for time-varying covariates.

# Purpose
Maps an outcome period (F) to the appropriate range of lag period indices (L) 
for extracting time-varying covariate values during balance calculations.

# Mathematical Mapping
For outcome period `outcome_period`, returns the range:
```
(outcome_period - min_outcome_period + 1):(outcome_period - min_outcome_period + matching_window_length)
```

# Arguments
- `outcome_period`: Outcome period index (from F range)
- `matching_window_length`: Length of matching window (number of L periods)
- `min_outcome_period`: Minimum outcome period value

# Returns
Range of indices into the lag period structure for covariate extraction

# Example
```julia
# If F = 51:60, L = -10:-1 (so matching_window_length = 10, min_outcome_period = 51)
# For outcome period 53:
outcome_to_lag_indices(53, 10, 51) # returns 3:12
# This maps to lag periods corresponding to outcome period 53
```
"""
outcome_to_lag_indices(outcome_period::Int, matching_window_length::Int, min_outcome_period::Int) = 
    (outcome_period - min_outcome_period + 1):(outcome_period - min_outcome_period + matching_window_length);

"""
    mean_fullbalance!(model)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function mean_fullbalance!(model::AbstractCICModel, bals::DataFrame)

    (; matches, ids, covariates, timevary, F, L) = refcalmodel

    fmin = minimum(F); mmlen = length(L);

    initialize_balance_storage!(refcalmodel);

    _meanbalance!(
        eachrow(refcalmodel.meanbalances), eachrow(bals),
        matches, ids, covariates, timevary, mmlen, fmin
    );

    return model
end

# meanbalances loop
function _meanbalance!(
    barrows, balrows,
    matches, ids, covariates, timevary, mmlen, fmin
)

    Threads.@threads :greedy for i in eachindex(balrows)
        br = balrows[i];
        mr = barrows[i];

        _, efsets = matchassignments(matches[i], ids);

        _covar_meanbalance!(
            mr, br, efsets, covariates, timevary, mmlen, fmin,
        );
    end
    return barrows
end

# inner functions
function _covar_meanbalance!(
    mr, br, efsets, covariates, timevary, mmlen, fmin,
)
    for covar in covariates
        # this is a [row, covar] for MB:
        #   the means across fs for a covariate (for a treated observation)
        if timevary[covar]
            row_covar_meanbalance!(
                mr[covar], br[covar], mr[:fs], efsets, mmlen, fmin
            );
        else
            row_covar_meanbalance!(mr[covar], br[covar], mr[:fs], efsets);
        end
    end
end

"""
    row_covar_meanbalance!(Holding, br_covar, fs, efsets, mmlen, fmin)

Calculate mean balance for a single covariate across outcome periods with time-varying support.

# Purpose
This function computes balance statistics by averaging covariate values across matched units
for each valid outcome period. For time-varying covariates, it accounts for temporal
structure in the matching process.

# Algorithm
1. **Outcome Period Iteration**: For each F period (outcome estimation window)
2. **Match Collection**: Gather all eligible matches for that outcome period  
3. **Temporal Indexing**: Map outcome period to appropriate lag period indices
4. **Balance Calculation**: Average covariate values across valid matches

# Mathematical Framework
For outcome period f and covariate c:
```
Balance[f][t] = (1/|M_f|) ∑_{m ∈ M_f} covariate_value[m][lag_time[t]]
```
where:
- M_f = set of matches eligible for outcome period f
- lag_time[t] = temporal lag corresponding to index t

# Arguments
- `Holding`: Output array storing mean balance values per outcome period
- `br_covar`: Balance results for the covariate across all matches
- `fs`: Boolean vector indicating valid outcome periods  
- `efsets`: Eligibility sets [match][outcome_period] indicating match validity
- `mmlen`: Length of matching window (L periods)
- `fmin`: Minimum outcome period value

# Implementation Notes
- Handles missing values through `skipmissing()` operations
- Pre-allocates temporary storage for efficient computation
- Uses `fidx()` to map outcome periods to lag period indices
"""
function row_covar_meanbalance!(
    Holding::Vector{Vector{Union{Missing, Float64}}},
    br_covar,
    fs, efsets, mmlen, fmin
)

    for (outcome_period_index, f) in enumerate(fs)
        if f # if there are any valid matches for this outcome period
            # Map outcome period to corresponding lag period indices for time-varying covariates
            ls = outcome_to_lag_indices(outcome_period_index + fmin - 1, mmlen, fmin);
            
            # Determine which matches are eligible for this specific outcome period
            matches_eligible_for_period = [efsets[m][outcome_period_index] for m in 1:length(br_covar)]
            
            # Pre-allocate storage for covariate values from eligible matches
            holding = Vector{Vector{Union{Float64, Missing}}}(
                undef, sum(matches_eligible_for_period)
            );
            
            # Collect covariate values from all eligible matches for this outcome period
            match_count = 0
            for m in eachindex(br_covar) 
                if efsets[m][outcome_period_index]  # if this match is valid for this outcome period
                    match_count += 1
                    # Extract covariate values over the lag periods for this match
                    holding[match_count] = br_covar[m][ls];
                end
            end
        end

        # Calculate mean balance across eligible matches for this outcome period
        __row_covar_meanbalance!(Holding, holding, outcome_period_index, ls);
    end
    return Holding
end

function __row_covar_meanbalance!(Holding, holding, outcome_period_index, ls)
    for l in eachindex(ls)
        lvec = Vector{Union{Float64, Missing}}(missing, length(holding));
        for (v, hold) in enumerate(holding)
            lvec[v] = hold[l]
        end
        
        Holding[outcome_period_index][l] = !isempty(skipmissing(lvec)) ? mean(skipmissing(lvec)) : missing
    end
    return Holding
end

"""
    row_covar_meanbalance!(Holding, br_covar, fs, efsets)

Calculate mean balance for a single covariate across outcome periods (time-invariant version).

# Purpose
Simplified version for time-invariant covariates that computes balance statistics
by averaging covariate values across matched units for each valid outcome period.

# Algorithm
1. **Outcome Period Iteration**: For each F period (outcome estimation window)
2. **Match Collection**: Gather all eligible matches for that outcome period
3. **Simple Averaging**: Calculate mean across valid matches (no temporal structure)

# Mathematical Framework
For outcome period f and time-invariant covariate c:
```
Balance[f] = (1/|M_f|) ∑_{m ∈ M_f} covariate_value[m]
```
where M_f = set of matches eligible for outcome period f

# Arguments
- `Holding`: Output vector storing mean balance values per outcome period
- `br_covar`: Balance results for the covariate across all matches
- `fs`: Boolean vector indicating valid outcome periods
- `efsets`: Eligibility sets [match][outcome_period] indicating match validity

# Performance Notes
- More efficient than time-varying version due to simpler indexing
- Uses `skipmissing()` for robust handling of missing values
"""
function row_covar_meanbalance!(Holding, br_covar, fs, efsets)
    for (outcome_period_index, f) in enumerate(fs)
        if f  # if there are any valid matches for this outcome period
            # Pre-allocate storage for covariate values from all potential matches
            holding = Vector{Union{Float64, Missing}}(undef, length(br_covar));
            
            # Collect covariate values from matches eligible for this outcome period
            for m in eachindex(br_covar)
                if efsets[m][outcome_period_index]  # if this match is valid for this outcome period
                    holding[m] = br_covar[m];  # time-invariant: single value per match
                end
            end
            
            # Calculate mean balance across all eligible matches, handling missing values
            Holding[outcome_period_index] = mean(skipmissing(holding))
        end
    end
    return Holding
end

function refinebalances(refcalmodel, model, balances)
    (; covariates, matches, observations) = model;
    bals = copy(balances);

    # remove tobs that do not exist in ref / cal model
    if length(observations) != length(refcalmodel.observations)
        keep = observations .∈ Ref(refcalmodel.observations);
        bals = bals[keep, :];
    end

    for covar in covariates
        for i in eachindex(refcalmodel.observations)
            rdx = refcalmodel.matches[i].eligible_matches[matches[i].eligible_matches] # 2096
            bals[i, covar] = balances[i, covar][rdx]
        end
    end
    return bals
end

"""
    balances!(refcalmodel, model, balances)

Complete balance calculation pipeline for refined/calipered models.

# Purpose
This is the main entry point for computing balance statistics when working with 
refined or calipered models. It orchestrates the full balance calculation workflow,
transforming raw balance data into final balance metrics suitable for assessment.

# Algorithm
1. **Balance Refinement**: Adjust balance data to match the refined/calipered model structure
2. **Mean Balance Calculation**: Compute mean balances across matches for each treated unit
3. **Grand Balance Computation**: Calculate overall balance statistics across all treated units

# Workflow Integration
This function is typically called after:
- Initial matching (`match!`)
- Model refinement (`refine!`) or caliper application (`caliper!`)
- Balance data generation

# Mathematical Framework
The function computes a hierarchy of balance statistics:
```
Raw Balances → Refined Balances → Mean Balances → Grand Balances
```

Where:
- **Refined Balances**: Subset to matches present in refined model
- **Mean Balances**: Average across eligible matches per treated unit
- **Grand Balances**: Population-level balance assessment

# Arguments
- `refcalmodel`: Refined or calipered model (output structure)
- `model`: Original full model (reference structure)  
- `balances`: Raw balance data from initial matching

# Returns
- `refcalmodel`: Updated model with complete balance statistics

# Usage Example
```julia
# After matching and refinement
refined_model = refine!(model, refinement_rules)
balances!(refined_model, original_model, initial_balances)
# refined_model now contains complete balance statistics
```

# Performance Notes
- Operates in-place on `refcalmodel` for memory efficiency
- Threading used internally for large datasets
- Balance refinement step filters data to avoid unnecessary computation
"""
function balances!(refcalmodel, model, balances)
    # Step 1: Refine balance data to match the structure of the refined/calipered model
    refined_balance_data = refinebalances(refcalmodel, model, balances);
    
    # Step 2: Calculate mean balances for each treated observation
    mean_fullbalance!(refcalmodel, refined_balance_data);
    
    # Step 3: Compute grand balance statistics across all treated units
    grandbalance!(refcalmodel);
    
    return refcalmodel
end
