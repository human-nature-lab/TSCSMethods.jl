# information.jl
# various functions that give information about the models, etc.

"""
    unitcounts(m)

count the total number of treated and the number of matched units for each F.
"""
function unitcounts(m)
    
    X = Vector{Vector{Int}}(undef, 0);

    Flen = length(m.F)
    
    for i in eachindex(m.matches)
        push!(
            X,
            vec(sum(m.matches[i].eligible_matches, dims = 1))
        )
    end

    strat = any(
        [
            typeof(m) == x for x in [
                CICStratified, RefinedCICStratified, RefinedCaliperCICStratified
            ]
        ]
    );

    if !strat
        # compute the number of treated units left (for each F)
        Y = zeros(Int, Flen)
        _treatedcount!(Y, X)
        
        # compute the total number of match units (for each F)
        U = sum(X)
        return Y, U
    else
        S = sort(unique(m.strata));
        # with stratification
        Ys = Dict{Int, Vector{Int}}();
        Us = Dict{Int, Vector{Int}}();

        # compute the number of treated units left (for each F)
        for (i, s) in enumerate(S)
            Ys[s] = zeros(Int, Flen)
            Us[s] = zeros(Int, Flen)
            c1 = m.strata .== s;
            _treatedcount!(Ys[s], X[c1])
        
            # compute the total number of match units (for each F)
            Us[i] = sum(X[c1])
        end
        return Ys, Us
    end
end

function _treatedcount!(Y, X)
    for x in X
        for (i, xi) in enumerate(x)
            if xi > 0
                Y[i] += 1
            end
        end
    end
end

"""
    eligibility(model::VeryAbstractCICModel)

Calculate unit eligibility across all treated observations.

# Arguments  
- `model`: Fitted TSCSMethods model

# Returns
- `Matrix{Int}`: Eligibility matrix (units × time periods F)

# Description
For each potential control unit, calculates how many times it is eligible
to serve as a match across all treated observations and time periods.
Higher eligibility indicates units that are consistently good potential matches
for multiple treated units.

Useful for identifying:
- Units that are consistently good matches across multiple treated observations
- Potential issues with match quality or overlap
- Whether certain units dominate the matching process
- Balance in the matching pool

# Examples
```julia
eligibility_matrix = eligibility(model)

# Find highly eligible units
total_eligibility = sum(eligibility_matrix, dims=2)
highly_eligible_indices = findall(x -> x > 5, vec(total_eligibility))
highly_eligible_ids = model.ids[highly_eligible_indices]

println("Units eligible for >5 matches: \$highly_eligible_ids")

# Check eligibility by time period
for (f_idx, f) in enumerate(model.F)
    eligible_count = sum(eligibility_matrix[:, f_idx])
    println("F=\$f: \$eligible_count total eligible matches")
end
```
"""
function eligibility(model::VeryAbstractCICModel)
    # Input validation
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    
    # Initialize eligibility matrix with same structure as match eligibility
    total_eligibility = zeros(Int, size(model.matches[1].eligible_matches))
    
    # Sum eligibility across all treated observations
    for match in model.matches
        total_eligibility += match.eligible_matches
    end
    
    return total_eligibility
end

"""
    whentreated(unit_id::Int, model::VeryAbstractCICModel)

Find when a specific unit received treatment.

# Arguments
- `unit_id`: ID of the unit to look up  
- `model`: Fitted TSCSMethods model

# Returns
- `Vector{Tuple{Int, Int}}`: Treatment events for this unit as (treatment_time, unit_id)
- `"Unit not found"`: If unit was never treated

# Description
Searches through all treatment observations to find when (if ever)
a specific unit received treatment. Returns all treatment events
for units that may have multiple treatments over time.

This function is useful for:
- Verifying treatment timing for specific units
- Checking if a unit received treatment multiple times
- Debugging treatment assignment issues
- Understanding the treatment history of particular units

# Examples
```julia
# Check when unit 1001 was treated
treatment_events = whentreated(1001, model)

if treatment_events isa String
    println("Unit 1001 was never treated")
else
    treatment_times = [event[1] for event in treatment_events]
    println("Unit 1001 treated at times: \$treatment_times")
    
    if length(treatment_events) > 1
        println("Note: Unit had multiple treatments!")
    end
end

# Check multiple units
for unit_id in [1001, 1002, 1003]
    events = whentreated(unit_id, model)
    if events isa String
        println("Unit \$unit_id: never treated")
    else
        println("Unit \$unit_id: treated \$(length(events)) times")
    end
end
```
"""
function whentreated(unit_id::Int, model::VeryAbstractCICModel)
    # Input validation
    if isempty(model.observations)
        throw(ArgumentError("Model has no treatment observations"))
    end
    
    # Find all treatment events for this unit using modern functional approach
    unit_events = filter(obs -> obs[2] == unit_id, model.observations)
    
    if isempty(unit_events)
        return "Unit not found"
    end
    
    return unit_events
end

"""
    showmatches(model::VeryAbstractCICModel, treatment_observation::Tuple{Int, Int})

Show ranked matches for a specific treatment observation.

# Arguments
- `model`: Fitted TSCSMethods model  
- `treatment_observation`: Tuple of (treatment_time, unit_id)

# Returns
- `Vector{Vector{Int}}`: Ranked matches for each time period F (sorted best to worst)
- `"Treatment observation not found"`: If observation doesn't exist in model

# Description
For a specific treatment event, shows the ranked list of matched control unit indices
for each time period in the post-treatment window F. Units are ranked from
best match (rank 1) to worst match based on the matching algorithm's distance
calculations and balancing constraints.

The returned structure is a vector where each element corresponds to a time period F,
containing the ranked list of control unit indices (internal to model.matches).
To get actual unit IDs, use `model.ids[index]`.

This function is useful for:
- Debugging matching quality for specific treatment events
- Understanding which units serve as matches across different time periods
- Assessing match consistency over time
- Manual inspection of matching results

# Examples
```julia
# Show matches for unit 1001 treated at time 15
matches = showmatches(model, (15, 1001))

if matches isa String
    println(matches)  # "Treatment observation not found"
else
    println("Matches across F periods:")
    for (f_idx, f) in enumerate(model.F)
        period_matches = matches[f_idx]
        if !isempty(period_matches)
            # Convert to actual unit IDs
            match_ids = [model.ids[idx] for idx in period_matches[1:min(3, length(period_matches))]]
            println("F=\$f: top 3 matches = \$match_ids")
        else
            println("F=\$f: no matches available")
        end
    end
end

# Check consistency of top match across periods
if matches isa Vector && all(length.(matches) .> 0)
    top_matches = [model.ids[period_matches[1]] for period_matches in matches]
    println("Top match by period: \$top_matches")
end
```
"""
function showmatches(model::VeryAbstractCICModel, treatment_observation::Tuple{Int, Int})
    # Input validation
    if isempty(model.observations)
        throw(ArgumentError("Model has no treatment observations"))
    end
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    
    # Find the index of this treatment observation using modern approach
    obs_index = findfirst(obs -> obs == treatment_observation, model.observations)
    
    if isnothing(obs_index)
        return "Treatment observation not found"
    end
    
    # Return sorted rankings for this observation across all F periods
    return [sort(ranking) for ranking in model.matches[obs_index].match_rankings]
end

"""
    matchinfo(model::Union{CIC, CICStratified}; maxrank = 5)

Generate match information DataFrame for standard CIC models.

# Purpose
Creates a detailed summary of matching results showing which control units
were matched to each treated observation across all forward periods (F).
Provides match unit IDs and their ranking within each time period.

# Arguments
- `model`: Fitted CIC or CICStratified model with completed matching
- `maxrank`: Maximum number of top-ranked matches to include per F period (default: 5)

# Returns
`DataFrame` with columns:
- `timetreated`: Time when treatment occurred
- `treatedunit`: ID of the treated unit
- `f`: Forward period within the outcome window
- `matchunits`: Vector of matched control unit IDs (ranked best to worst)
- `ranks`: Vector of rank positions (1 = best match, 2 = second best, etc.)

# Description
For each treatment observation and each forward period F, this function:
1. Extracts the top `maxrank` matched control units
2. Converts internal match indices to actual unit IDs
3. Filters out periods with no available matches
4. Returns results sorted by treatment time and unit ID

This is useful for:
- Inspecting match quality and consistency across time periods
- Understanding which units serve as controls for specific treatments
- Debugging matching algorithm performance
- Creating match summaries for reporting

# Examples
```julia
# Get match information with top 3 matches per period
matches_df = matchinfo(model; maxrank = 3)

# View matches for a specific treatment
unit_matches = filter(row -> row.treatedunit == 1001 && row.timetreated == 15, matches_df)
for row in eachrow(unit_matches)
    println("F=\$(row.f): matched to units \$(row.matchunits) with ranks \$(row.ranks)")
end

# Check if any treatment lacks matches
no_matches = filter(row -> isempty(row.matchunits), matches_df)
if !isempty(no_matches)
    println("Warning: Some treatments have no matches")
end
```
"""
function matchinfo(model::Union{CIC, CICStratified}; maxrank::Int = 5)
    # Input validation
    if maxrank <= 0
        throw(ArgumentError("maxrank must be positive, got $(maxrank)"))
    end
    if isempty(model.observations)
        throw(ArgumentError("Model has no treatment observations"))
    end
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    
    # Extract model components
    (; observations, matches, F, ids) = model
    
    # Pre-allocate result DataFrame with proper types
    match_results = DataFrame(
        timetreated = Int[],
        treatedunit = Int[],
        f = Int[],
        matchunits = Vector{Int}[],
        ranks = Vector{Int}[]
    )
    
    # Build match information for each treatment observation and F period
    for (obs_idx, treatment_obs) in enumerate(observations)
        treatment_time, treated_unit = treatment_obs
        
        for f in F
            # Calculate index for this F period
            period_index = f - minimum(F) + 1
            
            # Get ranked matches for this period (up to maxrank)
            period_matches = matches[obs_idx].match_rankings[period_index]
            match_count = min(length(period_matches), maxrank)
            
            # Skip periods with no matches
            if match_count == 0
                continue
            end
            
            # Extract top matches and convert indices to unit IDs
            top_match_indices = period_matches[1:match_count]
            match_unit_ids = [ids[idx] for idx in top_match_indices]
            match_ranks = collect(1:match_count)
            
            # Add row to results
            push!(match_results, (
                timetreated = treatment_time,
                treatedunit = treated_unit,
                f = f,
                matchunits = match_unit_ids,
                ranks = match_ranks
            ))
        end
    end
    
    # Sort results by treatment time and unit for consistent output
    sort!(match_results, [:timetreated, :treatedunit, :f])
    
    return match_results
end

"""
    matchinfo(refined_model, original_model; maxrank = 5)

Generate match information DataFrame for refined/calipered models.

# Purpose
Creates match information for models that have undergone refinement or
caliper operations, showing how refined matches relate to the original
full model's ranking system. Useful for understanding how refinement
affects match selection and quality.

# Arguments
- `refined_model`: Refined or calipered CIC model (result of `refine()` or `caliper()`)
- `original_model`: Original full CIC model before refinement
- `maxrank`: Maximum number of matches to include per F period (default: 5)

# Returns
`DataFrame` with columns:
- `timetreated`: Time when treatment occurred
- `treatedunit`: ID of the treated unit
- `f`: Forward period within the outcome window
- `matchunits`: Vector of matched control unit IDs from refined model
- `ranks`: Vector of rank positions from original model's full ranking

# Description
For each treatment observation in the refined model:
1. Finds the corresponding observation in the original model
2. Extracts eligible matches from the refined model
3. Determines their ranks in the original model's full ranking
4. Returns matches sorted by their original ranks

This allows comparison between refined and original matching, showing
whether refinement preserves the best matches or introduces changes.

# Examples
```julia
# Compare refined vs original matches
original_matches = matchinfo(original_model; maxrank = 10)
refined_matches = matchinfo(refined_model, original_model; maxrank = 10)

# Check if refinement preserved top matches
for (orig_row, ref_row) in zip(eachrow(original_matches), eachrow(refined_matches))
    if orig_row.treatedunit == ref_row.treatedunit && orig_row.f == ref_row.f
        original_top3 = orig_row.matchunits[1:min(3, length(orig_row.matchunits))]
        refined_units = ref_row.matchunits
        preserved = intersect(original_top3, refined_units)
        println("Unit \$(orig_row.treatedunit), F=\$(orig_row.f): \$(length(preserved))/3 top matches preserved")
    end
end
```
"""
function matchinfo(refined_model, original_model; maxrank::Int = 5)
    # Input validation
    if maxrank <= 0
        throw(ArgumentError("maxrank must be positive, got $(maxrank)"))
    end
    if isempty(refined_model.observations)
        throw(ArgumentError("Refined model has no treatment observations"))
    end
    if isempty(original_model.observations)
        throw(ArgumentError("Original model has no treatment observations"))
    end
    if isempty(refined_model.matches) || isempty(original_model.matches)
        throw(ArgumentError("Both models must have completed matching"))
    end
    
    # Pre-allocate result DataFrame
    match_results = DataFrame(
        timetreated = Int[],
        treatedunit = Int[],
        f = Int[],
        matchunits = Vector{Int}[],
        ranks = Vector{Int}[]
    )
    
    # Process each treatment observation in refined model
    for (refined_idx, treatment_obs) in enumerate(refined_model.observations)
        treatment_time, treated_unit = treatment_obs
        
        # Find corresponding observation in original model
        original_idx = findfirst(obs -> obs == treatment_obs, original_model.observations)
        
        if isnothing(original_idx)
            throw(ArgumentError("Treatment observation $(treatment_obs) not found in original model"))
        end
        
        # Process each F period
        for f in refined_model.F
            period_index = f - minimum(refined_model.F) + 1
            
            # Get eligible matches from refined model (as unit indices)
            eligible_indices = refined_model.matches[refined_idx].eligible_matches[:, period_index]
            
            # Convert to unit IDs and limit to maxrank
            eligible_unit_ids = original_model.ids[eligible_indices]
            selected_unit_ids = eligible_unit_ids[1:min(length(eligible_unit_ids), maxrank)]
            
            # Skip if no matches available
            if isempty(selected_unit_ids)
                continue
            end
            
            # Get full ranking from original model to determine ranks
            original_period_index = f - minimum(original_model.F) + 1
            full_ranking_indices = original_model.matches[original_idx].match_rankings[original_period_index]
            full_ranking_unit_ids = original_model.ids[full_ranking_indices]
            
            # Find ranks of selected units in the original full ranking
            ranks_in_original = [findfirst(==(unit_id), full_ranking_unit_ids) for unit_id in selected_unit_ids]
            
            # Filter out units not found in original ranking (shouldn't happen but be safe)
            valid_mask = .!isnothing.(ranks_in_original)
            final_unit_ids = selected_unit_ids[valid_mask]
            final_ranks = ranks_in_original[valid_mask]
            
            # Sort by original ranking to maintain rank order
            sort_order = sortperm(final_ranks)
            sorted_unit_ids = final_unit_ids[sort_order]
            sorted_ranks = final_ranks[sort_order]
            
            # Add to results if we have valid matches
            if !isempty(sorted_unit_ids)
                push!(match_results, (
                    timetreated = treatment_time,
                    treatedunit = treated_unit,
                    f = f,
                    matchunits = sorted_unit_ids,
                    ranks = sorted_ranks
                ))
            end
        end
    end
    
    # Sort results by treatment time and unit for consistent output
    sort!(match_results, [:timetreated, :treatedunit, :f])
    
    return match_results
end

"""
    obsinfo(
        match_info_df, data, variables;
        full_model_observations = nothing,
        time_column = :t, id_column = :id
    )

Extract covariate information for treated units at their treatment times.

# Purpose
Gathers covariate values for all treated units at the exact time of their
treatment events. For time-varying covariates, this ensures values are
captured at the moment of treatment rather than at arbitrary time points.
Useful for creating treatment group summaries and balance tables.

# Arguments
- `match_info_df`: DataFrame from `matchinfo()` containing treatment observations
- `data`: Original time-series cross-sectional dataset
- `variables`: Vector of column names to extract (e.g., `[:population, :gdp, :policy_index]`)
- `full_model_observations`: Optional vector of all treatment observations for comparison (default: `nothing`)
- `time_column`: Name of time variable in data (default: `:t`)
- `id_column`: Name of unit identifier in data (default: `:id`)

# Returns
`DataFrame` with columns:
- `timetreated`: Time when treatment occurred
- `treatedunit`: ID of the treated unit
- `removed`: Boolean indicating if observation was removed in refinement (if `full_model_observations` provided)
- Additional columns for each variable in `variables`

# Description
For each unique treatment observation:
1. Identifies the corresponding row in the original data
2. Extracts covariate values at the time of treatment
3. Optionally compares against full model to identify removed observations
4. Returns sorted results with all requested covariate information

This is essential for:
- Creating balance tables showing pre-treatment characteristics
- Comparing refined vs full model treatment groups
- Understanding treatment group composition
- Generating summary statistics for reporting

# Examples
```julia
# Get basic covariate info for treated units
match_df = matchinfo(model)
covariates = [:population, :gdp_per_capita, :unemployment_rate]
treated_info = obsinfo(match_df, data, covariates; 
                      time_column = :year, id_column = :state_id)

# Compare full vs refined model
full_obs = model.observations
refined_df = matchinfo(refined_model)
comparison = obsinfo(refined_df, data, covariates;
                    full_model_observations = full_obs,
                    time_column = :year, id_column = :state_id)

# View removed observations
removed_units = filter(row -> row.removed, comparison)
println("Removed \$(nrow(removed_units)) treatment observations in refinement")

# Create balance table
using Statistics
for var in covariates
    mean_val = mean(skipmissing(treated_info[!, var]))
    println("\$(var): mean = \$(round(mean_val, digits=3))")
end
```

# Notes
- Uses exact time-unit matching to ensure temporal accuracy
- Handles missing values gracefully through type system
- Maintains sort order for consistent output
- Column names are parameterized to avoid hardcoded dependencies
"""
function obsinfo(
    match_info_df::DataFrame, 
    data::DataFrame, 
    variables::Vector{Symbol};
    full_model_observations::Union{Nothing, Vector{Tuple{Int, Int}}} = nothing,
    time_column::Symbol = :t, 
    id_column::Symbol = :id
)
    # Input validation
    if isempty(match_info_df)
        throw(ArgumentError("Match info DataFrame cannot be empty"))
    end
    if isempty(data)
        throw(ArgumentError("Data DataFrame cannot be empty"))
    end
    
    required_columns = [:timetreated, :treatedunit]
    missing_match_cols = setdiff(required_columns, Symbol.(names(match_info_df)))
    if !isempty(missing_match_cols)
        throw(ArgumentError("Match info DataFrame missing required columns: $(missing_match_cols)"))
    end
    
    required_data_columns = [time_column, id_column]
    append!(required_data_columns, variables)
    missing_data_cols = setdiff(required_data_columns, Symbol.(names(data)))
    if !isempty(missing_data_cols)
        throw(ArgumentError("Data DataFrame missing required columns: $(missing_data_cols)"))
    end
    
    # Extract unique treatment observations from match info
    treatment_observations = unique(match_info_df[!, [:timetreated, :treatedunit]])
    
    # Initialize result DataFrame with proper column structure
    result_columns = [:timetreated, :treatedunit]
    if !isnothing(full_model_observations)
        push!(result_columns, :removed)
    end
    append!(result_columns, variables)
    
    # Create result DataFrame with appropriate types
    observation_info = DataFrame()
    observation_info[!, :timetreated] = treatment_observations[!, :timetreated]
    observation_info[!, :treatedunit] = treatment_observations[!, :treatedunit]
    
    # Add removal indicator if full model comparison requested
    if !isnothing(full_model_observations)
        observation_info[!, :removed] .= false
        
        # Identify observations present in full model but not in current match info
        current_obs_tuples = [(row.timetreated, row.treatedunit) for row in eachrow(treatment_observations)]
        removed_obs = setdiff(full_model_observations, current_obs_tuples)
        
        # Add removed observations to the DataFrame
        for (time_treated, unit_treated) in removed_obs
            push!(observation_info, Dict(
                :timetreated => time_treated,
                :treatedunit => unit_treated,
                :removed => true
            ))
        end
    end
    
    # Initialize covariate columns with appropriate types
    for var in variables
        var_type = eltype(data[!, var])
        observation_info[!, var] = Vector{var_type}(undef, nrow(observation_info))
    end
    
    # Extract covariate values for each treatment observation
    for (row_idx, obs_row) in enumerate(eachrow(observation_info))
        # Find matching row in original data using exact time-unit match
        data_mask = (data[!, time_column] .== obs_row.timetreated) .& 
                   (data[!, id_column] .== obs_row.treatedunit)
        
        matching_rows = findall(data_mask)
        
        if isempty(matching_rows)
            # Handle case where treatment observation not found in data
            for var in variables
                obs_row[var] = missing
            end
        else
            # Extract covariate values (use first match if multiple found)
            data_row_idx = matching_rows[1]
            for var in variables
                obs_row[var] = data[data_row_idx, var]
            end
        end
    end
    
    # Sort results for consistent output
    if !isnothing(full_model_observations)
        sort!(observation_info, [:removed, :timetreated, :treatedunit])
    else
        sort!(observation_info, [:timetreated, :treatedunit])
    end
    
    return observation_info
end