# match!.jl

"""
    default_treatmentcategories(x) -> Int

Default categorization function for treatment histories during matching eligibility.

# Purpose
Categorizes treatment exposure counts in the pre-crossover period to determine matching
eligibility. Units can only be matched if they fall into the same treatment category
during the pre-crossover window, ensuring similar treatment exposure patterns.

# Arguments
- `x`: Total count of treatment exposures in the pre-crossover period for a unit

# Returns
- `0`: Never treated (x == 0)
- `1`: Ever treated (x > 0)

# Usage in Matching
During matching, for each potential match pair:
1. Count treatments in pre-crossover period for treated unit → category A
2. Count treatments in pre-crossover period for potential match → category B  
3. If `treatmentcategories(A) == treatmentcategories(B)`, matching is allowed
4. Otherwise, the match is ineligible

# Default Logic
This binary categorization ensures that:
- Untreated units can only match with other untreated units
- Previously treated units can only match with other previously treated units
- Prevents contamination from units with different treatment exposure histories

# Custom Categories
Users can provide custom categorization functions for more nuanced matching:
```julia
# Example: Allow multiple treatment intensity levels
custom_categories(x) = x <= 1 ? 0 : (x <= 5 ? 1 : 2)
match!(model, data; treatcat = custom_categories)
```
"""
function default_treatmentcategories(x)
    return if x == 0
        0
    else 1
    end
end

"""
    match!(model::AbstractCICModel, dat::DataFrame; treatcat = default_treatmentcategories, exposure = nothing, variancesonly = true)

Perform matching for treatment events using Mahalanobis distance with crossover window restrictions.

# Purpose
Core matching algorithm that identifies eligible control units for each treated observation,
calculates distances using Mahalanobis metric, and ranks potential matches. Implements the
extended TSCS matching methodology with crossover window constraints to ensure proper
counterfactual identification.

# Arguments
- `model`: CIC model containing treated observations and match storage structures
- `dat`: DataFrame with time-series cross-sectional data
- `treatcat`: Function categorizing treatment histories (default: `default_treatmentcategories`)
- `exposure`: Optional exposure variable name for exposure-based matching
- `variancesonly`: Use diagonal covariance matrix (variances only) vs full covariance (default: true)

# Returns
Modified `model` with populated matches containing:
- Eligible match indicators for each treated observation
- Mahalanobis distances to potential matches
- Ranked match preferences for each forward period

# Matching Process
1. **Eligibility Determination**: Uses crossover window logic to determine which units can serve as matches
2. **Distance Calculation**: Computes Mahalanobis distances using sample covariance matrix
3. **Ranking**: Orders potential matches by distance for each forward period (F)

# Crossover Window Logic
For treatment at time `t` with forward period `f`:
- **Pre-crossover period**: `[t + f - max(F), t + f - min(F)]` 
- **Post-crossover period**: `[t + f - min(F), t + f - max(F)]`
- Matches must have similar treatment patterns in pre-crossover period
- Matches cannot be treated in post-crossover period

# Performance Notes
- Uses threaded processing for parallel distance calculations
- Memory-efficient SubArray views via groupindices system
- Handles missing data through multiple dispatch

# Example
```julia
model = makemodel(data, :time, :id, :treatment, :outcome, covariates, timevary, F, L)
match!(model, data)  # Standard matching
match!(model, data; exposure = :policy_exposure)  # With exposure matching
```
"""
function match!(
    model::AbstractCICModel, 
    dat::DataFrame;
    treatcat::Function = default_treatmentcategories,
    exposure::Union{Nothing, Symbol} = nothing,
    variancesonly::Bool = true
)::AbstractCICModel

    # Input validation
    if nrow(dat) == 0
        throw(ArgumentError("Input data cannot be empty"))
    end

    # Extract model components
    (; observations, matches, ids) = model;
    (; F, L, id, t, treatment, covariates) = model

    # Validate that required columns exist in data
    required_columns = [id, t, treatment]
    append!(required_columns, covariates)
    if !isnothing(exposure)
        push!(required_columns, exposure)
    end

    missing_columns = setdiff(required_columns, Symbol.(names(dat)))
    if !isempty(missing_columns)
        throw(ArgumentError("Missing required columns in data: $(missing_columns)"))
    end

    if length(observations) == 0
        @warn "No treatment observations found - matching may not be meaningful"
    end

    # Extract time window parameters
    forward_periods_count = length(F)

    # Bounds for outcome window (forward periods)
    min_forward_period, max_forward_period = extrema(F)

    # Bounds for covariate matching window (lag periods)
    min_lag_period, max_lag_period = extrema(L)

    covariate_count = length(covariates)

    # Extract covariate data matrix
    X = Matrix(dat[!, covariates])

    # Step 1: Determine which units are eligible to serve as matches
    covariate_groups, time_groups = eligibility!(
        matches, observations, X,
        ids, treatcat,
        dat[!, t], dat[!, id], dat[!, treatment],
        min_forward_period, max_forward_period, min_lag_period; exposure = exposure
    );

    # Step 2: Calculate Mahalanobis distances between treated units and potential matches
    distances_allocate!(matches, forward_periods_count, covariate_count)

    # Compute sample covariance matrix for distance calculations
    inverse_covariance_matrix = samplecovar(dat[!, t], X; variancesonly = variancesonly)

    distances_calculate!(
        matches, observations, ids, covariates, covariate_groups, time_groups, 
        min_forward_period, min_lag_period, max_lag_period, inverse_covariance_matrix
    )

    # Step 3: Remove infinite distances (ineligible matches) and rank remaining matches
    # Filter out units that were deemed ineligible (have infinite distances)
    for i in eachindex(matches)
        treatment_observation = @views matches[i]
        # Keep only finite distances (eligible matches)
        finite_distance_mask = .!isinf.(treatment_observation.distances[:, 1])
        matches[i] = @set treatment_observation.distances = treatment_observation.distances[finite_distance_mask, :]
    end

    # Step 4: Rank matches by distance for each forward period
    rank!(matches, forward_periods_count)

    return model
end

"""
    eligibility!(matches, observations, X::Matrix{Float64}, ids, treatcat, dat_t, dat_id, dat_trt, fmin, fmax, Lmin; exposure = nothing)

Determine matching eligibility for all treated observations using crossover window constraints.

# Purpose
Core eligibility determination that implements the TSCS matching methodology's crossover
window logic. Creates efficient indexed views of the data and determines which units
can serve as matches for each treated observation.

# Arguments
- `matches`: Vector of match objects to populate with eligibility information
- `observations`: Vector of treated observations (time, unit_id) tuples
- `X`: Covariate data matrix (Float64 version for non-missing data)
- `ids`: Vector of unique unit identifiers
- `treatcat`: Treatment categorization function for crossover period matching
- `dat_t, dat_id, dat_trt`: Time, unit ID, and treatment vectors from original data
- `fmin, fmax`: Bounds of forward period range (F)
- `Lmin`: Minimum lag period (lower bound of L)
- `exposure`: Optional exposure variable for exposure-based matching

# Returns
- `tg`: Covariate group indices dictionary
- `rg`: Time group indices dictionary

# Process
1. **Group Index Creation**: Builds efficient (treatment_time, unit_id) indexed views
2. **Eligibility Determination**: Calls `eligiblematches!()` to apply crossover window logic
3. **Exposure Handling**: Supports both standard and exposure-based matching

# Crossover Window Logic
For each treated observation, potential matches are evaluated using:
- Pre-crossover treatment similarity (via `treatcat` function)
- Post-crossover treatment exclusion (no treatment allowed in specific windows)
- Time window constraints ensuring proper temporal alignment
"""
function eligibility!(
    matches, observations, X::Matrix{Float64},
    ids, treatcat,
    data_time, data_unit_id, data_treatment,
    min_outcome_period, max_outcome_period, min_matching_period;
    exposure = nothing
)
    # Create efficient indexed views of data organized by (treatment_time, unit_id)
    # This groups observations for fast lookup during matching
    if isnothing(exposure)
        # Standard matching: group by covariate data, time patterns, and treatment status
        covariate_groups, time_groups, treatment_groups = make_groupindices(
            data_time, data_treatment,
            data_unit_id, ids,
            max_outcome_period, min_matching_period,
            X
        )

        # Apply crossover window constraints to determine which units can be matched
        # This enforces the "no treatment in crossover window" constraint
        eligiblematches!(
            observations, matches,
            time_groups, treatment_groups, ids, min_outcome_period, max_outcome_period, treatcat
        )
    else
        # Exposure-based matching: includes exposure history patterns in grouping
        covariate_groups, time_groups, treatment_groups, exposure_groups = make_groupindices(
            data_time, data_treatment,
            data_unit_id, ids,
            max_outcome_period, min_matching_period,
            X;
            exvec = exposure
        )

        # Apply crossover window constraints plus exposure pattern matching
        # Units must match both treatment timing and exposure history
        eligiblematches!(
            observations, matches,
            time_groups, treatment_groups, ids, min_outcome_period, max_outcome_period, treatcat, exg = exposure_groups
        )
    end

    return covariate_groups, time_groups
end

"""
    eligibility!(matches, observations, X::Matrix{Union{Missing, Float64}}, ids, treatcat, dat_t, dat_id, dat_trt, fmin, fmax, Lmin; exposure = nothing)

Determine matching eligibility with support for missing covariate data.

# Purpose
Missing data version of the core eligibility function. Handles datasets where covariates
may contain missing values, using appropriate handling in the groupindices creation and
subsequent matching algorithms.

# Arguments
Same as standard `eligibility!()` but with:
- `X`: Covariate matrix allowing missing values (`Matrix{Union{Missing, Float64}}`)

# Missing Data Handling
- Creates appropriate SubArray types that can handle missing values
- Propagates missing data handling through the matching pipeline
- Distance calculations will appropriately handle missing values in downstream functions

# Implementation Notes
Identical logic to the Float64 version but uses different type dispatch to ensure
proper handling of missing values throughout the matching process.
"""
function eligibility!(
    matches, observations, X::Matrix{Union{Missing, Float64}},
    ids, treatcat,
    data_time, data_unit_id, data_treatment,
    min_outcome_period, max_outcome_period, min_matching_period; exposure = nothing
)

  # Create efficient indexed views of data organized by (treatment_time, unit_id)
  # Handles missing values in covariate matrix through type dispatch
  if isnothing(exposure)
    # Standard matching: group by covariate data, time patterns, and treatment status
    covariate_groups, time_groups, treatment_groups = make_groupindices(
        data_time, data_treatment,
        data_unit_id, ids,
        max_outcome_period, min_matching_period,
        X
    )

    # Apply crossover window constraints to determine which units can be matched
    # Missing values in covariates are handled in downstream distance calculations
    eligiblematches!(
        observations, matches,
        covariate_groups, time_groups, treatment_groups, ids, min_outcome_period, max_outcome_period, treatcat
    )
  else
    # Exposure-based matching: includes exposure history patterns in grouping
    covariate_groups, time_groups, treatment_groups, exposure_groups = make_groupindices(
        data_time, data_treatment,
        data_unit_id, ids,
        max_outcome_period, min_matching_period,
        X;
        exvec = exposure
    )

    # Apply crossover window constraints plus exposure pattern matching
    # Units must match both treatment timing and exposure history
    eligiblematches!(
        observations, matches,
        covariate_groups, time_groups, treatment_groups, ids, min_outcome_period, max_outcome_period, treatcat, exg = exposure_groups
    )
  end

  return covariate_groups, time_groups
end

"""
    eligibility!(matches, observations, X::Matrix{Union{}}, ids, treatcat, dat_t, dat_id, dat_trt, fmin, fmax, Lmin; exposure = nothing)

Handle matching eligibility when no covariates are specified.

# Purpose
Special case handler for models with no covariates (empty covariates array). Converts
the empty matrix to an appropriate Float64 matrix structure and delegates to the
standard eligibility function.

# Arguments
Same as standard `eligibility!()` but with:
- `X`: Empty covariate matrix (`Matrix{Union{}}`)

# Behavior
- Creates a properly sized `Matrix{Float64}` with 0 columns but correct row count
- Delegates to the standard Float64 version of `eligibility!()`
- Enables matching based purely on treatment history patterns without covariate distance

# Use Cases
- Treatment effect estimation based only on temporal patterns
- Matching when covariates are unavailable or not needed
- Sensitivity analysis excluding covariate matching
"""
function eligibility!(
    matches, observations, X::Matrix{Union{}},
    ids, treatcat,
    data_time, data_unit_id, data_treatment,
    min_outcome_period, max_outcome_period, min_matching_period; exposure = nothing
)
  # When no covariates provided, create proper Float64 matrix with correct row count but 0 columns
  # This allows temporal-only matching without covariate distance calculations
    observation_count = length(data_time)
    empty_X = Matrix{Float64}(undef, observation_count, 0)
    return eligibility!(
        matches, observations, empty_X,
        ids, treatcat,
        data_time, data_unit_id, data_treatment,
        min_outcome_period, max_outcome_period, min_matching_period; exposure = exposure
  )
end
