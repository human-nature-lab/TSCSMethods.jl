# balancing utilities

# Buffer pools for BalanceData optimization
const BALANCE_FLOAT_POOL = Ref{Union{Nothing, Channel{Vector{Float64}}}}(nothing)
const BALANCE_BIT_POOL = Ref{Union{Nothing, Channel{BitVector}}}(nothing)

function __init_balance_pools__()
    if BALANCE_FLOAT_POOL[] === nothing
        pool = Channel{Vector{Float64}}(Threads.nthreads() * 4)
        for _ in 1:(Threads.nthreads() * 4)
            put!(pool, Vector{Float64}(undef, 100))  # Conservative max size
        end
        BALANCE_FLOAT_POOL[] = pool
    end
    
    if BALANCE_BIT_POOL[] === nothing
        pool = Channel{BitVector}(Threads.nthreads() * 4)
        for _ in 1:(Threads.nthreads() * 4)
            put!(pool, BitVector(undef, 100))
        end
        BALANCE_BIT_POOL[] = pool
    end
end

"""
    get_balance_data(size::Int, fill_missing::Bool = true) -> (BalanceData, Bool)

Efficiently allocate BalanceData objects using memory pooling for performance optimization.

# Purpose
Creates BalanceData objects while reusing pre-allocated memory buffers when possible.
This reduces allocation overhead during intensive balance calculations across many
treated observations and matching periods.

# Arguments
- `size`: Number of elements needed in the BalanceData object
- `fill_missing`: Whether to initialize all values as missing (default: true)

# Returns
Tuple of:
- `BalanceData`: The allocated/reused balance data object
- `Bool`: Whether the object came from the pool (true) or was newly allocated (false)

# Details
- For small objects (≤100 elements): attempts to reuse pooled memory
- For large objects or when pool is empty: allocates new memory directly
- Pooled objects should be returned via `return_balance_data()` when no longer needed
- Thread-safe through channel-based pooling system

# Example
```julia
balance_data, is_pooled = get_balance_data(50, true)
# ... use balance_data ...
return_balance_data(balance_data, is_pooled)  # Return to pool when done
```
"""
function get_balance_data(size::Int, fill_missing::Bool = true)
    __init_balance_pools__()
    
    if size > 100 || !isready(BALANCE_FLOAT_POOL[]) || !isready(BALANCE_BIT_POOL[])
        # Pool empty or size too large, allocate directly
        return BalanceData(size, fill_missing), false
    end
    
    values = take!(BALANCE_FLOAT_POOL[])
    is_missing = take!(BALANCE_BIT_POOL[])
    
    resize!(values, size)
    resize!(is_missing, size)
    
    if fill_missing
        fill!(is_missing, true)
    end
    
    return BalanceData(values, is_missing), true
end

"""
    return_balance_data(balance_data::BalanceData, is_pooled::Bool)

Return a BalanceData object to the memory pool for reuse, if it came from the pool.

# Purpose
Completes the memory pooling cycle by returning used BalanceData objects back to
the pool for future reuse. This prevents memory fragmentation and reduces allocation
overhead in subsequent balance calculations.

# Arguments
- `balance_data`: The BalanceData object to potentially return to pool
- `is_pooled`: Whether this object originally came from the pool (from `get_balance_data`)

# Details
- Only returns objects to pool if they originally came from it (`is_pooled == true`)
- Only pools objects with ≤100 elements (larger ones are discarded)
- Thread-safe operation using channel-based pools
- Should be called when BalanceData objects are no longer needed

# Usage Pattern
```julia
balance_data, is_pooled = get_balance_data(50, true)
# ... use balance_data for calculations ...
return_balance_data(balance_data, is_pooled)  # Return when done
```
"""
function return_balance_data(balance_data::BalanceData, is_pooled::Bool)
    if is_pooled && length(balance_data) <= 100
        put!(BALANCE_FLOAT_POOL[], balance_data.values)
        put!(BALANCE_BIT_POOL[], balance_data.is_missing)
    end
end

"""
    compute_treated_std(model::VeryAbstractCICModel, dat::DataFrame) -> Dict{Tuple{Int64, Symbol}, Float64}

Compute standardization factors (1/σ) for covariates across all treated units' matching periods.

# Purpose
For balance calculations, we need to standardize covariate differences by the standard deviation
of treated units during their matching windows. This ensures balance statistics are comparable
across covariates with different scales.

# Returns
Dictionary with keys `(time_offset, covariate)` and values `1/σ`, where:
- `time_offset`: Days relative to treatment (negative = pre-treatment, positive = post-treatment)  
- `covariate`: Covariate name
- `1/σ`: Inverse standard deviation for standardization

# Details
For each treated unit at treatment time `t`, we collect covariate values from their matching
window `[t + min(L), t + max(F)]`. We then compute the standard deviation of each covariate
at each time offset across all treated units and return the inverse for efficient multiplication
during balance calculations.

# Example
```julia
std_factors = compute_treated_std(model, data)
# std_factors[(-10, :population)] = 0.05  # 1/σ for population 10 days before treatment
```
"""
function compute_treated_std(model::VeryAbstractCICModel, dat::DataFrame)

  treated_observations = unique(dat[dat[!, model.treatment] .== 1, [model.t, model.id]])
  sort!(treated_observations, [model.t, model.id])
  observation_indices = Int[];
  time_offsets = Int[];

  min_lag_period = minimum(model.L)
  max_forward_period = maximum(model.F)

  for treated_obs in eachrow(treated_observations)

    unit_condition = dat[!, model.id] .== treated_obs[2];
    time_condition = (dat[!, model.t] .>= (treated_obs[1] - 1 + min_lag_period)) .& (dat[!, model.t] .<= treated_obs[1] + max_forward_period);

    append!(observation_indices, findall((unit_condition .& time_condition)))
    # time_offsets relative to treatment time, not the actual time
    append!(time_offsets, dat[unit_condition .& time_condition, model.t] .- treated_obs[1])
  end

  covariate_values = @view dat[observation_indices, model.covariates];

  if any([Missing <: eltype(c) for c in eachcol(covariate_values)])
    unique_offsets = unique(time_offsets);
    standardization_factors = zeros(Float64, length(unique_offsets));
    # standardization_factors = Dict{Tuple{Int64, Symbol}, Union{Float64, Missing}}()
    standardization_factors = Dict{Tuple{Int64, Symbol}, Float64}()

    for offset in unique_offsets
      timepoint_values = @view covariate_values[time_offsets .== offset, :]
      for covar in model.covariates
        standardization_factors[(offset, covar)] = 1.0 / std(
          skipmissing(timepoint_values[!, covar]); corrected = true
        )
      end
    end
  else
    unique_offsets = unique(time_offsets);
    standardization_factors = zeros(Float64, length(unique_offsets));
    standardization_factors = Dict{Tuple{Int64, Symbol}, Float64}()
    
    for offset in unique_offsets
      timepoint_values = @view covariate_values[time_offsets .== offset, :]
      for covar in model.covariates
        standardization_factors[(offset, covar)] = 1.0 / std(timepoint_values[!, covar]; corrected = true)
      end
    end
  end
  
  return standardization_factors 
end

"""
    initialize_balance_storage!(model)

Prepare the meanbalances DataFrame storage for balance calculations.

# Purpose
Allocates and initializes the meanbalances DataFrame which stores balance statistics
for each treated observation. The DataFrame structure depends on whether covariates
are time-varying or static.

# Details
- Creates one row per treated observation
- For time-varying covariates: Vector{Vector{BalanceData}} (periods × time points)
- For static covariates: Vector{BalanceData} (just periods)
- Uses object pooling for memory efficiency
"""
function initialize_balance_storage!(model)

  (; observations, matches, ids, meanbalances,
   covariates, timevary,
   t, id, treatment,
   F, L) = model


  lag_periods_count = length(L); forward_periods_count = length(F);

  # meanbalances = DataFrame(
  #   treattime = [ob[1] for ob in model.observations],
  #   treatunit = [ob[2] for ob in model.observations],
  #   matchunitsets = [mm for mm in model.matchunits]
  # )

  # need to check observations to see if there are any matches left

  meanbalances[!, :fs] = Vector{Vector{Bool}}(undef, length(matches));

  for covar in covariates
    if timevary[covar]
      meanbalances[!, covar] = Vector{Vector{BalanceData}}(undef, length(matches))
    else 
      meanbalances[!, covar] = Vector{BalanceData}(undef, length(matches))
    end
  end

  _fill_meanbalances!(
    meanbalances, matches, lag_periods_count, covariates, timevary, forward_periods_count
  );

  return model
end

"""
    _fill_meanbalances!(meanbalances, matches, lag_periods_count, covariates, timevary, forward_periods_count)

Populate meanbalances DataFrame with properly sized BalanceData objects for each treated observation.

# Purpose
For each treated observation, determines which forward periods (F values) have at least one
eligible match, then allocates appropriate BalanceData storage structures. The allocation
pattern depends on whether covariates are time-varying or static.

# Arguments
- `meanbalances`: DataFrame to populate (one row per treated observation)
- `matches`: Vector of match objects containing eligible_matches matrices
- `lag_periods_count`: Number of lag periods (length of L range)
- `covariates`: Vector of covariate names
- `timevary`: Dict indicating which covariates are time-varying
- `forward_periods_count`: Number of forward periods (length of F range)

# Details
- Creates `:fs` column indicating which forward periods have matches
- For time-varying covariates: Vector{BalanceData} of length periods_with_matches
- For static covariates: Single BalanceData of length periods_with_matches
- Uses memory pooling via `get_balance_data()` for efficiency
- Each BalanceData object sized according to lag_periods_count for time-varying covariates

# Storage Structure
```
meanbalances[i, :covariate] = 
  - Time-varying: [BalanceData(lag_periods), BalanceData(lag_periods), ...]
  - Static: BalanceData(periods_with_matches)
```
"""
function _fill_meanbalances!(
  meanbalances, matches, lag_periods_count, covariates, timevary, forward_periods_count
)

  for (i, balance_row) in enumerate(eachrow(meanbalances))
    
    # we want to create a vector for f, so long as at least one match unit allows it
    balance_row[:fs] = Vector{Bool}(undef, forward_periods_count);
    find_periods_with_matches!(balance_row[:fs], matches[i].eligible_matches);
    periods_with_matches = sum(balance_row[:fs]);
  
    __fill_meanbalances!(
      balance_row, periods_with_matches, lag_periods_count, covariates, timevary
    )
  end
  return meanbalances
end

"""
    __fill_meanbalances!(balance_row, periods_with_matches, lag_periods_count, covariates, timevary)

Allocate BalanceData objects for a single treated observation with error handling and memory pooling.

# Purpose
Handles the actual allocation of BalanceData objects for one row of the meanbalances DataFrame.
Implements error-safe memory pooling with proper cleanup if allocation fails partway through.

# Arguments
- `balance_row`: Single row from meanbalances DataFrame to populate
- `periods_with_matches`: Number of forward periods that have eligible matches
- `lag_periods_count`: Number of lag periods (for time-varying covariate sizing)
- `covariates`: Vector of covariate names to allocate storage for
- `timevary`: Dict indicating which covariates are time-varying

# Details
- For time-varying covariates: Creates Vector{BalanceData} with `periods_with_matches` elements, each of size `lag_periods_count`
- For static covariates: Creates single BalanceData of size `periods_with_matches`
- Uses memory pooling via `get_balance_data()` for allocation efficiency
- Implements error handling: if allocation fails, returns all pooled objects to prevent memory leaks
- Tracks all pooled allocations for proper cleanup

# Error Handling
If any allocation fails, all previously allocated pooled objects are automatically
returned to their respective pools to prevent memory leaks.
"""
function __fill_meanbalances!(
  balance_row, periods_with_matches, lag_periods_count, covariates, timevary
)
  pooled_data = []  # Track pooled data for cleanup
  
  try
    for covar in covariates
      if timevary[covar]
        covariate_balance_data = BalanceData[]
        for _ in 1:periods_with_matches
          balance_data, is_pooled = get_balance_data(lag_periods_count, true)
          push!(covariate_balance_data, balance_data)
          push!(pooled_data, (balance_data, is_pooled))
        end
        balance_row[covar] = covariate_balance_data
      else
        balance_data, is_pooled = get_balance_data(periods_with_matches, true)
        balance_row[covar] = balance_data
        push!(pooled_data, (balance_data, is_pooled))
      end
    end
  catch e
    # Clean up any allocated data on error
    for (balance_data, is_pooled) in pooled_data
      return_balance_data(balance_data, is_pooled)
    end
    rethrow(e)
  end
  
  return balance_row
end

"""
    find_periods_with_matches!(has_matches_by_period, eligible_matches_matrix)

Determine which forward periods have at least one eligible match across all potential match units.

# Purpose
For a given treated observation, creates a boolean vector indicating which forward periods (F values)
have at least one unit eligible to serve as a match. This is used to determine which periods
need BalanceData storage allocation.

# Arguments
- `has_matches_by_period`: Output boolean vector to populate (length = number of F periods)
- `eligible_matches_matrix`: Boolean matrix where `[unit, period]` indicates if unit is eligible match for that period

# Details
- Performs column-wise `any()` operation across the eligible_matches_matrix
- `has_matches_by_period[i] = true` if any unit can match in forward period i
- `has_matches_by_period[i] = false` if no units are eligible for forward period i
- Used to optimize storage: only allocate BalanceData for periods with potential matches

# Example
```julia
# eligible_matches_matrix: 3 units × 4 periods
# [true  false true  true ]  # unit 1
# [false true  false false]  # unit 2  
# [false false false true ]  # unit 3
# 
# Result: [true, true, true, true] - all periods have ≥1 eligible match
```
"""
function find_periods_with_matches!(has_matches_by_period, eligible_matches_matrix)
  for period_index in eachindex(has_matches_by_period)
    has_matches_by_period[period_index] = any(eligible_matches_matrix[:, period_index])
  end
end
