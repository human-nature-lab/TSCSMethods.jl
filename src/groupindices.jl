# groupindices.jl

"""
Primary Purpose: Balance Calculations

The groupindices are heavily used in balance calculations:
- compute_treated_std() uses them to calculate standardization factors
- meanbalancing.jl uses tg, rg extensively for balance computations
- The indexed views enable efficient access to covariate values during balance assessment

Secondary Purpose: Matching Support

They also support the matching process:
- retrieve_matches.jl uses them to access treatment histories (trtg)
- The exposure indices (exg) are specifically for exposure-based matching
- Time window access during eligibility determination

Why Both?

The groupindices create efficient indexed views of the data organized by (treatment_time, 
unit_id) pairs. This structure is needed because:

1. Matching needs to check treatment histories and eligibility within time windows
2. Balance calculations need fast access to covariate values for treated vs matched units during
  the same time windows

Core Insight

Both matching and balance calculations operate on the same conceptual structure: unit-specific 
data within treatment-relevant time windows. The groupindices provide an efficient way to access
  this structure without repeatedly filtering the original data.

So they're a shared infrastructure that optimizes data access for both phases of the TSCS
methodology - first for determining eligible matches, then for assessing the quality of those
matches through balance statistics.
"""

# types for group indices object
STypeMat = SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false};

STypeMatMis = SubArray{Union{Missing, Float64}, 2, Matrix{Union{Missing, Float64}}, Tuple{Vector{Int64}, Base.Slice{Base.OneTo{Int64}}}, false}

STypeVec = SubArray{Int64, 1, Vector{Int64}, Tuple{Vector{Int64}}, false}
STypeVecBool = SubArray{Bool, 1, Vector{Bool}, Tuple{Vector{Int64}}, false}
# last entry probably depends on original datatype before conversion

"""
    make_groupindices(tvec, treatvec, idvec, uid, fmax, Lmin, X; exvec = nothing)

Create efficient lookup dictionaries for accessing unit-specific data during matching windows.

# Purpose
Builds indexed views of the data organized by (treatment_time, unit_id) pairs, covering
the relevant time windows for matching. This enables fast access to covariate values,
treatment histories, and time vectors for each unit during balance calculations.

# Arguments
- `tvec`: Time vector from the data
- `treatvec`: Treatment indicator vector (0/1 or boolean)
- `idvec`: Unit identifier vector  
- `uid`: Vector of unique unit IDs to process
- `fmax`: Maximum forward period (upper bound of F range)
- `Lmin`: Minimum lag period (lower bound of L range) 
- `X`: Covariate data matrix
- `exvec`: Optional exposure vector for exposure-based matching

# Returns
Tuple of dictionaries with keys `(treatment_time, unit_id)`:
- `tidx`: Covariate values for unit during time window `[tt + Lmin, tt + fmax)`
- `ridx`: Time values for unit during the window
- `tridx`: Treatment indicators for unit during the window
- `exidx`: Exposure values (only if `exvec` provided)

# Details
For each treated observation at time `tt`, creates indexed views for all units covering
the time window from `tt + Lmin` to `tt + fmax - 1`. This spans both the lag periods
(for baseline covariate measurement) and forward periods (for outcome measurement).

# Performance Notes
- Uses threaded processing for parallel index construction
- Returns memory-efficient SubArray views rather than copying data
- Pre-allocates dictionaries with appropriate size hints

# Example
```julia
tidx, ridx, tridx = make_groupindices(data.time, data.treatment, data.unit_id, 
                                      unique_units, 40, -30, covariate_matrix)
# Access unit 1's covariates during treatment time 100's window:
unit1_covariates = tidx[(100, 1)]  # SubArray view of relevant data
```
"""
function make_groupindices(
  tvec, treatvec, idvec, uid, fmax, Lmin, X;
  exvec = nothing
)

  treatvecbool = convert(Vector{Bool}, treatvec);
  tts = sort(unique(tvec[treatvecbool]));
  
  XT = if Missing <: eltype(X)
    STypeMatMis
  else
    STypeMat
  end

  tidx = Dict{Tuple{Int, Int}, XT}()

  ridx = Dict{Tuple{Int, Int}, STypeVec}();
  tridx = Dict{Tuple{Int, Int}, STypeVecBool}();
  exidx = Dict{Tuple{Int, Int}, STypeVec}();

  sizehint!(tidx, length(Iterators.product(tts, uid)));
  sizehint!(ridx, length(Iterators.product(tts, uid)));
  sizehint!(tridx, length(Iterators.product(tts, uid)));
  sizehint!(exidx, length(Iterators.product(tts, uid)));

  if isnothing(exvec)
    _makegroupindices(
      tidx, ridx, tridx,
      tts, uid, fmax, Lmin,
      tvec, idvec,
      treatvecbool, X
    )

    return tidx, ridx, tridx
  else
    _makegroupindices(
      tidx, ridx, tridx,
      tts, uid, fmax, Lmin,
      tvec, idvec,
      treatvecbool, X,
      exidx, exvec
    )
    return tidx, ridx, tridx, exidx
  end
  
end

"""
    _makegroupindices(tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool, X)

Core implementation for building group indices without exposure data.

# Purpose
Populates the index dictionaries by creating SubArray views for each (treatment_time, unit)
combination. Uses parallel processing to efficiently handle large datasets.

# Arguments
- `tidx`, `ridx`, `tridx`: Pre-allocated dictionaries to populate
- `tts`: Vector of treatment times (when treatments occurred)  
- `uid`: Vector of unique unit IDs
- `fmax`, `Lmin`: Time window bounds
- `tvec`, `idvec`, `treatvecbool`: Original data vectors
- `X`: Covariate data matrix

# Implementation Details
- Uses `Threads.@threads :greedy` for parallel processing across (treatment_time, unit) pairs
- For each pair, calls `getyes!()` to determine which data rows fall within the time window
- Creates SubArray views using `@views` macro for memory efficiency
- Thread-safe due to non-overlapping dictionary key assignments
"""
function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool,
  X
)

  # Pre-compute length once to avoid repeated calls
  tvec_len = length(tvec)
  
  # Pre-allocate thread-local buffers for modern threading
  let products = collect(Iterators.product(tts, uid))
    Threads.@threads :greedy for (tt, unit) in products
      yesrows = Vector{Bool}(undef, tvec_len)
      getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)

      tidx[(tt, unit)] = @views X[yesrows, :];
      ridx[(tt, unit)] = @views tvec[yesrows];
      tridx[(tt, unit)] = @views treatvecbool[yesrows];
    end
  end
  return tidx, ridx, tridx
end

"""
    _makegroupindices(tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool, X, exidx, exvec)

Core implementation for building group indices with exposure data support.

# Purpose  
Extended version that additionally creates indexed views of exposure data alongside
the standard covariate, time, and treatment indices. Used for exposure-based matching
where treatment assignment depends on external exposure variables.

# Arguments
- `tidx`, `ridx`, `tridx`: Standard index dictionaries to populate
- `exidx`: Additional dictionary for exposure data views
- `exvec`: Exposure data vector
- Other arguments same as standard version

# Details
Identical logic to the standard version but additionally populates `exidx[(tt, unit)]`
with SubArray views of the exposure data for each (treatment_time, unit) combination.
This enables efficient access to exposure histories during matching algorithms that
consider exposure patterns.
"""
function _makegroupindices(
  tidx, ridx, tridx, tts, uid, fmax, Lmin, tvec, idvec, treatvecbool,
  X, exidx, exvec
)
  # Pre-compute length once to avoid repeated calls
  tvec_len = length(tvec)
  
  # Pre-allocate thread-local buffers for modern threading
  let products = collect(Iterators.product(tts, uid))
    Threads.@threads :greedy for (tt, unit) in products
      yesrows = Vector{Bool}(undef, tvec_len)
      getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)

      tidx[(tt, unit)] = @views X[yesrows, :];
      ridx[(tt, unit)] = @views tvec[yesrows];
      tridx[(tt, unit)] = @views treatvecbool[yesrows];
      exidx[(tt, unit)] = @views exvec[yesrows];
    end
  end
  return tidx, ridx, tridx, exidx
end

"""
    getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)

Determine which data rows fall within the matching time window for a specific unit and treatment time.

# Purpose
Populates a boolean mask indicating which rows of the original data should be included
in the indexed view for a given (treatment_time, unit) combination. This implements the
core logic for defining matching windows in TSCS designs.

# Arguments
- `yesrows`: Output boolean vector to populate (same length as `tvec`)
- `tvec`: Time vector from original data
- `idvec`: Unit ID vector from original data  
- `tt`: Treatment time (when treatment occurred)
- `fmax`: Maximum forward period (upper bound of F range)
- `Lmin`: Minimum lag period (lower bound of L range, typically negative)
- `unit`: Unit ID to extract data for

# Time Window Logic
Selects rows where:
1. `tvec[k] >= tt + Lmin` (at or after start of lag period)
2. `tvec[k] < tt + fmax` (before end of forward period)  
3. `idvec[k] == unit` (belongs to the specified unit)

# Example Time Windows
```
Treatment at tt=100, Lmin=-30, fmax=40:
- Window: [70, 139] (from 100-30 to 100+40-1)
- Includes: lag periods for covariate measurement + forward periods for outcomes
- Excludes: data outside this window or from other units
```

# Performance Notes
- Operates in-place on `yesrows` for memory efficiency
- Used within threaded loops, so must be thread-safe (which it is)
"""
function getyes!(yesrows, tvec, idvec, tt, fmax, Lmin, unit)
  # Core matching window logic: select data rows within the time window for this unit
  for k in eachindex(tvec)
    yesrows[k] = ((tvec[k] < tt + fmax) & (tvec[k] >= tt + Lmin)) & (idvec[k] == unit)
  end
  return yesrows
end
