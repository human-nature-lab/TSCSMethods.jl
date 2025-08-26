# balancing.jl

"""
    balance!(model::VeryAbstractCICModel, dat::DataFrame) -> VeryAbstractCICModel

Perform covariate balancing on a matched model to ensure treated and control units
are comparable on observed characteristics.

# Arguments
- `model::VeryAbstractCICModel`: A CIC model that has been matched (using `match!`)
- `dat::DataFrame`: Input data containing all model variables

# Returns
- The input model with updated balance statistics

# Description
This function performs two types of balancing:
1. Mean balancing: Computes balance statistics for each covariate across treatment periods
2. Grand balancing: Aggregates balance statistics across the entire model

Balancing assesses how well the matching procedure achieved covariate balance between
treated and control groups. Poor balance may indicate the need for refinement via
calipers or other restrictions.

# Throws
- `ArgumentError`: If input data is empty
- `ErrorException`: If balance calculations fail

# Examples
```julia
# After constructing and matching a model
model = makemodel(data, :time, :id, :treatment, :outcome, 
                 [:covar1, :covar2], timevary_dict, F_periods, L_periods)
match!(model, data)

# Perform balancing
balance!(model, data)

# Check balance results
checkbalances(model, data)
```

# See Also
- [`checkbalances`](@ref): Assess balance quality
- [`autobalance`](@ref): Automatic balance improvement
- [`match!`](@ref): Matching procedure that should precede balancing
"""
function balance!(model::VeryAbstractCICModel, dat::DataFrame)::VeryAbstractCICModel
  # Input validation
  if nrow(dat) == 0
    throw(ArgumentError("Input data cannot be empty"))
  end

  try
    meanbalance!(model, dat);
    grandbalance!(model);
  catch e
    throw(ErrorException("Balance calculation failed: $e"))
  end

  return model
end
