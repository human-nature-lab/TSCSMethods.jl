# construction.jl

"""
    makemodel(
        dat::DataFrame, 
        t::Symbol, 
        id::Symbol, 
        treatment::Symbol, 
        outcome::Union{Symbol, Vector{Symbol}},
        covariates::Vector{Symbol}, 
        timevary::Dict{Symbol, Bool},
        F::UnitRange{Int}, 
        L::UnitRange{Int};
        title::String = "model",
        estimator::String = "ATT"
    ) -> CIC

Construct a CIC (Changes-in-Changes) model for causal inference analysis.

# Arguments
- `dat::DataFrame`: Input data containing all required variables
- `t::Symbol`: Column name for time variable
- `id::Symbol`: Column name for unit identifier
- `treatment::Symbol`: Column name for treatment variable (0/1 coding expected)
- `outcome::Union{Symbol, Vector{Symbol}}`: Column name(s) for outcome variable(s)
- `covariates::Vector{Symbol}`: Vector of column names for matching covariates
- `timevary::Dict{Symbol, Bool}`: Dictionary mapping each covariate to whether it's time-varying
- `F::UnitRange{Int}`: Time periods for treatment effect estimation (post-treatment)
- `L::UnitRange{Int}`: Time periods for pre-treatment matching window
- `title::String`: Optional title for the model (default: "model")
- `estimator::String`: Estimator type (default: "ATT" for Average Treatment Effect on Treated)

# Returns
- `CIC`: A constructed CIC model ready for matching, balancing, and estimation

# Throws
- `ArgumentError`: If input data is empty, required columns are missing, or parameters are invalid

# Examples
```julia
using TSCSMethods, DataFrames

# Create sample data
data = DataFrame(
    time = repeat(1:100, 50),
    unit_id = repeat(1:50, inner=100),
    treated = repeat([fill(1, 10); fill(0, 40)], inner=100),
    outcome = randn(5000),
    covar1 = randn(5000),
    covar2 = randn(5000)
)

# Construct model
model = makemodel(
    data, :time, :unit_id, :treated, :outcome,
    [:covar1, :covar2],
    Dict(:covar1 => false, :covar2 => true),
    51:60, 40:49  # F and L periods
)
```

# References
Based on the methodology from Imai, Kim, and Wang (2021) for matching methods 
with time-series cross-sectional data.
"""
function makemodel(
  dat::DataFrame, 
  t::Symbol, 
  id::Symbol, 
  treatment::Symbol, 
  outcome::Union{Symbol, Vector{Symbol}},
  covariates::Vector{Symbol}, 
  timevary::Dict{Symbol, Bool},
  F::UnitRange{Int}, 
  L::UnitRange{Int};
  title::String = "model",
  estimator::String = "ATT"
)::CIC

  # Input validation
  if nrow(dat) == 0
    throw(ArgumentError("Input data cannot be empty"))
  end
  
  required_cols = [t, id, treatment]
  if outcome isa Symbol
    push!(required_cols, outcome)
  else
    append!(required_cols, outcome)
  end
  append!(required_cols, covariates)
  
  missing_cols = setdiff(required_cols, Symbol.(names(dat)))
  if !isempty(missing_cols)
    throw(ArgumentError("Missing required columns: $(missing_cols)"))
  end
  
  if isempty(F)
    throw(ArgumentError("F (treatment period) cannot be empty"))
  end
  
  if isempty(L)
    throw(ArgumentError("L (pre-treatment period) cannot be empty"))
  end
  
  # Validate timevary dict matches covariates
  if Set(keys(timevary)) != Set(covariates)
    throw(ArgumentError("timevary dict keys must exactly match covariates"))
  end
  
  observations, ids = observe(dat[!, t], dat[!, id], dat[!, treatment]);
  
  # tobs = make_tobsvec(length(observations), length(ids));
  tobs = make_matches(length(observations), length(ids), length(F));

  return CIC(
    title = title,
    id = id,
    t = t,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    timevary = timevary,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobs,
    treatednum = length(observations),
    estimator = estimator
  )
end
