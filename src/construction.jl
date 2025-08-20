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

  # Input validation with user-friendly error messages
  if nrow(dat) == 0
    throw(ArgumentError("Input data cannot be empty. Please provide a DataFrame with your panel data."))
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
    available_cols = join(string.(names(dat)), ", ")
    throw(ArgumentError("""
    Missing required columns in your data: $(join(string.(missing_cols), ", "))
    
    Available columns in your DataFrame: $available_cols
    
    Make sure your data includes:
    - Time variable: $t
    - Unit identifier: $id  
    - Treatment variable: $treatment
    - Outcome variable(s): $(outcome isa Symbol ? outcome : join(string.(outcome), ", "))
    - Covariates: $(join(string.(covariates), ", "))
    """))
  end
  
  if isempty(F)
    throw(ArgumentError("""
    F (post-treatment periods) cannot be empty. 
    
    F specifies which periods after treatment to estimate effects for.
    Example: F = 1:5 estimates effects 1-5 periods after treatment.
    """))
  end
  
  if isempty(L)
    throw(ArgumentError("""
    L (pre-treatment periods) cannot be empty.
    
    L specifies which pre-treatment periods to use for matching.
    IMPORTANT: L periods should be NEGATIVE (e.g., L = -10:-1 for 10 periods before treatment).
    Example: L = -15:-5 uses periods 15-5 before treatment for matching.
    """))
  end
  
  # Check for common mistake: positive L periods
  if any(L .> 0)
    throw(ArgumentError("""
    L periods should be NEGATIVE (pre-treatment periods).
    
    You specified L = $L, but L should represent periods BEFORE treatment.
    Correct examples:
    - L = -10:-1  (uses 10 periods before treatment)
    - L = -20:-5  (uses periods 20-5 before treatment)
    
    This is a common mistake when adapting from other software!
    """))
  end
  
  # Validate timevary dict matches covariates
  if Set(keys(timevary)) != Set(covariates)
    missing_keys = setdiff(Set(covariates), Set(keys(timevary)))
    extra_keys = setdiff(Set(keys(timevary)), Set(covariates))
    
    error_msg = "The timevary dictionary must specify time-varying status for all covariates.\n\n"
    
    if !isempty(missing_keys)
      error_msg *= "Missing from timevary dict: $(join(string.(missing_keys), ", "))\n"
    end
    
    if !isempty(extra_keys)  
      error_msg *= "Extra keys in timevary dict: $(join(string.(extra_keys), ", "))\n"
    end
    
    error_msg *= "\nExample: timevary = Dict(:pop_dens => false, :cumul_death_rate => true)"
    
    throw(ArgumentError(error_msg))
  end
  
  # Validate treatment variable format
  treatment_values = unique(dat[!, treatment])
  if !all(v in [0, 1] for v in treatment_values)
    throw(ArgumentError("""
    Treatment variable '$treatment' must contain only 0 and 1 values.
    
    Found values: $(sort(unique(treatment_values)))
    
    For staggered treatment designs:
    - Use 1 ONLY on the day treatment occurs for each unit
    - Use 0 for all other time periods
    - Do NOT use continuous treatment (1 for all post-treatment periods)
    """))
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
