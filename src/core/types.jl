#=
model object
* the data is kept separately
=#

const DEFAULT_BOOTSTRAP_ITERATIONS = 500

MatchDist = Vector{Vector{Vector{Float64}}} # = Vector{Vector{Vector{Float64}}}

"""
    TreatmentObservationMatches

Stores matching results for each treated observation in basic CIC models.

Contains the core data structures for matching:
- `eligible_matches`: Boolean matrix indicating which control units are eligible matches
- `distances`: Computed distances between treated and control units  
- `match_rankings`: Ranked preferences for matches by time period

This is the fundamental data structure for storing match relationships in TSCS designs.
"""
@with_kw struct TreatmentObservationMatches
  eligible_matches::Matrix{Bool}
  distances::Union{Vector{Matrix{Float64}}, Matrix{Float64}} = Matrix{Float64}(undef, 0, 0)
  match_rankings::Dict{Int, Vector{Int}}
end

"""
    TreatmentObservationCaliperMatches

Stores matching results for caliper-constrained CIC models.

Similar to TreatmentObservationMatches but for models where matches are restricted
by caliper constraints (maximum allowable distance thresholds).
"""
@with_kw struct TreatmentObservationCaliperMatches
  eligible_matches::Matrix{Bool}
  match_rankings::Dict{Int, Vector{Int}}
end

"""
    TreatmentObservationRefinedMatches  

Stores matching results for refined CIC models.

Used in refined matching procedures where the initial match set is 
iteratively improved through additional matching criteria.
"""
@with_kw struct TreatmentObservationRefinedMatches
  eligible_matches::Matrix{Bool}
end

# Vector{Vector{MVector{C}{Float64}}}([[[1],[1]],[[1],[1]]])

StratDict = Dict{Int64, Dict{Symbol, Float64}};
GrandDictNoStrat = Dict{Symbol, Union{Vector{Float64}, Float64}};
GrandDictStrat = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}};

abstract type VeryAbstractCICModel end

abstract type AbstractCICModel <: VeryAbstractCICModel end

abstract type AbstractCICModelStratified <: VeryAbstractCICModel end

"""
    CIC <: AbstractCICModel

Main model type for causal inference with time-series cross-sectional data.

Contains all information needed for matching, balancing, and estimation in TSCS designs.
Created by `makemodel()` and used throughout the analysis workflow.

# Fields
- `title::String`: Model title for identification
- `id::Symbol`: Column name for unit identifier 
- `t::Symbol`: Column name for time variable
- `outcome::Union{Symbol,Vector{Symbol}}`: Outcome variable(s)
- `treatment::Symbol`: Treatment variable (binary 0/1)
- `covariates::Vector{Symbol}`: Covariates used for matching
- `timevary::Dict{Symbol, Bool}`: Whether each covariate is time-varying
- `reference::Int`: Reference time period (default: -1)
- `F::UnitRange{Int}`: Post-treatment periods for estimation
- `L::UnitRange{Int}`: Pre-treatment periods for matching (negative values)
- `observations::Vector{Tuple{Int, Int}}`: Treated observations (time, unit)
- `ids::Vector{Int}`: All unit identifiers
- `matches::Vector{TreatmentObservationMatches}`: Matching results for each treated observation
- `meanbalances::DataFrame`: Covariate balance statistics
- `grandbalances::Dict`: Overall balance measures
- `iterations::Int`: Bootstrap iterations (default: 500)
- `results::DataFrame`: Treatment effect estimates
- `treatednum::Int`: Number of treated observations
- `estimator::String`: Estimator type (default: "ATT")

# Examples
```julia
using TSCSMethods, DataFrames

data = example_data()
model = makemodel(data, :day, :fips, :gub, :death_rte, 
                 [:pop_dens], Dict(:pop_dens => false),
                 1:5, -10:-1)
```
"""
@with_kw struct CIC <: AbstractCICModel
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}# Vector{Tuple{Int, Int}}
  ids::Vector{Int} #Vector{Int}
  matches::Vector{TreatmentObservationMatches}
  # balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Int64
  estimator::String = "ATT"
end

"""
    CICStratified <: AbstractCICModelStratified

Stratified model type for causal inference with time-series cross-sectional data.

Extends the basic CIC model to handle stratified analysis where matching and estimation
are performed separately within subgroups defined by a stratifying variable.

# Fields
- `title::String`: Model title for identification
- `id::Symbol`: Column name for unit identifier
- `t::Symbol`: Column name for time variable  
- `outcome::Union{Symbol,Vector{Symbol}}`: Outcome variable(s)
- `treatment::Symbol`: Treatment variable (binary 0/1)
- `covariates::Vector{Symbol}`: Covariates used for matching
- `timevary::Dict{Symbol, Bool}`: Whether each covariate is time-varying
- `stratifier::Symbol`: Variable used for stratification
- `strata::Vector{Int}`: Values of stratifying variable
- `reference::Int`: Reference time period (default: -1)  
- `F::UnitRange{Int}`: Post-treatment periods for estimation
- `L::UnitRange{Int}`: Pre-treatment periods for matching (negative values)
- `observations::Vector{Tuple{Int, Int}}`: Treated observations (time, unit)
- `ids::Vector{Int}`: All unit identifiers
- `matches::Vector{TreatmentObservationMatches}`: Matching results for each treated observation
- `meanbalances::DataFrame`: Covariate balance statistics by stratum
- `grandbalances::Dict`: Overall balance measures by stratum
- `iterations::Int`: Bootstrap iterations (default: 500)
- `results::DataFrame`: Treatment effect estimates by stratum
- `treatednum::Dict{Int64, Int64}`: Number of treated observations per stratum
- `estimator::String`: Estimator type (default: "ATT")
- `labels::Dict{Int64, String}`: Human-readable labels for strata

# Examples
```julia
# Create stratified model
model = makemodel(data, :day, :fips, :gub, :death_rte,
                 [:pop_dens], Dict(:pop_dens => false),
                 1:5, -10:-1)
strat_model = stratify(model, data, :region)
```
"""
@with_kw struct CICStratified <: AbstractCICModelStratified
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  strata::Vector{Int} = Vector{Int}()
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationMatches}
  # balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Dict{Int64, Int64}
  estimator::String = "ATT"
  labels::Dict{Int64, String}
end

@with_kw struct CaliperCIC <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationCaliperMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Int64
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CIC}
end

@with_kw struct CaliperCICStratified <: AbstractCICModelStratified
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  strata::Vector{Int} = Vector{Int}()
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationCaliperMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Dict{Int64, Int64}
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CICStratified}
  labels::Dict{Int64, String}
end

@with_kw struct RefinedCIC <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  refinementnum::Int
  timevary::Dict{Symbol, Bool}
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationRefinedMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Int64
  estimator::String = "ATT"
  fullmod = Base.RefValue{CIC}
end

@with_kw struct RefinedCICStratified <: AbstractCICModelStratified
  title::String = "refined"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  refinementnum::Int
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  strata::Vector{Int} = Vector{Int}()
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationRefinedMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Int64
  estimator::String = "ATT"
  fullmod = Base.RefValue{CIC}
  labels::Dict{Int64, String}
end

@with_kw struct RefinedCaliperCIC <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  refinementnum::Int
  timevary::Dict{Symbol, Bool}
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationRefinedMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Int64
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CaliperCIC}
end

@with_kw struct RefinedCaliperCICStratified <: AbstractCICModelStratified
  title::String = "refined"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  refinementnum::Int
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  strata::Vector{Int} = Vector{Int}()
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matches::Vector{TreatmentObservationRefinedMatches}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = DEFAULT_BOOTSTRAP_ITERATIONS
  results::DataFrame = DataFrame()
  treatednum::Dict{Int64, Int64}
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CIC}
  labels::Dict{Int64, String}
end

# save a small-size copy of the model essentials
@with_kw struct ModelRecord
  title::String = "model record"
  id::Symbol
  t::Symbol
  outcome::Union{Symbol,Vector{Symbol}}
  treatment::Symbol
  covariates::Vector{Symbol}
  stratifier::Symbol = Symbol()
  strata::Vector{Int} = Vector{Int}()
  reference::Int
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  grandbalances::Union{GrandDictNoStrat, GrandDictStrat}
  iterations::Int64
  results::DataFrame
  treatednum::Union{Int64, Dict{Int64, Int64}}
  estimator::String
  caliper::Dict{Symbol, Float64} = Dict{Symbol, Float64}()
  labels::Dict{Int64, String} = Dict{Int64, String}()
  refinementnumber::Int = Int()
end

"""
        Fblock

Holds the relevant information for bootstrapping and estimation,
for a specific f in the outcome window (an f in a stratum when the model
is stratified).
"""
struct Fblock
    f::Int
    matchunits::Vector{Int}
    weightedoutcomes::Vector{Float64}
    weightedrefoutcomes::Vector{Float64}
    treatment::Vector{Bool}
end

"""
    BalanceData

Efficient storage for balance computation data with separate arrays for values and missing indicators.
This replaces Vector{Union{Missing, Float64}} with better performance and memory pooling capability.
"""
struct BalanceData
    values::Vector{Float64}
    is_missing::BitVector
end

# Constructor that mimics the old Vector{Union{Missing, Float64}} interface
function BalanceData(size::Int, fill_missing::Bool = true)
    values = Vector{Float64}(undef, size)
    is_missing = BitVector(undef, size)
    if fill_missing
        fill!(is_missing, true)  # Start with all missing
    end
    return BalanceData(values, is_missing)
end

# Indexing interface to maintain compatibility
Base.getindex(bd::BalanceData, i::Int) = bd.is_missing[i] ? missing : bd.values[i]
function Base.setindex!(bd::BalanceData, val::Union{Missing, Float64}, i::Int)
    if ismissing(val)
        bd.is_missing[i] = true
    else
        bd.is_missing[i] = false
        bd.values[i] = val
    end
end

Base.length(bd::BalanceData) = length(bd.values)
Base.eachindex(bd::BalanceData) = eachindex(bd.values)

# Iteration interface for skipmissing and Statistics functions
Base.iterate(bd::BalanceData) = length(bd) == 0 ? nothing : (bd[1], 1)
Base.iterate(bd::BalanceData, state::Int) = state >= length(bd) ? nothing : (bd[state + 1], state + 1)

# Size and axes for array-like behavior
Base.size(bd::BalanceData) = size(bd.values)
Base.axes(bd::BalanceData) = axes(bd.values)

# Arithmetic operations for balance calculations
Base.:/(bd::BalanceData, n::Number) = 
    BalanceData(bd.values ./ n, copy(bd.is_missing))

Base.:+(bd1::BalanceData, bd2::BalanceData) = 
    BalanceData(bd1.values .+ bd2.values, bd1.is_missing .| bd2.is_missing)

Base.:-(bd1::BalanceData, bd2::BalanceData) = 
    BalanceData(bd1.values .- bd2.values, bd1.is_missing .| bd2.is_missing)

Base.:*(bd::BalanceData, n::Number) = 
    BalanceData(bd.values .* n, copy(bd.is_missing))

# Support for Statistics functions
Base.sum(bd::BalanceData) = 
    isempty(bd.values) ? 0.0 : sum(bd.values[.!bd.is_missing])
