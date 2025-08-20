#=
model object
* the data is kept separately
=#

MatchDist = Vector{Vector{Vector{Float64}}} # = Vector{Vector{Vector{Float64}}}

@with_kw struct Tob
  mus::Matrix{Bool}
  distances::Union{Vector{Matrix{Float64}}, Matrix{Float64}} = Matrix{Float64}(undef, 0, 0)# Vector{Matrix{Float64}}(undef, 0)
  # ranks::Dict{Int, SubArray{Bool, 1, Matrix{Bool}, Tuple{Vector{Int64}, Int64}, false}}
  ranks::Dict{Int, Vector{Int}}
end

@with_kw struct TobC
  mus::Matrix{Bool}
  ranks::Dict{Int, Vector{Int}}
end

@with_kw struct TobR
  mus::Matrix{Bool}
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
- `matches::Vector{Tob}`: Matching results for each treated observation
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
  matches::Vector{Tob}
  # balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = 500
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
- `matches::Vector{Tob}`: Matching results for each treated observation
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
  matches::Vector{Tob}
  # balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobC}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobC}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobR}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobR}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobR}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictNoStrat = GrandDictNoStrat()
  iterations::Int64 = 500
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
  matches::Vector{TobR}
  meanbalances::DataFrame = DataFrame()
  grandbalances::GrandDictStrat = GrandDictStrat()
  iterations::Int64 = 500
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
