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

"""
    MissingData{T}

Efficient storage for data with missing values using separate arrays for values and missing indicators.
This replaces Vector{Union{Missing, T}} with better performance and memory pooling capability.
"""
struct MissingData{T}
    values::Vector{T}
    is_missing::BitVector
end

# Type aliases for clarity
const BalanceData = MissingData{Float64}
const DistanceData = MissingData{Float64}

# Constructor that mimics the old Vector{Union{Missing, T}} interface
function MissingData{T}(size::Int, fill_missing::Bool = true) where T
    values = Vector{T}(undef, size)
    is_missing = BitVector(undef, size)
    if fill_missing
        fill!(is_missing, true)  # Start with all missing
    end
    return MissingData{T}(values, is_missing)
end

# Convenience constructors for type aliases (use generic constructor)
# Note: Using function definitions to avoid method overwriting during precompilation

# Indexing interface to maintain compatibility
Base.getindex(md::MissingData{T}, i::Int) where T = md.is_missing[i] ? missing : md.values[i]
function Base.setindex!(md::MissingData{T}, val::Union{Missing, T}, i::Int) where T
    if ismissing(val)
        md.is_missing[i] = true
    else
        md.is_missing[i] = false
        md.values[i] = val
    end
end

Base.length(md::MissingData) = length(md.values)
Base.eachindex(md::MissingData) = eachindex(md.values)

# Iteration interface - critical for skipmissing and other operations
Base.iterate(md::MissingData) = length(md) == 0 ? nothing : (md[1], 1)
Base.iterate(md::MissingData, state::Int) = state >= length(md) ? nothing : (md[state + 1], state + 1)

# Size and axes for array-like behavior
Base.size(md::MissingData) = size(md.values)
Base.axes(md::MissingData) = axes(md.values)

# Arithmetic operations needed for Statistics functions
Base.:/(md::MissingData{T}, n::Number) where T = 
    MissingData{T}(md.values ./ n, copy(md.is_missing))

Base.:+(md1::MissingData{T}, md2::MissingData{T}) where T = 
    MissingData{T}(md1.values .+ md2.values, md1.is_missing .| md2.is_missing)

Base.:-(md1::MissingData{T}, md2::MissingData{T}) where T = 
    MissingData{T}(md1.values .- md2.values, md1.is_missing .| md2.is_missing)

Base.:*(md::MissingData{T}, n::Number) where T = 
    MissingData{T}(md.values .* n, copy(md.is_missing))

Base.zero(::Type{MissingData{T}}) where T = 
    MissingData{T}(T[], BitVector())

# Support for sum and mean operations
Base.sum(md::MissingData{T}) where T = 
    isempty(md.values) ? zero(T) : sum(md.values[.!md.is_missing])

# Mean function (will be extended in overallbalancing.jl where Statistics is available)
function _missing_data_mean(md::MissingData{T}) where T
    valid_indices = .!md.is_missing
    n_valid = sum(valid_indices)
    n_valid == 0 ? convert(T, NaN) : sum(md.values[valid_indices]) / n_valid
end
