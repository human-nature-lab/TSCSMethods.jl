#=
model object
* the data is kept separately
=#

MatchDist = Vector{Vector{Vector{Float64}}} # = Vector{Vector{Vector{Float64}}}

@with_kw struct Tob
  mus::Matrix{Bool}
  fs::Vector{Vector{Bool}}
  distances::Vector{Matrix{Float64}} = Vector{Matrix{Float64}}(undef, 0)
  ranks::Dict{Int, SubArray{Bool, 1, Vector{Bool}, Tuple{Vector{Int64}}, false}}
end

@with_kw struct TobC
  mus::Matrix{Bool}
  fs::Vector{Vector{Bool}}
  ranks::Dict{Int, SubArray{Bool, 1, Vector{Bool}, Tuple{Vector{Int64}}, false}}
end

@with_kw struct TobR
  mus::Matrix{Bool}
  fs::Vector{Vector{Bool}}
end

# Vector{Vector{MVector{C}{Float64}}}([[[1],[1]],[[1],[1]]])

StratDict = Dict{Int64, Dict{Symbol, Float64}};
GrandDictNoStrat = Dict{Symbol, Union{Vector{Float64}, Float64}};
GrandDictStrat = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}};

abstract type AbstractCICModel end

abstract type AbstractCICModelStratified <: AbstractCICModel end

@with_kw struct CIC <: AbstractCICModel
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Symbol
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

@with_kw struct CICStratified <: AbstractCICModelStratified
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Symbol
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
  outcome::Symbol
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
  treatedleft::Int64
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CIC}
end

@with_kw struct CaliperCICStratified <: AbstractCICModelStratified
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Symbol
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
  treatedleft::Dict{Int64, Int64}
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CICStratified}
  labels::Dict{Int64, String}
end

@with_kw struct RefinedCIC <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Symbol
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
  outcome::Symbol
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
  outcome::Symbol
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
  treatedleft::Int64
  estimator::String = "ATT"
  caliper::Dict{Symbol, Float64}
  fullmod = Base.RefValue{CaliperCIC}
end

@with_kw struct RefinedCaliperCICStratified <: AbstractCICModelStratified
  title::String = "refined"
  id::Symbol
  t::Symbol
  outcome::Symbol
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
  treatedleft::Dict{Int64, Int64}
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
  outcome::Symbol
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
  treatedleft::Union{Int64, Dict{Int64, Int64}} = Int64()
  estimator::String
  caliper::Dict{Symbol, Float64} = Dict{Symbol, Float64}()
  labels::Dict{Int64, String} = Dict{Int64, String}()
  refinementnumber::Int = Int()
end
