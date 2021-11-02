#=
model object
* the data is kept separately
=#

MatchDistances = Vector{Union{Vector{Union{Vector{Union{Float64, Missing}}, Missing}}, Missing}}; # move elsewhere?
MatchUnits = Vector{Vector{Vector{Int}}};

StratDict = Dict{Int64, Dict{Symbol, Float64}};
GrandDictNoStrat = Dict{Symbol, Union{Vector{Float64}, Float64}};
GrandDictStrat = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}};

abstract type AbstractCICModel end

@with_kw struct CIC <: AbstractCICModel
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  reference::Int = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matchdistances::Vector{MatchDistances} = Vector{MatchDistances}()
  matchunits::MatchUnits = MatchUnits()
  balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::Union{GrandDictNoStrat, GrandDictStrat} = GrandDictNoStrat()
  iterations::Int64 = 500
  results::DataFrame = DataFrame()
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  estimator::String = "ATT"
end

@with_kw struct CaliperCIC <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  reference::Int64 = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matchunits::MatchUnits
  balances::SubDataFrame # = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::Union{GrandDictNoStrat, GrandDictStrat} = GrandDictNoStrat()
  iterations::Int64 = 500
  results::DataFrame = DataFrame()
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  treatedleft::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  estimator::String = "ATT"
  caliper::Union{Dict{Symbol, Float64}, StratDict} = Dict{Symbol, Float64}()
  fullmod = Base.RefValue{cicmodel}
end


@with_kw struct RefinedCIC <: AbstractCICModel
  title::String = "refined"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  stratifier::Symbol = Symbol()
  reference::Int64 = -1
  F::UnitRange{Int}
  L::UnitRange{Int}
  observations::Vector{Tuple{Int, Int}}
  ids::Vector{Int}
  matchunits::MatchUnits
  balances::SubDataFrame
  meanbalances::DataFrame = DataFrame()
  grandbalances::Union{GrandDictNoStrat, GrandDictStrat} = GrandDictNoStrat()
  iterations::Int64 = 500
  results::DataFrame = DataFrame()
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  treatedleft::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  estimator::String = "ATT"
  refinementnum::Int64 = 5
  caliper::Union{Dict{Symbol, Float64}, StratDict} = Dict{Symbol, Float64}()
  fullmod = Base.RefValue{cicmodel}
end
