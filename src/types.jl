#=
model object
* the data is kept separately
=#

StratDict = Dict{Int64, Dict{Symbol, Float64}};
GrandDictNoStrat = Dict{Symbol, Union{Vector{Float64}, Float64}};
GrandDictStrat = Dict{Int64, Dict{Symbol, Union{Float64, Vector{Float64}}}};

abstract type AbstractCICModel end

@with_kw mutable struct cicmodel <: AbstractCICModel
  title::String = ""
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  fmin::Int64
  fmax::Int64
  stratifier::Symbol = Symbol()
  reference::Int64 = -1
  mmin::Int64
  mmax::Int64
  matches::DataFrame = DataFrame()
  balances::DataFrame = DataFrame()
  meanbalances::DataFrame = DataFrame()
  grandbalances::Union{GrandDictNoStrat, GrandDictStrat} = GrandDictNoStrat()
  iterations::Int64 = 500
  results::DataFrame = DataFrame()
  treatednum::Union{Int64, Dict{Int64, Int64}} = Int64(0)
  estimator::String = "ATT"
end

@with_kw mutable struct calipercicmodel <: AbstractCICModel
  title::String = "caliper"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  fmin::Int64
  fmax::Int64
  stratifier::Symbol = Symbol()
  reference::Int64 = -1
  mmin::Int64
  mmax::Int64
  matches::SubDataFrame # = DataFrame()
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


@with_kw mutable struct refinedcicmodel <: AbstractCICModel
  title::String = "refined"
  id::Symbol
  t::Symbol
  outcome::Symbol
  treatment::Symbol
  covariates::Vector{Symbol}
  timevary::Dict{Symbol, Bool}
  fmin::Int64
  fmax::Int64
  stratifier::Symbol = Symbol()
  reference::Int64 = -1
  mmin::Int64
  mmax::Int64
  matches::SubDataFrame
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
