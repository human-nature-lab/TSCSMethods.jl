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

"""
    make_caliper(cc::cicmodel, caliper)

only ever applied to the full model
"""
function make_caliper(cc::cicmodel, caliper)
  # caliper::Dict{Symbol, Float64}()
  inc = applycaliper(cc, caliper);

  cal = calipercicmodel(
    title = cc.title * " caliper",
    id = cc.id,
    t = cc.t,
    outcome = cc.outcome,
    treatment = cc.treatment,
    covariates = cc.covariates,
    timevary = cc.timevary,
    fmin = cc.fmin,
    fmax = cc.fmax,
    stratifier = cc.stratifier,
    reference = cc.reference,
    mmin = cc.mmin,
    mmax = cc.mmax,
    matches = @view(cc.matches[inc, :]), # caliper here
    balances = @view(cc.balances[!, :]),
    # meanbalances = , # do func
    # grandbalances = , # do func
    iterations = cc.iterations,
    results = DataFrame(),
    treatednum = Int64(0),
    treatedleft = Int64(0),
    estimator = "ATT",
    caliper = caliper,
    fullmod = Base.RefValue{cicmodel}
  )

  meanbalance!(cal);

  # meanbalance! does not carry over the stratum label, so reassign
  if cc.stratifier != Symbol("") 
    uset = unique(
      cc.meanbalances[!, [:treattime, :treatunit, :stratum]],
      view = true
    )
    cal.meanbalances = leftjoin(
      cal.meanbalances, uset, on = [:treattime, :treatunit]
    )
  end

  grandbalance!(cal);

  cal.treatednum = cc.treatednum;
  if cc.stratifier == Symbol("")
    cal.treatedleft = nrow(unique(cal.meanbalances[!, [:treattime, :treatunit]]));
  else
    treatedlefts!(cal) # left
  end

  return cal
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

function make_refined(cc::AbstractCICModel; refinementnum = 5)
  # caliper::Dict{Symbol, Float64}()

  # EXTEND THIS FOR STRATIFICATION
  if typeof(cc) == calipercicmodel
    calip = cc.caliper
  else
    calip = Dict{Symbol, Float64}()
  end;

  rf = refinedcicmodel(
    title = cc.title * " refined",
    id = cc.id,
    t = cc.t,
    outcome = cc.outcome,
    treatment = cc.treatment,
    covariates = cc.covariates,
    timevary = cc.timevary,
    fmin = cc.fmin,
    fmax = cc.fmax,
    stratifier = cc.stratifier,
    reference = cc.reference,
    mmin = cc.mmin,
    mmax = cc.mmax,
    matches = refine(cc, refinementnum),
    balances = @view(cc.balances[!, :]),
    # meanbalances = , # do func
    # grandbalances = , # do func
    iterations = cc.iterations,
    results = DataFrame(),
    treatednum = Int64(0),
    treatedleft = Int64(0),
    estimator = "ATT",
    caliper = calip,
    fullmod = Base.RefValue{cicmodel}
  )

  meanbalance!(rf)
  
  # meanbalance! does not carry over the stratum label, so reassign
  if cc.stratifier != Symbol("") 
    uset = unique(
      cc.meanbalances[!, [:treattime, :treatunit, :stratum]],
      view = true
    )
    rf.meanbalances = leftjoin(
      rf.meanbalances, uset, on = [:treattime, :treatunit]
    )
  end
  
  grandbalance!(rf)

  rf.treatednum = cc.treatednum;
  if cc.stratifier == Symbol("")
    rf.treatedleft = nrow(unique(cc.meanbalances[!, [:treattime, :treatunit]]));
  else
    treatedlefts!(rf)
  end

  return rf
end
