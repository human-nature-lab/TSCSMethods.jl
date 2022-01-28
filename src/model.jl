# model.jl

"""
    name_model(cc::AbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(model::VeryAbstractCICModel)
  strat = typeof(model) <: AbstractCICModelStratified ? string(model.stratifier) : ""
  return model.title * "_" * string(model.outcome) * "_" * strat
end


@with_kw struct CICTrace
  title::String
  observationsinfo::Dict{String, DataFrame}
  matchinginfo = Dict{String, DataFrame}
  results::Dict{String, DataFrame}
  balances::Dict{String, Union{GrandDictStrat, GrandDictNoStrat}}
  outcome::Symbol
  covariates::Vector{Symbol}
  ids::Vector{Int}
  t::Symbol
  id::Symbol
  stratifier::Symbol = Symbol()
  iterations::Int
  refinementnum::Int
end

function makerecord(
  model::Union{AbstractCICModel, AbstractCICModelStratified},
  models;
  matchinfos = nothing, obsinfos = nothing
)
  
  typedict = Dict(
    CIC => "all",
    CICStratified => "all",
    CaliperCIC => "caliper",
    CaliperCICStratified => "caliper",
    RefinedCIC => "refined",
    RefinedCICStratified => "refined",
    RefinedCaliperCIC => "refined caliper",
    RefinedCaliperCICStratified => "refined caliper"
  )

  refdict = Dict(
    RefinedCIC => true,
    RefinedCICStratified => true,
    RefinedCaliperCIC => true,
    RefinedCaliperCICStratified => true
  )

  rfd(m, refdict) = get(refdict, typeof(m), false)

  ftd(m) = typedict[typeof(m)];

  if typeof(model) == AbstractCICModelStratified
    stratvar = model.stratifier
  else stratvar = Symbol()
  end

  refnum = 0
  for m in enumerate(models)
    if rfd(m, refdict)
      refnum = m.refinementnum
    end
  end

  if refnum == 0
    refnum = length(model.ids)
  end

  cr = CICTrace(
    title = model.title,
    observationsinfo = Dict{String, DataFrame}(),
    matchinginfo = Dict{String, DataFrame}(),
    results = Dict{String, DataFrame}(),
    balances = Dict{String, Union{GrandDictStrat, GrandDictNoStrat}}(),
    outcome = model.outcome,
    covariates = model.covariates,
    ids = model.ids,
    t = model.t,
    id = model.id,
    stratifier = stratvar,
    iterations = model.iterations,
    refinementnum = refnum
  )
  
  cr.results[ftd(model)] = model.results; 

  for m in models
    cr.results[ftd(m)] = m.results; 
  end

  if !isnothing(obsinfos)
    for obi in obsinfos
      cr.observationsinfo[ftd(obi[1])] = obi[2]
    end
  end
  
  if !isnothing(matchinfos)
    for mtc in matchinfos
      cr.matchinginfo[ftd(mtc[1])] = mtc[2]
    end
  end

  return cr
end
