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
  matchinfos = nothing, obsinfos = nothing,
  savepath = nothing
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

  if !isnothing(savepath)
    save_object(
      savepath * name_model(model) * ".jld2",
      cr
    )
  end
  return cr
end

#### 

function save_record(savepath, modelrecord::ModelRecord)
  save_object(savepath * name_model(modelrecord), modelrecord(modelrecord))
end

function save_records_separate(savepath, models...)
  for model in models
    save_object(
      savepath * name_model(model) * "_" * string(typeof(model)) * ".jld2",
      modelrecord(models[1])
    )
  end
end

function save_records(savepath, models...)
  records = Vector{ModelRecord}(undef, length(models))
  for (i, model) in enumerate(models)
    if i < 5
      records[i] = modelrecord(model)
    end
  end
  if length(models) < 5
    save_object(
      savepath * name_model(models[1]) * ".jld2",
      records
    )
  else
    save_object(
      savepath * name_model(models[1]) * ".jld2",
      [records, models[5:end]...]
    )
  end
end

function modelrecord(model::CIC)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator
  )
end

function modelrecord(model::RefinedCIC)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, refinementnum = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator, refinementnumber = refinementnum
  )
end

function modelrecord(model::CaliperCIC)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, caliper = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results,
    treatednum = treatednum,
    estimator = estimator, caliper = caliper
  )
end

function modelrecord(model::RefinedCaliperCIC)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, caliper, refinementnum = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results,
    treatednum = treatednum,
    estimator = estimator,
    caliper = caliper,
    refinementnumber = refinementnum
  )
end

function modelrecord(model::CICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels
  )
end

function modelrecord(model::RefinedCICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, refinementnum = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator, refinementnumber = refinementnum,
    stratifier = stratifier, strata = strata, labels = labels
  )
end

function modelrecord(model::CaliperCICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, caliper = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results,
    treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    caliper = caliper
  )
end

function modelrecord(model::RefinedCaliperCICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, refinementnum = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    refinementnumber = refinementnum
  )
end
