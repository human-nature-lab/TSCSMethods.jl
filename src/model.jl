# model.jl

"""
    name_model(cc::AbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(model::AbstractCICModel)
  strat = typeof(model) <: AbstractCICModelStratified ? string(model.stratifier) : ""
  return model.title * "_" * string(model.outcome) * "_" * strat
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
    results = results, treatednum = treatednum,
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
    results = results, treatednum = treatednum,
    estimator = estimator,
    caliper = caliper, refinementnumber = refinementnum
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
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    refinementnumber = refinementnum
  )
end

function save_record(savepath, modelrecord::ModelRecord)
  JLD2.save_object(savepath * name_model(modelrecord), modelrecord(modelrecord))
end

function save_records_separate(savepath, models...)
  for model in models
    JLD2.save_object(
      savepath * name_model(model) * "_" * string(typeof(model)) * ".jld2",
      modelrecord(model)
    )
  end
end

function save_records(savepath, models...)
  records = Vector{ModelRecord}(undef, length(models))
  for (i, model) in enumerate(models)
    records[i] = modelrecord(model)
  end
  JLD2.save_object(
    savepath * name_model(model) * ".jld2",
    records
  )
end