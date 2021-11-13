# model.jl

"""
    name_model(cc::AbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(cc::AbstractCICModel)
  strat = Symbol("") == cc.stratifier ? "" : string(cc.stratifier)
  return cc.title * "_" * string(cc.outcome) * "_" * strat
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

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, refinementnumber = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator, refinementnumber = refinementnumber
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

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, caliper, refinementnumber = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    caliper = caliper, refinementnumber = refinementnumber
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

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, refinementnumber = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    refinementnumber = refinementnumber
  )
end

function save_record(model::CaliperCICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, caliper = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    caliper = caliper
  )
end

function save_record(model::CaliperCICStratified)

  @unpack title, id, t, outcome, treatment, covariates, reference, F, L, observations, ids, grandbalances, iterations, results, treatednum, estimator, stratifier, strata, labels, refinementnumber, caliper = model

  return ModelRecord(
    title = title, id = id, t = t, outcome = outcome,
    treatment = treatment,
    covariates = covariates, reference = reference,
    F = F, L = L, observations = observations, ids = ids,
    grandbalances = grandbalances, iterations = iterations,
    results = results, treatednum = treatednum,
    estimator = estimator,
    stratifier = stratifier, strata = strata, labels = labels,
    caliper = caliper, refinementnumber = refinementnumber
  )
end
