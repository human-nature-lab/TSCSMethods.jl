# model.jl

"""
    name_model(model::VeryAbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(model::VeryAbstractCICModel)
  strat = typeof(model) <: AbstractCICModelStratified ? string(model.stratifier) : ""
  return model.title * "_" * string(model.outcome) * "_" * strat
end

# should add treatment
struct CICRecord
  title::String
  type::DataType
  results::DataFrame
  balances::Union{GrandDictStrat, GrandDictNoStrat}
  outcome::Symbol
  covariates::Vector{Symbol}
  ids::Vector{Int}
  observations::Vector{Tuple{Int, Int}}
  F::UnitRange{Int}
  L::UnitRange{Int}
  t::Symbol
  id::Symbol
  stratifier::Symbol
  strata::Vector{Int}
  iterations::Int
  refinementnum::Int
  labels::Dict{Int, String}
end

struct CICRecords
  model::Union{CICRecord, Nothing}
  refinedmodel::Union{CICRecord, Nothing}
  calmodel::Union{CICRecord, Nothing}
  refcalmodel::Union{CICRecord, Nothing}
  matchinfo::Union{DataFrame, Nothing}
  obsinfo::Union{DataFrame, Nothing}
end

function makerecord(m::VeryAbstractCICModel)

  if typeof(m) .∈ Ref(
    [
      CICStratified, RefinedCICStratified,
      CaliperCICStratified, RefinedCaliperCICStratified
    ])
    stratifier = m.stratifier
    strata = m.strata
    labels = m.labels
  else
    stratifier = Symbol()
    strata = Vector{Int}()
    labels = Dict{Int, String}()
  end

  if typeof(m) .∈ Ref(
    [
      RefinedCIC,
      RefinedCICStratified,
      RefinedCaliperCIC,
      RefinedCaliperCICStratified
      ])
    refinementnum = m.refinementnum
  else refinementnum = length(m.ids)
  end
    
  mrecord = CICRecord(
    m.title,
    typeof(m),
    m.results,
    m.grandbalances,
    m.outcome,
    m.covariates,
    m.ids,
    m.observations,
    m.F,
    m.L,
    m.t,
    m.id,
    stratifier,
    strata,
    m.iterations,
    refinementnum,
    labels
  )

  return mrecord
end

function makerecords(dat, savepath, models; obscovars = nothing)

  mr = nothing; mr2 = nothing; mr3 = nothing; mr4 = nothing
  m1 = []; mrf = []
  for model in models
    if typeof(model) <: Union{CIC, CICStratified}
      mr = makerecord(model)
      m1 = model
    end
    if typeof(model) <: Union{RefinedCIC, RefinedCICStratified}
      mr2 = makerecord(model)
    end
    if typeof(model) <: Union{CaliperCIC, CaliperCICStratified}
      mr3 = makerecord(model)
    end
    if typeof(model) <: Union{RefinedCaliperCIC, RefinedCaliperCICStratified}
      mr4 = makerecord(model)
      mrf = model
    end
  end

  # models = [model, refinedmodel, calmodel, refcalmodel]

  covarset = if isnothing(obscovars)
    m1.covariates
  else union(m1.covariates, obscovars)
  end
  
  if !isnothing(mr) & !isnothing(mr4)
    rcinfo = matchinfo(mrf, m1);
    obinfo = obsinfo(
      rcinfo, dat, covarset;
      fullmodobs = m1.observations,
      t = m1.t, id = m1.id
    );
  else
    rcinfo = nothing
    obinfo = nothing
  end
  
  records = CICRecords(
    mr, mr2,
    mr3, mr4,
    rcinfo,
    obinfo
  );

  if !isnothing(savepath)
    save_object(
      savepath * name_model(m1) * ".jld2",
      records
    )
  end

  return records
end

"""
    variable_filter(
      model, variable, dat;
      mn = nothing, mx = nothing
    )

Remove treated observations according to some variable values.
"""
function variable_filter(
  model, variable, dat;
  mn = nothing, mx = nothing
)

  @unpack t, id, treatment = model

  # remove elections prior to March 10

  dt = @subset(dat, $treatment .== 1);

  tple = [(dt[i, t], dt[i, id]) for i in eachindex(dt[!, t])];
  dict = Dict(tple .=> dt[!, variable]);

  obinclude = fill(false, length(model.observations));
  for (i, ob) in enumerate(model.observations)
    cond = true
    if !isnothing(mn)
      cond = cond & (dict[ob] >= mn)
    end
    if !isnothing(mx)
      cond = cond & (dict[ob] <= mx)
    end
      obinclude[i] = cond
  end

  @reset model.observations = model.observations[obinclude];
  @reset model.matches = model.matches[obinclude];
  # @reset model.results = tscsmethods.DataFrame();

  @reset model.treatednum = length(model.observations)
  
  return model
end

"""
    function treatedinfo(
      model, variables, dat;
    )

Gives variable values for treated observations present in the model, for
the chosen set of variables. Order is the same as model.observations.
"""
function treatedinfo(
  model, variables, dat;
)

  @unpack t, id, treatment = model

  # remove elections prior to March 10

  dt = @subset(dat, $treatment .== 1);

  tple = [(dt[i, t], dt[i, id]) for i in eachindex(dt[!, t])];
  
  obout = DataFrame(
    obs = model.observations,    
    )
    
    for variable in variables
      obout[!, variable] = Vector{eltype(dat[!, variable])}(undef, nrow(obout))
      dict = Dict(tple .=> dt[!, variable]);
      
      for (i, e) in enumerate(obout.obs)
        obout[i, variable] = dict[e]
      end
    end
  return obout
end

function relabel!(
  calmodel, refcalmodel, dat; stratifier = nothing, digits = 2
)

  stratifier = if isnothing(stratifier)
    calmodel.stratifier
  else
    stratifier
  end

  # the cal model and refcalmodel have the same stratum ranges
  calinfo = treatedinfo(
    calmodel, [stratifier], dat;
  )
  calinfo[!, :stratum] = calmodel.strata

  calinfo = @chain calinfo begin
    groupby(:stratum)
    combine(stratifier => extrema => stratifier)
  end

  exts = Dict(calinfo.stratum .=> calinfo[!, stratifier])

  relabels = Dict{Int, String}();
  for s in sort(collect(keys(exts)))
    mn, mx = exts[s]
    if typeof(mn) == Float64
      mn = round(mn; digits = digits)
      mx = round(mx; digits = digits)
    end
    relabels[s] = if ismissing(mn) & ismissing(mx)
      "Missing Values"
    else
      if mn < mx
        string(mn) * " to " * string(mx)
      elseif mn == mx
        string(mn)
      end
    end
  end

  @reset refcalmodel.labels = deepcopy(relabels)
  @reset calmodel.labels = deepcopy(relabels)

  # for (k,v) in relabels; calmodel.labels[k] = v end # add labels
  # for (k,v) in relabels; refcalmodel.labels[k] = v end # add labels
  return calmodel, refcalmodel
end

function relabel!(
  m, dat; stratifier = nothing, digits = 2
)

  stratifier = if isnothing(stratifier)
    m.stratifier
  else
    stratifier
  end

  # the cal model and refcalmodel have the same stratum ranges
  calinfo = treatedinfo(
    m, [stratifier], dat;
  )
  calinfo[!, :stratum] = m.strata

  calinfo = @chain calinfo begin
    groupby(:stratum)
    combine(stratifier => extrema => stratifier)
  end

  exts = Dict(calinfo.stratum .=> calinfo[!, stratifier])

  relabels = Dict{Int, String}();
  for s in sort(collect(keys(exts)))
    mn, mx = exts[s]
    if typeof(mn) == Float64
      mn = round(mn; digits = digits)
      mx = round(mx; digits = digits)
    end
    relabels[s] = if ismissing(mn) & ismissing(mx)
      "Missing Values"
    else
      if mn < mx
        string(mn) * " to " * string(mx)
      elseif mn == mx
        string(mn)
      end
    end
  end

  @reset m.labels = deepcopy(relabels)

  # for (k,v) in relabels; calmodel.labels[k] = v end # add labels
  # for (k,v) in relabels; refcalmodel.labels[k] = v end # add labels
  return m
end

#=
function makerecord(
  model::VeryAbstractCICModel,
  matchinfos = nothing,
  obsinfos = nothing,
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
=#