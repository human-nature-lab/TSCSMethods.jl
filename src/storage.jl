"""
    name_model(model::VeryAbstractCICModel)

Generate the filename for a set of models.
"""
function name_model(model::VeryAbstractCICModel)
  strat = typeof(model) <: AbstractCICModelStratified ? string(model.stratifier) : ""

  nms = if typeof(model.outcome) == Symbol
    string(model.outcome)
  else
    "multiple_" * string(model.outcome[1])
  end

  return model.title * "_" * nms * "_" * strat
end

# should add treatment
struct CICRecord
  title::String
  type::DataType
  results::DataFrame
  balances::Union{GrandDictStrat, GrandDictNoStrat}
  outcome::Union{Symbol, Vector{Symbol}}
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