# refine.jl

# construction and execution

"""
    refine(
        model::CIC, 
        dat::DataFrame;
        refinementnum::Int = 5, 
        dobalance::Bool = true,
        doestimate::Bool = true
    ) -> RefinedCIC

Refine a CIC model by keeping only the best `refinementnum` matches for each treated unit.

# Arguments
- `model::CIC`: A matched CIC model
- `dat::DataFrame`: Input data containing all model variables
- `refinementnum::Int`: Number of best matches to retain per treated unit (default: 5)
- `dobalance::Bool`: Whether to perform balancing on the refined model (default: true)
- `doestimate::Bool`: Whether to perform estimation on the refined model (default: true)

# Returns
- `RefinedCIC`: A refined version of the input model with reduced matches

# Description
Refinement improves match quality by keeping only the closest matches for each treated unit,
based on the distance metrics computed during the matching stage. This typically improves
covariate balance and estimation precision at the cost of statistical power.

# Examples
```julia
# After matching a model
model = makemodel(data, :time, :id, :treatment, :outcome, covariates, timevary, F, L)
match!(model, data)

# Refine to top 3 matches per treated unit
refined_model = refine(model, data; refinementnum = 3)

# Refine without automatic balancing and estimation
refined_model = refine(model, data; refinementnum = 5, dobalance = false, doestimate = false)
```

# See Also
- [`caliper`](@ref): Alternative approach using distance thresholds
- [`match!`](@ref): Initial matching procedure
- [`balance!`](@ref): Balance assessment
"""
function refine(
  model::CIC, dat;
  refinementnum = 5, dobalance = true,
  doestimate = true)
  
  tobscr = _refine(model, refinementnum);

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator = model;

  modelref = RefinedCIC(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    refinementnum = refinementnum,
    timevary = timevary,
    reference = reference,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobscr,
    meanbalances = DataFrame(),
    grandbalances = GrandDictNoStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = length(tobscr),
    estimator = estimator,
    fullmod = Ref(model)
  );

  if dobalance
    meanbalance!(modelref, dat);
    grandbalance!(modelref);
  end

  if doestimate
    estimate!(modelref, dat)
  end

  return modelref
end

# import TSCSMethods:_refine,RefinedCaliperCICStratified,@unpack,meanbalance!,grandbalance!,TobR

function refine(
  model::CICStratified, dat;
  refinementnum = 5, dobalance = true, doestimate = true
)
  
  tobscr = _refine(model, refinementnum);

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator, labels = model;

  @unpack stratifier, strata = model # no. treated obs don't change

  modelref = RefinedCICStratified(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    refinementnum = refinementnum,
    timevary = timevary,
    stratifier = stratifier,
    strata = strata,
    reference = reference,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobscr,
    meanbalances = DataFrame(),
    grandbalances = GrandDictStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = length(tobscr),
    estimator = estimator,
    fullmod = Ref(model),
    labels = deepcopy(labels)
  );

  if dobalance
    meanbalance!(modelref, dat);
    grandbalance!(modelref);
  end

  if doestimate
    estimate!(modelref, dat)
  end

  return modelref
end

function refine(
  calmodel::CaliperCIC, dat;
  refinementnum = 5, dobalance = true, doestimate = true
)
  
  tobscr = _refine(calmodel, refinementnum)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator, caliper = calmodel;

  @unpack treatednum = calmodel;
  
  modelcalref = RefinedCaliperCIC(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    refinementnum = refinementnum,
    timevary = timevary,
    reference = reference,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobscr,
    meanbalances = DataFrame(),
    grandbalances = GrandDictNoStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = treatednum,
    estimator = estimator,
    caliper = caliper,
    fullmod = Ref(calmodel)
  );

  if dobalance
    meanbalance!(modelcalref, dat);
    grandbalance!(modelcalref);
  end

  if doestimate
    estimate!(modelcalref, dat)
  end

  return modelcalref
end

function refine(
  calmodel::CaliperCICStratified, dat;
  refinementnum = 5, dobalance = true, doestimate = true)
  
  tobscr = _refine(calmodel, refinementnum)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator, labels, caliper = calmodel;

  @unpack stratifier, strata = calmodel
  @unpack treatednum = calmodel;

  modelcalref = RefinedCaliperCICStratified(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    refinementnum = refinementnum,
    timevary = timevary,
    stratifier = stratifier,
    strata = strata,
    reference = reference,
    F = F, L = L,
    observations = observations,
    ids = ids,
    matches = tobscr,
    meanbalances = DataFrame(),
    grandbalances = GrandDictStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = treatednum,
    estimator = estimator,
    caliper = caliper,
    fullmod = Ref(calmodel),
    labels = deepcopy(labels)
  );

  if dobalance
    meanbalance!(modelcalref, dat);
    grandbalance!(modelcalref);
  end

  if doestimate
    estimate!(modelcalref, dat)
  end

  return modelcalref
end

# mechanics

function _refine(model, refinementnum)
  @unpack observations, matches, ids, F = model;
  
  # new objects
  tobscr = Vector{TobR}(undef, length(observations));
  _refine_assign!(tobscr, matches, refinementnum, length(ids), length(F));

  return tobscr
end

function _refine_assign!(tobscr, matches, refinementnum, idlen, flen)
  Threads.@threads for i in eachindex(tobscr)
    tobscr[i] = TobR(
      mus = fill(false, idlen, flen)
    )
    for φ in 1:flen
      rlen = length(matches[i].ranks[φ])
      if rlen > 0
        for n in 1:min(refinementnum, rlen)
          idx = matches[i].ranks[φ][n]
          tobscr[i].mus[idx, φ] = true
        end
      end
    end
  end
  return tobscr
end
