# refine.jl

# construction and execution

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

# import tscsmethods:_refine,RefinedCaliperCICStratified,@unpack,meanbalance!,grandbalance!,TobR

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
    labels = labels
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

  @unpack treatednum, treatedleft = calmodel;
  
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
    treatedleft = treatedleft,
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
  @unpack treatednum, treatedleft = calmodel;

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
    treatedleft = treatedleft,
    estimator = estimator,
    caliper = caliper,
    fullmod = Ref(calmodel),
    labels = labels
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
