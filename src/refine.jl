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

  @unpack stratifier, strata = model
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

function refinesetup!(tobscr, idlen)
  for i in eachindex(tobscr)
    # positions
    tobscr[i] = TobR(
      mus = fill(false, idlen),
      fs = Vector{Vector{Bool}}(undef, idlen)
    )
  end
  return tobscr
end

function _refine(model, refinementnum)
  @unpack observations, matches, ids, F = model;
  flen = length(F);
  
  tobscr = Vector{TobR}(undef, length(observations));
  refinesetup!(tobscr, length(ids));
  _refine!(tobscr, matches, flen, refinementnum)
  return tobscr
end

function _refine!(tobscr, matches, flen, refinementnum)
  Threads.@threads for i in eachindex(matches)
    @unpack ranks = matches[i];
    @unpack mus, fs = tobscr[i];
    refinetob!(mus, fs, ranks, 1:flen, refinementnum)
  end
  return tobscr
end

function refinetob!(mus, fs, ranks, Φ, refinementnum)
  for φ in Φ
    # the caliper may have left fewer than refinementnum matches
    for r in 1:min(refinementnum, length(ranks[φ]))
      m = ranks[φ][r]
      as = isassigned(fs, m);
      if !as
        fs[m] = fill(false, length(Φ))
      end
      fs[m][φ] = true
      mus[m] = true
    end
  end
end
