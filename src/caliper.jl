# caliper.jl

# construction and execution

function caliper(model::CIC, acaliper, dat; dobalance = true)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator = model;

  tobscr = calipermatches(model, acaliper);

  # only include observations that have matches post-caliper
  obsleft = [isassigned(tobscr, i) for i in 1:length(tobscr)];

  modelcal = CaliperCIC(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    timevary = timevary,
    reference = reference,
    F = F, L = L,
    observations = observations[obsleft],
    ids = ids,
    matches = tobscr[obsleft],
    meanbalances = DataFrame(),
    grandbalances = GrandDictNoStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = length(tobscr),
    treatedleft = length(obsleft[obsleft]),
    estimator = estimator,
    caliper = acaliper,
    fullmod = Ref(model)
  );

  if dobalance
    meanbalance!(modelcal, dat);
    grandbalance!(modelcal);
  end

  return modelcal
end

function caliper(model::CICStratified, acaliper, dat; dobalance = true)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator = model;

  @unpack stratifier, strata = model

  tobscr = calipermatches(model, acaliper);

  # only include observations that have matches post-caliper
  obsleft = [isassigned(tobscr, i) for i in 1:length(tobscr)];

  modelcal = CaliperCICStratified(
    title = title,
    id = id,
    t = t,
    outcome = outcome,
    treatment = treatment,
    covariates = covariates,
    timevary = timevary,
    stratifier = stratifier,
    strata = strata[obsleft],
    reference = reference,
    F = F, L = L,
    observations = observations[obsleft],
    ids = ids,
    matches = tobscr[obsleft],
    meanbalances = DataFrame(),
    grandbalances = GrandDictStrat(),
    iterations = iterations,
    results = DataFrame(),
    treatednum = length(tobscr),
    treatedleft = length(obsleft[obsleft]),
    estimator = estimator,
    caliper = acaliper,
    fullmod = Ref(model)
  );

  if dobalance
    meanbalance!(modelcal, dat);
    grandbalance!(modelcal);
  end

  return modelcal
end

# mechanics

function calipervars(caliper, covariates)
  calipers = fill(Inf, length(covariates)+1);
  for (c, covar) in enumerate(vcat(:Mahalanobis, covariates))
    calipers[c] = get(caliper, covar, Inf);
  end
  return calipers
end

function calipermatches(model, acaliper)
  @unpack observations, matches, covariates, F = model;

  calipers = calipervars(acaliper, covariates);

  tobscr = Vector{TobC}(undef, length(observations));
  _caliper!(tobscr, matches, calipers, length(F))

  return tobscr
end

function _caliper!(tobscr, matches, calipers, flen)
  Threads.@threads for i in eachindex(matches)
    tob = @views matches[i]
    @unpack mus, fs, mudistances, ranks = tob;
    tobcr = TobC(
      deepcopy(mus),
      deepcopy(fs),
      Dict{Int, Vector{Int}}()
    ); # probably memory intensive
    calipertobs!(tobcr, mus, fs, mudistances, calipers);
    if !(sum(tobcr.mus) == 0)
      tobscr[i] = tobcr;
      tobcr.mus
      
      for φ in 1:flen
        # tobcr.mus[mus][tobcr.mus[mus]]
        # this should convert the indices
        # tobcr.mus[mus] sends tobcr (calipered) to coords in space of original
        # model true mus. Then, we get only the true ones in that space
        # to get the ranks that are left
        tobcr.ranks[φ] = tob.ranks[φ][tobcr.mus[mus]] # mu[mu] order
        # sorting is preserved, with rejects removed
      end
    end
  end
  return tobscr
end

function calipertobs!(tobcr, mus, fs, mudistances, calipers)
  cnt = 0
  for (m, mu) in enumerate(mus)
    if mu
      fsmu = @views fs[m]
      cnt += 1
      cntfs = 0
      for (φ, fb) in enumerate(fsmu)
        if fb
          cntfs += 1
          dists = mudistances[cnt][cntfs] # an f for an mu
          for c in 1:length(calipers)
            if abs(dists[c]) > calipers[c]
              # change correct fs in tobCR to false
              tobcr.fs[m][φ] = false
              break
            end
          end
        end
      end
      if !any(tobcr.fs[m])
        # if no fs are true, knock the whole match unit out
        tobcr.mus[m] = false
      end
    end
  end
  return tobcr
end
