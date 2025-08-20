# caliper.jl

# construction and execution

"""
    caliper(
        model::CIC, 
        acaliper::Dict{Symbol, Float64}, 
        dat::DataFrame; 
        dobalance::Bool = true
    ) -> CaliperCIC

Apply caliper restrictions to a CIC model by excluding matches beyond specified distance thresholds.

# Arguments
- `model::CIC`: A matched CIC model
- `acaliper::Dict{Symbol, Float64}`: Dictionary mapping covariates to maximum allowed distances
- `dat::DataFrame`: Input data containing all model variables
- `dobalance::Bool`: Whether to perform balancing on the calipered model (default: true)

# Returns
- `CaliperCIC`: A calipered version of the input model with restricted matches

# Description
Caliper restrictions improve match quality by excluding matches where the distance
on specified covariates exceeds the given thresholds. This helps ensure that
matched units are similar on the most important dimensions, potentially improving
balance and reducing bias.

# Examples
```julia
# After matching a model
model = makemodel(data, :time, :id, :treatment, :outcome, covariates, timevary, F, L)
match!(model, data)

# Apply calipers - exclude matches with distance > 0.25 on covar1 or > 0.5 on covar2
caliper_specs = Dict(:covar1 => 0.25, :covar2 => 0.5)
calipered_model = caliper(model, caliper_specs, data)

# Apply calipers without automatic balancing
calipered_model = caliper(model, caliper_specs, data; dobalance = false)
```

# Notes
- Stricter calipers improve match quality but reduce the number of matches
- Use [`autobalance`](@ref) to automatically determine appropriate caliper values

# See Also
- [`refine`](@ref): Alternative approach using match ranking
- [`autobalance`](@ref): Automatic caliper selection
- [`match!`](@ref): Initial matching procedure
"""
function caliper(model::CIC, acaliper, dat; dobalance = true)

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator = model;

  tobscr = calipermatches(model, acaliper);

  obsleft = Vector{Bool}(undef, length(tobscr));
  get_observations_left!(obsleft, tobscr);

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

  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator, labels = model;

  @unpack stratifier, strata = model;

  tobscr = calipermatches(model, acaliper);

  obsleft = Vector{Bool}(undef, length(tobscr));
  get_observations_left!(obsleft, tobscr);

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
    treatednum = model.treatednum,
    estimator = estimator,
    caliper = acaliper,
    fullmod = Ref(model),
    labels = deepcopy(labels)
  );

  if dobalance
    meanbalance!(modelcal, dat);
    grandbalance!(modelcal);
  end

  return modelcal
end

function get_observations_left!(obsleft, tobscr)
  # only include observations that have matches post-caliper
  for (i, e) in enumerate(tobscr)
    obsleft[i] = sum(e.mus) > 0
  end
  return obsleft
end

# mechanics

function calipermatches(model, acaliper)
  @unpack observations, matches, covariates = model;

  calipers = calipervars(acaliper, covariates);

  tobscr = Vector{TobC}(undef, length(observations));
  _fillcal!(tobscr, matches);

  _caliper!(tobscr, matches, calipers);

  return tobscr
end

function calipervars(caliper, covariates)
  calipers = fill(Inf, length(covariates)+1);
  for (c, covar) in enumerate(vcat(:Mahalanobis, covariates))
    calipers[c] = get(caliper, covar, Inf);
  end
  return calipers
end

function _fillcal!(tobscr, matches)
  Threads.@threads for i in eachindex(matches)
    tob = @views matches[i]
    @unpack mus, distances, ranks = tob;

    tobscr[i] = TobC(
      deepcopy(mus),
      deepcopy(ranks)
    );
  end
  return tobscr
end

function _caliper!(tobscr, matches, calipers)
  
  Threads.@threads for m in eachindex(tobscr)

    tob = @views tobscr[m];
    @unpack distances = @views matches[m];

    @unpack mus, ranks = tob;

    anymus = Vector{Bool}(undef, size(mus)[1]);
    get_anymus!(anymus, mus);
    validmus = @views mus[anymus, :];

    _inner_caliper!(validmus, permutedims(distances), calipers);
    rerank!(ranks, mus)

  end
end

function _inner_caliper!(validmus, distances, calipers)
  for i in 1:size(distances[1])[1]
    for j in 1:size(distances[1])[2]
      if validmus[j, i]
        for k in 1:length(calipers)
          if abs(distances[k][i, j]) > calipers[k]
            validmus[j, i] = false
          end
        end
      end
    end
  end
  return validmus
end

function _inner_caliper!(validmus, distances_tr::Matrix{Float64}, calipers)
  # whole rows get blasted, since the covariate window is fixed
  for (i, r) in enumerate(eachrow(validmus))
    for j in eachindex(calipers)
      if abs(distances_tr[j, i]) > calipers[j]
        r .= false
      end
    end
  end
  return validmus
end

function rerank!(ranks, mus)
  for key in keys(ranks)
    ranks[key] = ranks[key][mus[:, key][ranks[key]]]
  end
  return ranks
end
