# new functions

import tscsmethods:mean,weuclidean,mahalanobis,samplecovar,cicmodel,match

treatga!(dat);

model = deathmodel("ga test", :gaspecial, :epi);

dat = dataprep(dat, model);

@time match!(model, dat);

# balancing
@time fullbalance!(model, dat);

@time meanbalance!(model);

#time grandbalance!(model::AbstractCICModel)

# typeof MatchUnits
fslen = 31
@time refined = refine(matchunits, refinementnum, fslen);

cals = [Inf, 0.15, 0.025, 0.0001]

# typeof MatchUnits
calipers = Dict{Symbol, Float64}()
[calipers[covar] = 0.5 for covar in model.covariates];
@time calipered =  caliper(matchunits, matchdistances, ids, model, calipers);
@time refinedcalipered = refine(calipered, refinementnum, fslen);

cals = [Inf, 0.5, 0.5, 0.5]

##

# meanbalance
"""
simple implementation from

https://discourse.julialang.org/t/how-to-calculate-a-weighted-mean-with-missing-observations/19281
"""
function weightedmean(x, w)
  wxsum = wsum = 0
  for (x,w) in zip(x,w)
      wx = w*x
      if !ismissing(wx)
          wxsum += wx
          wsum += w
      end
  end
  return wxsum / wsum
end

# another way, crate a 31xmus matrix, populate accordingly with distances, and sort
# then construct model, rather than copy
# (this may be much better actually)


function refine(model::CaliperCIC; refinementnum = 5)
  @unpack title, id, t, outcome, treatment, covariates, timevary, reference, F, L, observations, ids, iterations, estimator = model;
  
  @unpack matches = model;
  tobscr = Vector{TobCR}(undef, length(matches));
  _refine!(tobscr, matches, F, refinementnum);

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
    treatednum = model.treatednum,
    treatedleft = model.treatedleft,
    estimator = estimator,
    caliper = caliper,
    fullmod = Ref(model)
  );

  meanbalance!(modelref, dat);
  grandbalance!(modelref);

  return modelcalref
end

###
vn = VariableNames();
acaliper = Dict(
  # :Mahalanobis => 100,
  vn.pd => 0.5,
  vn.cdr => 0.5,
  vn.fc => 0.5
);

@time calmod = caliper(model, acaliper, dat);

"""
Since model ids are stored in the same order always, create a map from an id to an index. Used mainly for caliper.
"""
function getidmap(uid)
  idmap = Dict{Int, Int}();
  sizehint!(idmap, length(uid))
  [idmap[u] = ux for (ux, u) in enumerate(uid)];
  return idmap
end

function refine(matchunits, refinementnum, fslen)
  refined = similar(matchunits);
  _dorefine!(refined, matchunits, refinementnum, fslen);
  return refined
end

function _dorefine!(refined, matchunits, refinementnum, fslen)
  for i in eachindex(matchunits)
    refined[i] = Vector{Vector{Int}}(undef, fslen);
    for φ in 1:fslen
      refined[i][φ] = if length(matchunits[i][φ]) >= refinementnum
         matchunits[i][φ][1:refinementnum]
      else matchunits[i][φ][1:end]
      end
    end
  end
  return refined
end
