# autobalancing.jl

"""
    autobalance(
      model, dat;
      refinementnum = 5,
      calmin = 0.08, step = 0.05,
      initial_bals = nothing
    )

Automatically balance via a simple algorithm. Start with initial caliper of 1.0, and subtract `step` whenever the grand mean balance threshold (0.1) is not met.

initial_bals is specified, work downward from the initial specified caliper for one or more selected variables. Unspecified variables start at a caliper of 1.
"""
function autobalance(
  model, dat;
  threshold = 0.1,
  min_treated_obs = 10,
  refinementnum = 5, calmin = 0.1, step = 0.05,
  initial_bals = nothing,
  doestimate = true,
  verbose = true,
  dooverall = false,
  bayesfactor = false
)

  @unpack ids = model;
  @unpack t, id, treatment, covariates = model;
  @unpack F, L = model;

  fmax = maximum(F);
  mmin = minimum(L);
  
  # import TSCSMethods:make_groupindices,caliper,refine,calipermatches,calipervars,TobC,_fillcal!,_caliper!,get_anymus!,_inner_caliper!

  tg, rg, _ = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmax, mmin,
    Matrix(dat[!, covariates]);
  );

  if isnothing(initial_bals)
    acaliper = Dict{Symbol, Float64}();
    for c in covariates
      acaliper[c] = 1.0
    end

    if verbose
      println(acaliper)
    end
  else
    acaliper = Dict{Symbol, Float64}();
    for c in covariates
      acaliper[c] = 1.0
    end
    for (k, v) in initial_bals
      acaliper[k] = v
    end
  end

  calmodel = caliper(model, acaliper, dat; dobalance = false);
  refcalmodel = refine(
    calmodel, dat;
    refinementnum = refinementnum, dobalance = false, doestimate = false
  );
  meanbalance!(refcalmodel, dat, tg, rg);
  grandbalance!(refcalmodel)

  bc = checkbalances(refcalmodel; threshold = threshold);

  while checkwhile(refcalmodel, acaliper, min_treated_obs, calmin, bc)

    covset = [k for k in keys(bc)];
    for covar in sample(covset, length(covset); replace = false)
      if bc[covar] & (acaliper[covar] > calmin)
        acaliper[covar] = acaliper[covar] - step
      end
    end
    calmodel = caliper(model, acaliper, dat; dobalance = false);
    refcalmodel = refine(
      calmodel, dat;
      refinementnum = refinementnum, dobalance = false, doestimate = false
    );
    
    meanbalance!(refcalmodel, dat, tg, rg)
    grandbalance!(refcalmodel)

    bc = checkbalances(refcalmodel; threshold = threshold)

    if verbose
      println(refcalmodel.caliper)
    end
  end

  meanbalance!(calmodel, dat, tg, rg)
  grandbalance!(calmodel)

  if doestimate & !dooverall
    estimate!(refcalmodel, dat; bayesfactor = bayesfactor)
    estimate!(calmodel, dat; bayesfactor = bayesfactor)
  elseif doestimate & dooverall
      overall = estimate!(refcalmodel, dat; overall = true, bayesfactor = bayesfactor)
      estimate!(calmodel, dat; bayesfactor = bayesfactor)
  end

  if !dooverall
    return calmodel, refcalmodel
  else
    return calmodel, refcalmodel, overall
  end
end

"""
    checkbalances(
      m::Dict{Symbol, Union{Float64, Vector{Float64}}};
      threshold = 0.1, stratareduce = true
    )

Simply check whether the grand means are above the std. balance threshold. Returns a Bool for each covariate.
If Stratareduce is true, then the strata balances will be agggregated to the covariate level, such that a violation in any caliper triggers a violation in the aggregated output.
"""
function checkbalances(
  model::AbstractCICModel;
  threshold = 0.1,
)

  if isempty(model.grandbalances)
    error("no grandbalances present")
  end

  chk = Dict{Symbol, Bool}();
  for (k, v) in model.grandbalances
    _balancecheck!(chk, v, k, threshold)
  end

  return chk
end

function checkbalances(
    model::AbstractCICModelStratified;
  threshold = 0.1,
  stratareduce = true
)

  if model.stratifier == Symbol("")
    chk = Dict{Symbol, Bool}();
    for (k, v) in model.grandbalances
      _balancecheck!(chk, v, k, threshold)
    end
  else
    chk = Dict{Int, Dict{Symbol, Bool}}();
    [chk[s] = Dict{Symbol, Bool}() for s in keys(model.grandbalances)]
    for (k, v) in model.grandbalances
      for (ki, vi) in v
      _balancecheck!(chk[k], vi, ki, threshold)
      end
    end
  end

  if stratareduce & (model.stratifier != Symbol(""))
    bc2 = Dict{Symbol, Bool}()
    for covar in model.covariates
      bc2[covar] = false
    end
    
    for v in values(chk)
      for (covar, booléen) in v
        if booléen
          bc2[covar] = true
        end
      end
    end
    return bc2
  end

  return chk
end

function _balancecheck!(chki, v, k, threshold)
  if any(abs.(v) .> threshold)
    chki[k] = true
  else
    chki[k] = false
  end
  return chki
end

function checkwhile(
  refcalmodel::RefinedCaliperCIC, acaliper, min_treated_obs, calmin, bc
)

  treated, _ = unitcounts(refcalmodel);

  # the minimum no. of treated observations left over all Fs, strata
  mintreated = minimum(treated);

  return any(values(bc)) & (mintreated >= min_treated_obs) & all(values(acaliper) .>= calmin)
end

function checkwhile(
  refcalmodel::RefinedCaliperCICStratified, acaliper, min_treated_obs, calmin, bc
)

  treated, _ = unitcounts(refcalmodel);

  # the minimum no. of treated observations left over all Fs, strata
  mintreated = minimum([minimum(treated[i]) for i in keys(treated)]);

  return any(values(bc)) & (mintreated >= min_treated_obs) & all(values(acaliper) .>= calmin)
end
