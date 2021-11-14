# autobalancing.jl

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

"""
    autobalance(
      model;
      refinementnum = 5,
      calmin = 0.08, step = 0.05, initial_bals = false
    )

Automatically balance via a simple algorithm. Start with initial caliper of 1.0, and subtract `step` whenever the grand mean balance threshold (0.1) is not met.
"""
function autobalance(
  model;
  threshold = 0.1,
  min_treated_obs = 10,
  refinementnum = 5, calmin = 0.1, step = 0.05, initial_bals = false
)

if !initial_bals
  acaliper = Dict{Symbol, Float64}();
  for c in model.covariates
    acaliper[c] = 1.0
  end
end

  calmodel = caliper(model, acaliper, dat; dobalance = false);
  refcalmodel = refine(
    calmodel; refinementnum = refinementnum, dobalance = true
  );

  # check calr only
  bc = checkbalances(refcalmodel; threshold = threshold)

  while any(values(bc)) & (refcalmodel.treatedleft >= min_treated_obs) & all(values(acaliper) .>= calmin)

    covset = [k for k in keys(bc)];
    for covar in sample(covset, length(covset); replace = false)
      if bc[covar] & (acaliper[covar] > calmin)
        acaliper[covar] = acaliper[covar] - step
      end
    end
    calmodel = acaliper(model, acaliper, dat; dobalance = false);
    refcalmodel = refine(
      calmodel; refinementnum = refinementnum, dobalance = true
    );
    bc = checkbalances(refcalmodel; threshold = threshold)
  end

  return calmodel, refcalmodel
end