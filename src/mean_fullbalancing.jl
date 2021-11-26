# mean_fullbalancing.jl

## scratch

# @time balances = fullbalance(model, dat);

# names(model.meanbalances)

# size(model.meanbalances)
# length(model.meanbalances[1, vn.cdr])
# length(model.meanbalances[1, vn.cdr][1])
# length(model.meanbalances[1, vn.pd])

# MBCOPY = deepcopy(model.meanbalances);

# mean_fullbalance!(model, balances);

##

"""
    fidx(f::Int, mlen::Int, fmin::Int)

get indices from f
"""
fidx(f::Int, mlen::Int, fmin::Int) = (f - fmin + 1):(f - fmin + mlen);

"""
    mean_fullbalance!(model)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function mean_fullbalance!(model::AbstractCICModel, balances::DataFrame)

  @unpack matches, ids = refinedmodel;
  @unpack covariates, timevary = refinedmodel;
  @unpack F, L = model;

  fmin = minimum(F); mmlen = length(L);

  allocate_meanbalances!(refinedmodel);

  _meanbalance!(
    eachrow(model.meanbalances), eachrow(balances),
    matches, ids, covariates, timevary, mmlen, fmin
  );

  return model
end

### meanbalances loop
function _meanbalance!(
  barrows, balrows,
  matches, ids, covariates, timevary, mmlen, fmin
)
  
  # for i in eachindex(balrows)
  Threads.@threads for i in eachindex(balrows)
    br = balrows[i];
    mr = barrows[i];

    mus, efsets = matchassignments(matches[i], ids);

    rdx = refinedmodel.matches[1].mus[model.matches[1].mus] # 2096

    _covar_meanbalance!(
      mr, br, efsets, covariates, timevary, mmlen, fmin, rdx,
    );
  end
  return barrows
end
  
### inner functions
function _covar_meanbalance!(
  mr, br, efsets, covariates, timevary, mmlen, fmin, rdx
)
  for covar in covariates
    # this is a [row, covar] for MB:
    #   the means across fs for a covariate (for a treated observation)
    if timevary[covar]
      row_covar_meanbalance!(
        mr[covar], br[covar][rdx], mr[:fs], efsets, mmlen, fmin
      );
    else
      row_covar_meanbalance!(mr[covar], br[covar], mr[:fs], efsets);
    end
  end
end

function row_covar_meanbalance!(
  Holding::Vector{Vector{Union{Missing, Float64}}},
  br_covar,
  fs, efsets, mmlen, fmin
)

  # cntf = 0
  for (φ, f) in enumerate(fs)
    # cntf += 1
    if f # if there are any valid matches for that f
      ls = fidx(φ + fmin - 1, mmlen, fmin);
      # holding = Vector{Union{Missing, Vector{Union{Float64, Missing}}}}(
      #   missing, length(br_covar)
      # );
      efmφ  = [efsets[m][φ] for m in 1:length(br_covar)]
      holding = Vector{Vector{Union{Float64, Missing}}}(
        undef, sum(efmφ)
      );
      mcnt = 0
      for m in eachindex(br_covar) # this is the site of averaging
        # for a single match
        if efsets[m][φ]
          mcnt += 1
          # if that match is valid for that f
          # else leave as missing
          holding[mcnt] = br_covar[m][ls];
        end
      end
    end
    
    # already missing
    # Holding[φ] = Vector{Union{Missing, Float64}}(undef, mmlen);
    __row_covar_meanbalance!(Holding, holding, φ, ls);
  end
  return Holding
end

function __row_covar_meanbalance!(Holding, holding, φ, ls)
  for l in eachindex(ls)
    lvec = Vector{Union{Float64, Missing}}(missing, length(holding));
    for (v, hold) in enumerate(holding)
      lvec[v] = hold[l]
    end
    # [isdefined(holding, lx) for lx in eachindex(holding)]
    Holding[φ][l] = !isempty(skipmissing(lvec)) ? mean(skipmissing(lvec)) : missing
  end
  return Holding
end

function row_covar_meanbalance!(
  Holding, br_covar, fs, efsets,
)
  for (φ, f) in enumerate(fs)
    if f
      holding = Vector{Union{Float64, Missing}}(undef, length(br_covar));
      for m in eachindex(br_covar) # this is the site of averaging
        if efsets[m][φ]
          holding[m] = br_covar[m];
        end
      end
      Holding[φ] = mean(skipmissing(holding))
    end
  end
  return Holding
end

function refinebalances(refinedmodel, model, balances)
  @unpack covariates, matches, observations = model;
  bals = copy(balances);

  # remove tobs that do not exist in ref / cal model
  if length(observations) != length(refinedmodel.observations)
    keep = observations .∈ Ref(calmodel.observations);
    bals = bals[keep, :];
  end

  for covar in covariates
    for i in eachindex(observations)
      rdx = refinedmodel.matches[i].mus[matches[i].mus] # 2096
      bals[i, covar] = balances[i, covar][rdx]
    end
  end
  return bals
end

function balances!(refcalmodel, model, balances)
  bals = refinebalances(refcalmodel, model, balances);
  mean_fullbalance!(refcalmodel, bals);
  grandbalance!(refcalmodel);
  return refcalmodel
end