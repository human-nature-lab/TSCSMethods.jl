# mean_fullbalancing.jl

## scratch

@time balances = fullbalance(model, dat);

names(model.meanbalances)

size(model.meanbalances)
length(model.meanbalances[1, vn.cdr])
length(model.meanbalances[1, vn.cdr][1])
length(model.meanbalances[1, vn.pd])

MBCOPY = deepcopy(model.meanbalances);

mean_fullbalance!(model, balances);

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

  @unpack observations, matches, ids = model;
  @unpack covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L, reference = model;

  fmin = minimum(F); mmlen = length(L);

  model.meanbalances = setup_meanbalance(model);

  _meanbalance!(
    eachrow(model.meanbalances), eachrow(balances),
    matches, ids, covariates, timevary, mmlen, fmin
  )

  return model
end

### setup meanbalances
function setup_meanbalance(model::AbstractCICModel)
  
  MB = @chain model.matches[!, [:treattime, :treatunit, :matchunit, :f]] begin
  groupby([:treattime, :treatunit, :f])
    combine(
      # :f => Ref => :fs,
      :matchunit => Ref => :matchunits,
    )
    groupby([:treattime, :treatunit])
    combine(
      :f => Ref => :fs,
      :matchunits => Ref => :matchunitsets,
    )
  end;
  
  for covar in model.covariates
    if model.timevary[covar]
      MB[!, covar] = Vector{Vector{Vector{Union{Float64, Missing}}}}(undef, nrow(MB))
    else 
      MB[!, covar] = Vector{Vector{Union{Float64, Missing}}}(undef, nrow(MB))
    end
  end
  
  return MB
end

###

### meanbalances loop
function _meanbalance!(
  barrows, balrows,
  matches, ids, covariates, timevary, mmlen, fmin
)
  
  Threads.@threads for i in eachindex(balances)
    br = balrows[i]; mr = barrows[i];

    _, efsets = matchassignments(matches[i], ids);

    _covar_meanbalance!(
      mr, br, efsets, covariates, timevary, mmlen, fmin
    )
  end
  return barrows
end
  
### inner functions
function _covar_meanbalance!(
  mr, br, efsets, covariates, timevary, mmlen, fmin
)
  for covar in covariates
    # this is a [row, covar] for MB:
    #   the means across fs for a covariate (for a treated observation)
    if timevary[covar]
      row_covar_meanbalance!(
        mr[covar], br[covar], mr[:fs], efsets, mmlen, fmin
      );
    else
      row_covar_meanbalance!(
        mr[covar], br[covar], mr[:fs], efsets
      );
    end
  end
end

function row_covar_meanbalance!(
  Holding::Vector{Vector{Union{Missing, Float64}}},
  br_covar,
  fs, efsets, mmlen, fmin
)
  for (i, f) in enumerate(fs)
    if f # if there are any valid matches for that f
      ls = fidx(i + fmin - 1, mmlen, fmin);
      holding = Vector{Vector{Union{Float64, Missing}}}(
        undef, length(br_covar)
      );
      for m in eachindex(br_covar) # this is the site of averaging
        # for a single match
        if efsets[m][i]
          # if that match is valid for that f
          # else leave as missing
          holding[m] = br_covar[m][ls];
        end
      end
    end
    
    Holding[i] = Vector{Union{Missing, Float64}}(undef, mmlen);
    for l in eachindex(ls)
      lvec = skipmissing([vec[l] for vec in holding]);
      Holding[i][l] = !isempty(lvec) ? mean(lvec) : missing
    end
  end
  return Holding
end

function row_covar_meanbalance!(
  Holding, br_covar, fs, efsets,
)
  for (i, f) in eachindex(fs)
    if f
      holding = Vector{Union{Float64, Missing}}(undef, length(br_covar));
      for m in eachindex(br_covar) # this is the site of averaging
        if efsets[m][i]
          holding[m] = br_covar[m];
        end
      end
      Holding[i] = mean(skipmissing(holding))
    end
  end
  return Holding
end
