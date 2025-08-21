# meanbalancing.jl

"""
    meanbalance!(model)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(model::VeryAbstractCICModel, dat)

  # DEV
  # import TSCSMethods:allocate_meanbalances!,std_treated,make_groupindices,_meanbalance!
  # using Parameters

  @unpack meanbalances, observations, matches, ids, covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L, reference = model;

  fmax = maximum(F);
  Lmin, Lmax = extrema(L);

  tmin = minimum(dat[!, t]);
  
  allocate_meanbalances!(model);
  Lσ = std_treated(model, dat);

  tg, rg, _ = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmax, Lmin,
    Matrix(dat[!, covariates]);
  );

  _meanbalance!(
    meanbalances,
    observations,
    matches,
    ids,
    F, Lmin, Lmax,
    tg, rg,
    covariates, timevary,
    reference, Lσ, tmin
  );

  return model
end

"""
    meanbalance!(model, dat, tg, rg)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(model::VeryAbstractCICModel, dat, tg, rg)

  @unpack meanbalances, observations, matches, ids, covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L, reference = model;

  Lmin, Lmax = extrema(L);

  tmin = minimum(dat[!, t]);
  
  allocate_meanbalances!(model);
  Lσ = std_treated(model, dat);

  _meanbalance!(
    meanbalances,
    observations,
    matches,
    ids,
    F, Lmin, Lmax,
    tg, rg,
    covariates, timevary,
    reference, Lσ, tmin
  );

  return model
end

## calculation

function _meanbalance!(
  meanbalances,
  observations,
  matches,
  ids,
  F, Lmin, Lmax,
  tg, rg,
  covariates, timevary,
  reference, Lσ, tmin
)

  # import TSCSMethods:getfunion!
  # for i in eachindex(observations)
  @inbounds Threads.@threads :greedy for i in eachindex(observations)
  
    matches_i = @views matches[i]
    @unpack mus = matches_i;

    # if any(mus)
    #   continue
    # end

    balrw = @view meanbalances[i, :];
    ob = (tt, _) = observations[i];
  
    # this should be recycled with @init and floops
    bwpres = Vector{Bool}(undef, length(F));
    getfunion!(bwpres, mus);

    efsum = Vector{Int}(undef, length(F));
    for (cnt, mc) in enumerate(eachcol(mus)); efsum[cnt] = sum(mc) end

    # we need the difference between the same covariates at each time point
    # (don't need each row)
    # we know which fs are allowable already, so we don't need to worry about lengths, etc. (since the panels are necessarily balanced)

    Ys = eachcol(tg[ob]);

    # γrs = eachrow(tg[ob]);
    Yt = rg[ob];

    # total match window across all fs
    # not needed, calculate fw locally
    # ttmatch = mmin + ob[1] + fmin : mmax + ob[1] + fmax;

    # if the data vector is shorter than ttmatch
    # because of late data start
    # can ignore early end case (there won't be a
    # treatment known past the data...)
    # tadjustment = if minimum(ttmatch) < tmin
    #   abs(minimum(ttmatch))
    # else 0
    # end
    
    _addmatches!(
      balrw, tt, Ys, Yt, ids, eachrow(mus),
      Lσ, covariates, timevary,
      reference, Lmin, Lmax,
      tg, tmin
    );

    # testing: resolved -- was undef instead of missing
    # if any(isnan.(balrw[vn.pd]))
    #   error("nan found in addition " * string(i))
    # end

    # _divbalance
    # alternative with broadcasting: time balrw[covar] ./ sum(efsets);
    _meanmatch!(balrw, covariates, efsum[bwpres]);

    # if any(isnan.(balrw[vn.pd]))
    #   error("nan found after division " * string(i))
    # end

  end
  return meanbalances
end

function _addmatches!(
  balrw, tt, Ys, Yt, ids, musrows,
  Lσ, covariates, timevary,
  reference, Lmin, Lmax,
  tg, tmin
)

  # musrows = eachrow(mus)
  # (m, murow) = collect(enumerate(musrows))[327];
  for (m, murow) in enumerate(musrows)
    if any(murow)
      mu = ids[m]; # matches[1].mus
      Xs = eachcol(tg[(tt, mu)]);

      # cnt: since φ will track 1:31, and we will have only those that exist 
      _addmatch!(
        balrw, Ys, Xs, Yt, Lσ, murow,
        tt, covariates, timevary, reference,
        Lmin, Lmax, tmin
      );
    end
  end
  return balrw
end

function _addmatch!(
  balrw, Ys, Xs, Yt, Lσ,
  murow,
  tt, covariates, timevary,
  reference,
  Lmin, Lmax, tmin
)

  # fmin would be needed for sliding, crossover-
  # window-based matching
  # just leave input in place

  # (c, (Y, X, covar)) = collect(enumerate(zip(Ys, Xs, covariates)))[1]
  for (c, (Y, X, covar)) in enumerate(zip(Ys, Xs, covariates))
    addmatch_f!(
      balrw[c+1], # first col is fs
      Y, X, Yt, Lσ, murow,
      tt, covar, timevary, reference,
      Lmin, Lmax, tmin
    )
  end
  return balrw
end

function addmatch_f!(
  balrwc,
  Y, X, Yt, Lσ, murow,
  tt, covar, timevary, reference,
  Lmin, Lmax, tmin
)
  # balrwc = balrw[c+1]
  fbcnt = 0 # index balrwc, since it will only have vectors for included fs

  # fb = true
  for fb in murow
    if fb
      fbcnt += 1
      
      # sliding
      # implement later
      # fw = matchwindow(φ + fmin - 1, tt, mmin, mmax);

      # static
      fw = Lmin + tt : Lmax + tt

      fwadjustment = minimum(fw) < tmin ? abs(minimum(fw)) : 0;
      if timevary[covar]
        _timevarybalance!(
          balrwc[fbcnt],
          Y, X, Yt, fw, Lσ, tt, covar, fwadjustment
        )
      else
        if ismissing(balrwc[fbcnt])
          balrwc[fbcnt] = 0.0
        else
          balrwc[fbcnt] += (Y[1] - X[1]) * Lσ[reference, covar]
        end
        continue
      end
    end
  end
  return balrwc
end

function _timevarybalance!(
  balrwcφ, Y, X, Yt, fw, Lσ, tt, covar, fwadjustment
)
  # balrwcφ = balrwc[fbcnt]

  downcount = 0 # when the thing is longer, we need to adjust the other way

  # (l, (y, x, τ)) = collect(enumerate(zip(Y, X, Yt)))[1]

  for (l, (y, x, τ)) in enumerate(zip(Y, X, Yt))
  # (l, (y, x, τ)) = collect(enumerate(zip(Y, X, Yt)))[1]
    # timepoint must exist and be in window
    # if in window and does not exist, will remain missing
    if τ < minimum(fw)
      downcount += 1 # substract one for every period before the window
    end
    if τ ∈ fw
      # if the timepoint in Yt is in the window
      # (recall we have a Yt for the period over whole window range across fs)
      # this is going to be an addition
      # mcountscφ[l+fwadjustment-downcount] += 1
      if ismissing(balrwcφ[l + fwadjustment - downcount])
        balrwcφ[l + fwadjustment - downcount] = 0.0
      end
      balrwcφ[l + fwadjustment - downcount] += (y - x) * Lσ[τ - tt, covar]
    end
  end
  return balrwcφ #, mcountscφ
end

function _meanmatch!(balrw, covariates, efsum)
  for covar in covariates
      __meanmatch!(balrw[covar], efsum)
  end
  return balrw
end

function __meanmatch!(
  balrwcovar::Vector{BalanceData},
  efsum
)
  for (φ, balφ) in enumerate(balrwcovar)
    for μ in eachindex(balφ)
      if !balφ.is_missing[μ]
        balφ.values[μ] = balφ.values[μ] / efsum[φ]
      end
    end
  end
  return balrwcovar
end

function __meanmatch!(balrwcovar::BalanceData, efsum)
  for μ in eachindex(balrwcovar)
    if !balrwcovar.is_missing[μ]
      balrwcovar.values[μ] = balrwcovar.values[μ] / efsum[μ]
    end
  end
  return balrwcovar
end
