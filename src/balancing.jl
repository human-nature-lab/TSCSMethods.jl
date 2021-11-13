# balancing functions

"""
gives (1 / (std of treated units)) for each l in the matching period

outputs dict[t-l, covariate] where t-l is l days prior to the treatment

(adjust this for f, in sliding window F-defined match period case, which is the case...)
"""
function std_treated(model::AbstractCICModel, dat::DataFrame)

  trtobs = unique(dat[dat[!, model.treatment] .== 1, [model.t, model.id]])
  sort!(trtobs, [model.t, model.id])
  idx = Int[];
  L = Int[];

  mmin = model.L[begin]
  fmax = model.F[end]

  for i in 1:nrow(trtobs)
    to = trtobs[i, :]

    c1 = dat[!, model.id] .== to[2];
    ct = (dat[!, model.t] .>= (to[1] - 1 + mmin)) .& (dat[!, model.t] .<= to[1] + fmax);

    append!(idx, findall((c1 .& ct)))
    append!(L, dat[c1 .& ct, model.t] .- to[1]) # L relative to tt, not the actual time
  end

  allvals = @view dat[idx, model.covariates];

  Lset = unique(L);
  Lstd = zeros(Float64, length(Lset));
  Lstd = Dict{Tuple{Int64, Symbol}, Float64}()
  for l in Lset
    lvals = @view allvals[L .== l, :]
    for covar in model.covariates
      Lstd[(l, covar)] = inv(std(lvals[!, covar]; corrected = true))
    end
  end
  return Lstd 
end

function allocate_balances!(
  balances, tobs, covariates, timevary, matches, ids, Lrnge
);
  obslen = length(matches)

  for covar in covariates
    if timevary[covar]
      balances[!, covar] = Vector{
        Vector{Vector{Union{Missing, Float64}}}
        }(undef, obslen);
    else
      balances[!, covar] = Vector{Vector{Float64}}(undef, obslen);
    end
  end
  
  _fill_balances!(balances, tobs, ids, Lrnge, covariates, timevary)

  # don't bother storing these again
  # better to grab from matches, since it is fast enough
  # for i in eachindex(observations)
  #   balances.matchunits[i] = ids[matches[i].mus]
  # end;
  return balances
end

function _fill_balances!(balances, tobs, ids, Lrnge, covariates, timevary)
  for (i, balrw) in enumerate(eachrow(balances))
    emus = matchassignments(tobs[i], ids; returnefsets = false)
    __fill_balances!(
      balrw, emus, Lrnge, covariates, timevary
    )
  end
  return balances
end

function __fill_balances!(
  balrw, emus, Lrnge, covariates, timevary
)
  for covar in covariates
    if timevary[covar]
      balrw[covar] = [Vector{Union{Missing, Float64}}(missing, Lrnge) for i in 1:length(emus)]
    else
      balrw[covar] = Vector{Float64}(undef, length(emus))
    end
  end
  return balrw
end

"""
    fullbalance!(model::AbstractCICModel, dat::DataFrame)

calculate the full set of standardized balance scores, for each treated observation - matched unit pair.

cf. meanbalance for the averages by treated observation
"""
function fullbalance!(model::AbstractCICModel, dat::DataFrame)

  @unpack observations, matches, ids = model;
  @unpack F, L, id, t, treatment, covariates, timevary = model;
  fmin = minimum(F); fmax = maximum(F); mmin = minimum(L); mmax = maximum(L);
  Lrnge = length((fmin + mmin):(fmax + mmax));
  
  tmin = minimum(dat[!, t])
  cdat = Matrix(dat[!, covariates]);
  @unpack balances, reference = model;

  allocate_balances!(
    balances, matches, covariates, timevary, matches, ids, Lrnge
  );
  
  Lσ = std_treated(model, dat);

  tg, rg, trtg = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmin, fmax, mmin,
    cdat
  );

  _balance!(
    balances,
    observations,
    matches,
    ids,
    fmin, fmax, mmin, mmax,
    tmin,
    tg, rg,
    covariates, timevary,
    Lσ, reference
  );

  return model
end

function _balance!(
  balances,
  observations,
  tobs,
  ids,
  fmin, fmax, mmin, mmax,
  tmin,
  tg, rg,
  covariates, timevary,
  Lσ, reference
)

  @inbounds Threads.@threads for i in eachindex(observations)
    balrwi = @view balances[i, :];
    ob = observations[i];
    emus = matchassignments(tobs[i], ids; returnefsets = false);

    if length(emus) == 0
      continue
    end

    # we need the difference between the same covariates at each time point
    # (don't need each row)
    # we know which fs are allowable already, so we don't need to worry about lengths, etc. (since the panels are necessarily balanced)

    # fw = matchwindow(φ + fmin - 1, tt, mmin, mmax);
    ttmatch = mmin + ob[1] + fmin : mmax + ob[1] + fmax;

    γcs = eachcol(tg[ob]);

    # γrs = eachrow(tg[ob]);
    # γtimes = rg[ob];
    firstγtime = rg[ob][begin];

    # if the data vector is shorter than ttmatch
    # because of late data start
    # can ignore early end case (there won't be a
    # treatment known past the data...)
    tadjustment = if minimum(ttmatch) < tmin
      abs(minimum(ttmatch))
    else 0
    end

    balance_atreated!(
      balrwi, emus, ob[1], tg, covariates, timevary, γcs, Lσ,
      tadjustment, reference, firstγtime
    );

  end
  return tobs
end

function _timevarybalance!(
  distlocs, γc, gc, Lσ, tt, covar, tadjustment, xtfirst
)
  for (l, (x, y)) in enumerate(zip(γc, gc))
    xtfirst += 1
    distlocs[l+tadjustment] = (x - y) * Lσ[xtfirst - tt, covar]
  end
end

function balance_amatch!(
  balrwi, covariates, timevary,
  γcs, gcs, Lσ, tadjustment, reference, firstγtime, tt, em
)
  for (covar, γc, gc) in zip(covariates, γcs, gcs)
    # N.B. eachrow(balances)[i][covar]
    balic = balrwi[covar]
    if timevary[covar]
      _timevarybalance!(
        balic[em], γc, gc, Lσ, tt, covar, tadjustment, firstγtime
      )
    else
      balic[em] = (γc[1] - gc[1]) * Lσ[reference, covar]
      continue
    end
  end
  return balrwi
end

function balance_atreated!(
  balrwi, emus, tt, tg, covariates, timevary, γcs, Lσ,
  tadjustment, reference, firstγtime
)
  for (em, emu) in enumerate(emus)
    # emu = emus[e];
    # efs = efsets[em]; # logical rep. of fs that exist for a match    
    gcs = eachcol(tg[(tt, emu)]);

    # do all 80, for each covar, here
    balance_amatch!(
      balrwi, covariates, timevary,
      γcs, gcs, Lσ, tadjustment, reference, firstγtime, tt, em
    )
  end
  return balrwi
end

### setup meanbalances

function allocate_meanbalances!(model)

  @unpack F, L, covariates, timevary, matches, ids, observations, meanbalances = model;

  Len = length(L); Flen = length(F);

  # meanbalances = DataFrame(
  #   treattime = [ob[1] for ob in model.observations],
  #   treatunit = [ob[2] for ob in model.observations],
  #   matchunitsets = [mm for mm in model.matchunits]
  # )

  # need to check observations to see if there are any matches left

  meanbalances[!, :fs] = Vector{Vector{Bool}}(undef, length(matches));

  for covar in covariates
    if timevary[covar]
      meanbalances[!, covar] = Vector{Vector{Vector{Union{Float64, Missing}}}}(undef, length(matches))
    else 
      meanbalances[!, covar] = Vector{Vector{Union{Float64, Missing}}}(undef, length(matches))
    end
  end

  _fill_meanbalances!(
    meanbalances, matches, ids, Len, covariates, timevary, Flen
  );

  return model
end

function getfunion!(funion, efsets)
  for j in 1:length(funion)
    funion[j] = any([efset[j] for efset in efsets])
  end
end

function _fill_meanbalances!(
  meanbalances, matches, ids, Len, covariates, timevary, Flen
)
  
  for (i, balrw) in enumerate(eachrow(meanbalances))
    emus, efsets = matchassignments(matches[i], ids);
    # we want to create a vector for f, so long as at least one match unit allows it
    
    balrw[:fs] = Vector{Bool}(undef, Flen);    
    getfunion!(balrw[:fs], efsets)
    fpresent = sum(balrw[:fs])
  
    __fill_meanbalances!(
      balrw, fpresent, Len, covariates, timevary
    )
  end
  return meanbalances
end

function __fill_meanbalances!(
  balrw, fpresent, Len, covariates, timevary
)
  for covar in covariates
    if timevary[covar]
      balrw[covar] = [Vector{Union{Missing, Float64}}(missing, Len) for _ in 1:fpresent]
    else
      balrw[covar] = Vector{Float64}(undef, fpresent)
    end
  end
  return balrw
end

"""
    meanbalance!(model)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(model::AbstractCICModel, dat)

  @unpack meanbalances, observations, matches, ids, covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L, reference = model;

  fmin = minimum(F); fmax = maximum(F)
  mmin = minimum(L); mmax = maximum(L);

  tmin = minimum(dat[!, t]);
  
  allocate_meanbalances!(model);
  Lσ = std_treated(model, dat);

  tg, rg, trtg = make_groupindices(
    dat[!, t], dat[!, treatment],
    dat[!, id], ids,
    fmin, fmax, mmin,
    Matrix(dat[!, covariates]);
  );

  _meanbalance!(
    meanbalances,
    observations,
    matches,
    ids,
    fmin, mmin, mmax,
    tg, rg,
    covariates, timevary,
    reference, Lσ, tmin
  );

  return model
end

function _meanbalance!(
  meanbalances,
  observations,
  matches,
  ids,
  fmin, mmin, mmax,
  tg, rg,
  covariates, timevary,
  reference, Lσ, tmin
)

 # cnti = 0
  #for i in eachindex(observations)
   #cnti += 1
  @inbounds Threads.@threads for i in eachindex(observations)

    balrw = @view meanbalances[i, :];
    ob = (tt, tu) = observations[i];
    emus, efsets = matchassignments(matches[i], ids);

    if length(emus) == 0
      continue
    end

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
      balrw, tt, Ys, Yt, emus, efsets,
      Lσ, covariates, timevary,
      reference, fmin, mmin, mmax,
      tg, tmin
    );

    # _divbalance
    # alternative with broadcasting: time balrw[covar] ./ sum(efsets);
    efsum = sum(efsets);
    _meanmatch!(balrw, covariates, efsum);

  end
  return meanbalances
end

function _addmatches!(
  balrw, tt, Ys, Yt, emus, efsets,
  Lσ, covariates, timevary,
  reference, fmin, mmin, mmax,
  tg, tmin
)
#  emcount = 0
  for (em, emu) in enumerate(emus)
 #   emcount += 1
    # emu = emus[e];
    efs = efsets[em]; # logical rep. of fs that exist for a match
    Xs = eachcol(tg[(tt, emu)]);

    # cnt: since φ will track 1:31, and we will have only those that exist 
    _addmatch!(
      balrw, Ys, Xs, Yt, Lσ, efs, tt, covariates, timevary,
      reference, fmin, mmin, mmax, tmin
    );
  end
  return balrw
end

function _meanmatch!(balrw, covariates, efsum)
  for covar in covariates
      __meanmatch!(balrw[covar], efsum)
  end
  return balrw
end

function __meanmatch!(balrwcovar::Vector{Vector{Union{Missing, Float64}}}, efsum)
  for (φ, balφ) in enumerate(balrwcovar)
    for μ in eachindex(balφ)
      balφ[μ] = balφ[μ] / efsum[φ]
    end
  end
  return balrwcovar
end

function __meanmatch!(balrwcovar::Vector{Union{Missing, Float64}}, efsum)
  for μ in eachindex(balrwcovar)
    balrwcovar[μ] = balrwcovar[μ] / efsum[μ]
  end
  return balrwcovar
end

function _addmatch!(
  balrw, Ys, Xs, Yt, Lσ, efs, tt, covariates, timevary,
  reference, fmin, mmin, mmax, tmin
)
  for (c, (Y, X, covar)) in enumerate(zip(Ys, Xs, covariates))
    addmatch_f!(
      balrw[c+1], Y, X, Yt, Lσ, efs, tt, covar, timevary,
      reference, fmin, mmin, mmax, tmin
    )
  end
  return balrw
end

function addmatch_f!(
  balrwc, Y, X, Yt, Lσ, efs, tt, covar, timevary, reference,
  fmin, mmin, mmax, tmin
)
  φcnt = 0; # COUNTING ISSUE?
  # ocnt = 0
  for (φ, fb) in enumerate(efs)
    # ocnt += 1
    if fb
      φcnt += 1;
      fw = matchwindow(φ + fmin - 1, tt, mmin, mmax);
      fwadjustment = minimum(fw) < tmin ? abs(minimum(fw)) : 0;
      if timevary[covar]
        _timevarybalance!(
          balrwc[φcnt], Y, X, Yt, fw, Lσ, tt, covar, fwadjustment
        )
      else
        if ismissing(balrwc[φcnt])
          balrwc[φcnt] = 0.0
        else
          balrwc[φcnt] += (Y[1] - X[1]) * Lσ[reference, covar]
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
  downcount = 0 # when the thing is longer, we need to adjust the other way, or something
  for (l, (y, x, τ)) in enumerate(zip(Y, X, Yt))
    # timepoint must exist and be in window
    # if in window and does not exist, will remain missing
    if τ < minimum(fw)
      downcount += 1 # substract one for every period before the window
    end
    if τ ∈ fw
      # cntl += 1
      if ismissing(balrwcφ[l+fwadjustment-downcount])
        balrwcφ[l+fwadjustment-downcount] = 0.0
      else
        balrwcφ[l+fwadjustment-downcount] += (y - x) * Lσ[τ - tt, covar]
      end
    end
  end
  return balrwcφ
end

####################

"""
    grandbalance!(model::AbstractCICModel)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(model::AbstractCICModel)

  @unpack F, L, covariates, timevary = model;
  @unpack grandbalances, meanbalances = model;

  __grandbalance!(
    grandbalances, meanbalances,
    covariates, timevary, length(L)
  );

  return model
end

"""
    grandbalance!(model::AbstractCICModelStratified)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(stratmodel::AbstractCICModelStratified)

  @unpack F, L, covariates, timevary, strata = stratmodel;
  @unpack grandbalances, meanbalances = stratmodel;

  Len = length(L);

  if isempty(strata)
    println("N.B. stratum not defined in meanbalances. so, calculating nonstratified grandbalances.")
    @reset grandbalances = GrandDictNoStrat();
    __grandbalance!(
      grandbalances, meanbalances,
      covariates, timevary, length(L)
    );
  else
    us = sort(unique(strata));
    [grandbalances[s] = GrandDictNoStrat() for s in us];
    
    for s in us
      mb_s = @view meanbalances[strata .== s, :];
      __grandbalance!(
        grandbalances, mb_s,
        covariates, timevary, Len,
        s
      )
    end
  end
  return stratmodel
end

function __grandbalance!(
  grandbalances, meanbalances, covariates, timevary, Len, s
)
  for covar in covariates
    if timevary[covar]
      gbc = disallowmissing(
        _grandbalance(
          meanbalances[!, covar],
          Len
        )
      )
      grandbalances[s][covar] = gbc
    else
      grandbalances[s][covar] = _grandbalance(
        meanbalances[!, covar]
      )
    end
  end
  return grandbalances
end

function __grandbalance!(
  grandbalances, meanbalances, covariates, timevary, mmlen
)
  for covar in covariates
    if timevary[covar]
      gbc = disallowmissing(
        _grandbalance(
          meanbalances[!, covar],
          mmlen
        )
      )
      grandbalances[covar] = gbc
    else
      grandbalances[covar] = _grandbalance(
        meanbalances[!, covar]
      )
    end
  end
  return grandbalances
end

function _grandbalance(covec, mmlen)
  reduced = reduce(vcat, covec);
  means = Vector{Union{Missing, Float64}}(undef, mmlen);
  for l in eachindex(1:mmlen)
    means[l] = mean(skipmissing([vec[l] for vec in reduced]))
  end
  return means
end

function _grandbalance(covec)
  return mean(skipmissing(reduce(vcat, covec)))
end

function balance!(model::AbstractCICModel, dat)

  meanbalance!(model, dat);
  grandbalance!(model);

  return model
end
