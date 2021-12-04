# balancing.jl

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

### setup meanbalances

function allocate_meanbalances!(model)

  @unpack observations, matches, ids, meanbalances = model;
  @unpack covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L = model;

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
      balrw[covar] = Vector{Union{Missing, Float64}}(missing, fpresent)
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

  tg, rg, _ = make_groupindices(
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
    mmin, mmax,
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
function meanbalance!(model::AbstractCICModel, dat, tg, rg)

  @unpack meanbalances, observations, matches, ids, covariates, timevary = model;
  @unpack t, id, treatment = model;
  @unpack F, L, reference = model;

  mmin = minimum(L); mmax = maximum(L);

  tmin = minimum(dat[!, t]);
  
  allocate_meanbalances!(model);
  Lσ = std_treated(model, dat);

  _meanbalance!(
    meanbalances,
    observations,
    matches,
    ids,
    mmin, mmax,
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
  mmin, mmax,
  tg, rg,
  covariates, timevary,
  reference, Lσ, tmin
)

  @inbounds Threads.@threads for i in eachindex(observations)

    balrw = @view meanbalances[i, :];
    ob = (tt, _) = observations[i];
    emus, efsets = matchassignments(matches[i], ids);
    efsum = sum(efsets);
    bwpres = efsum .> 0
    fset = F[bwpres]

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
      balrw, tt, Ys, Yt, emus, efsets, fset, bwpres,
      Lσ, covariates, timevary,
      reference, mmin, mmax,
      tg, tmin,
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
  balrw, tt, Ys, Yt, emus, efsets, fset, bwpres,
  Lσ, covariates, timevary,
  reference, mmin, mmax,
  tg, tmin
)

  for (em, emu) in enumerate(emus)
    efs = efsets[em][bwpres];
       # logical rep. of fs that exist for a match, then
       # reindex by bwpres, since balrw only contains bwpres 
       # entries (out of flen)
    Xs = eachcol(tg[(tt, emu)]);

    # cnt: since φ will track 1:31, and we will have only those that exist 
    _addmatch!(
      balrw, Ys, Xs, Yt, Lσ, efs, fset,
      tt, covariates, timevary,
      reference, mmin, mmax, tmin
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

function __meanmatch!(
  balrwcovar::Vector{Vector{Union{Missing, Float64}}},
  efsum
)

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
  balrw, Ys, Xs, Yt, Lσ, efs, fset,
  tt, covariates, timevary,
  reference, mmin, mmax, tmin
)
  for (c, (Y, X, covar)) in enumerate(zip(Ys, Xs, covariates))
    addmatch_f!(
      balrw[c+1], Y, X, Yt, Lσ, efs, fset, tt, covar, timevary,
      reference, mmin, mmax, tmin
    )
  end
  return balrw
end

function addmatch_f!(
  balrwc,
  Y, X, Yt, Lσ, efs, fset, tt, covar, timevary, reference,
  mmin, mmax, tmin
)

  for ((φ, fb), f) in zip(enumerate(efs), fset)
    # ANOTHER COUNTING/INDEXING ISSUE
    # we restrict to the cases where there is at least one f, so we need to index accordingly
    # efs[bwpres] is the correct object to index
    
    if fb
      fw = tscsmethods.matchwindow(f, tt, mmin, mmax);
      fwadjustment = minimum(fw) < tmin ? abs(minimum(fw)) : 0;
      if timevary[covar]
        _timevarybalance!(
          balrwc[φ], # mcounts[c][φ],
          Y, X, Yt, fw, Lσ, tt, covar, fwadjustment
        )
      else
        if ismissing(balrwc[φ])
          balrwc[φ] = 0.0
        else
          balrwc[φ] += (Y[1] - X[1]) * Lσ[reference, covar]
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
  downcount = 0 # when the thing is longer, we need to adjust the other way
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

# stratified
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

# not stratified
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

#####

#= 
covar = vn.pd
covec = model.meanbalances[!, covar];

X = [false for _ in 1:length(covec)];
for (i, cove) in enumerate(covec)
  if any(isnan.(cove))
    X[i] = true
  end
end

findall(X)

# this is random

 pop density NAN values
4-element Vector{Int64}:
  716
 1951
 2407
 2413

corresponding to observations
 (2, 48257)
 (93, 42083)
 (114, 36053)
 (114, 36065)
=#