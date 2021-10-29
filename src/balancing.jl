# new balancing functions

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

  for i in 1:nrow(trtobs)
    to = trtobs[i, :]

    c1 = dat[!, model.id] .== to[2];
    ct = (dat[!, model.t] .>= (to[1] - 1 + model.mmin)) .& (dat[!, model.t] .<= to[1] + model.fmax);

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

# sgi = gi[r, sv]
# Ls =  Lσ[model.reference, sv]
function svcalc(iusv1, jusv1, Ls)
  return (iusv1 - jusv1) * Ls
end

function tvcalc(iul, jul, Ls)
  return (iul - jul) * Ls;
end

"""
    fullbalance!(model::AbstractCICModel, dat::DataFrame)

calculate the full set of standardized balance scores, for each treated observation - matched unit pair.

cf. meanbalance for the averages by treated observation
"""
function fullbalance!(model::AbstractCICModel, dat::DataFrame)

  varset = vcat([model.t, model.id], model.covariates)
  
  sdat = dat[!, varset]
  tmin = minimum(dat[!, model.t]);
  tmax = maximum(dat[!, model.t]);

  Lσ = std_treated(model, dat);

  model.balances = unique(model.matches[!, [:treattime, :treatunit, :matchunit]])
  
  Lrnge = length((model.fmin + model.mmin):(model.fmax + model.mmax))
  
  timevar = Vector{Symbol}();
  staticvar = Vector{Symbol}();
  for covar in model.covariates
    if model.timevary[covar]
      push!(timevar, covar)
      model.balances[!, covar] = [Vector{Union{Missing, Float64}}(missing, Lrnge) for i in 1:nrow(model.balances)]
    else
      push!(staticvar, covar)
      model.balances[!, covar] = zeros(Float64, nrow(model.balances));
    end
  end

  gdf = groupby(model.balances, [:treattime, :treatunit]);

  _balance!(
    model, sdat, Lσ,
    staticvar, timevar,
    tmin, tmax,
    gdf
  );

  return model
end

function _balance!(
  model::AbstractCICModel, sdat::DataFrame, Lσ,
  staticvar, timevar,
  tmin, tmax,
  gdf
)
  
  Threads.@threads for gi in gdf
    
    iu = gi[1, :treatunit]
    tt = gi[1, :treattime]

    ttl = (tt + model.fmin + model.mmin);
    ttu = (tt + model.fmax - 1);
    sdf = @view sdat[(sdat[!, model.t] .>= ttl) .& (sdat[!, model.t] .<= ttu), :];

    iudat = @view sdf[sdf[!, model.id] .== iu, :];

    _row_balance!(
      gi, sdf, model, staticvar, timevar, iudat, tt, ttl, ttu, tmin, tmax, Lσ
    )
  end
  return model
end

function _row_balance!(
  gi, sdf, model, staticvar, timevar, iudat, tt, ttl, ttu, tmin, tmax, Lσ
)
  # @eachrow! model.balances[parentindices(gi)[1], :] begin
  # cgi = @view model.balances[parentindices(gi)[1], :];
  for r in eachrow(gi)
    judat = @view sdf[sdf[!, model.id] .== r[:matchunit], :];

    # if unit series aren't all of same length, this would be an issue
    for sv in staticvar
      r[sv] = svcalc(iudat[1, sv], judat[1, sv], Lσ[model.reference, sv])
      # model.reference (tt - 1): a value that must exist, doesn't matter since it's static...
    end

    lcount = 0
    for l ∈ ttl:ttu
      if (l >= tmin) & (l <= tmax) # should be for iu, ju?
        lcount += 1
        # iul = iudat[lcount, :]
        # jul = judat[lcount, :]

        for tv in timevar
          r[tv][l - ttl + 1] = tvcalc(
              iudat[lcount, tv], judat[lcount, tv], Lσ[l - tt, tv]
            )
        end
      end
    end
  end
  return gi
end

"""
    fidx(f::Int, mlen::Int, fmin::Int)

get indices from f
"""
fidx(f::Int, mlen::Int, fmin::Int) = (f - fmin + 1):(f - fmin + mlen);

"""
    meanbalance!(model)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in model.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(model::AbstractCICModel)
  model.meanbalances, groupedbalances = setup_meanbalance(model);
  mmlen = length(model.mmin:model.mmax);
  _meanbalance!(
    model.meanbalances, groupedbalances,
    model.covariates, model.timevary, mmlen, model.fmin
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
  
  GB = groupby(model.balances, [:treattime, :treatunit, :matchunit])

  return MB, GB
end

### meanbalances loop
function _meanbalance!(
  meanbalances, groupedbalances, covariates, timevary, mmlen, fmin
)
  # for r in eachrow(meanbalances)
  Threads.@threads for r in eachrow(meanbalances)
    fslen = length(r[:fs]);
    _covar_meanbalance!(
      r, fslen, groupedbalances, covariates, timevary, mmlen, fmin
    )
  end
  return meanbalances
end

### inner functions
function _covar_meanbalance!(
  r, fslen, groupedbalances, covariates, timevary, mmlen, fmin
)
  for covar in covariates
    # this is a [row, covar] for MB:
    #   the means across fs for a covariate (for a treated observation)
    if timevary[covar]
      r[covar] = Vector{Vector{Union{Float64, Missing}}}(undef, fslen);
      row_covar_meanbalance!(
        r[covar], r[:treattime], r[:treatunit],
        r[:fs], r[:matchunitsets], groupedbalances, mmlen, fmin, covar
      );
    else
      r[covar] = Vector{Union{Float64, Missing}}(undef, fslen);
      row_covar_meanbalance!(
        r[covar], r[:treattime], r[:treatunit],
        r[:fs], r[:matchunitsets], groupedbalances, covar
      );
    end
  end
end

function row_covar_meanbalance!(
  Holding::Vector{Vector{Union{Missing, Float64}}},
  treattime, treatunit,
  fs, muset, GB, mmlen, fmin, covar
)
  for (i, f) in enumerate(fs)

    ls = fidx(f, mmlen, fmin);
    mus = muset[i];
    
    holding = Vector{Vector{Union{Float64, Missing}}}(undef, length(mus));
    for m in eachindex(mus) # this is the site of averaging
      g = get(
        GB,
        (treattime = treattime, treatunit = treatunit, matchunit = mus[m]),
        nothing
      );
      holding[m] = g[1, covar][ls] # for a single match
    end
    
    Holding[i] = Vector{Union{Missing, Float64}}(undef, mmlen);
    for l in eachindex(ls)
      # lvec = [vec[l] for vec in holding]; # place the values for match period l for all matches to a treated obs into a single vector.
      # assign to part of MB row, variable (for an f)
      # Holding[i][l] = try mean(skipmissing([vec[l] for vec in holding])) catch; missing end
      # Holding[i][l] = all(ismissing.(lvec)) ? missing : mean(skipmissing(lvec))
      lvec = skipmissing([vec[l] for vec in holding]);
      Holding[i][l] = !isempty(lvec) ? mean(lvec) : missing
    end
  end
  return Holding
end
  
function row_covar_meanbalance!(
  Holding, treattime, treatunit, fs, muset, GB, covar
)
  for i in eachindex(fs)

    mus = muset[i];
    
    holding = Vector{Union{Float64, Missing}}(undef, length(mus));
    for m in eachindex(mus) # this is the site of averaging
      g = get(
        GB,
        (treattime = treattime, treatunit = treatunit, matchunit = mus[m]),
        nothing
      );
      holding[m] = g[1, covar] # for a single match
    end

    Holding[i] = mean(skipmissing(holding))
  end
  return Holding
end

"""
    grandbalance!(model::AbstractCICModel)

Calculate the overall mean covariate balance for a model.
"""
function grandbalance!(model::AbstractCICModel)

  mmlen = length(model.mmin:model.mmax)

  cs1 = model.stratifier != Symbol("")
  cs2 = "stratum" ∈ names(model.meanbalances)

  if (cs1 & !cs2)
    println("N.B. stratum not defined in meanbalances. so, calculating nonstratified grandbalances.")
  end

  if (cs1 & cs2)
    model.grandbalances = GrandDictStrat();
    us = sort(unique(model.meanbalances.stratum));
    [model.grandbalances[s] = GrandDictNoStrat() for s in us];
    
    for s in us
      mb_s = @view(model.meanbalances[model.meanbalances.stratum .== s, :]);
      __grandbalance!(
        model.grandbalances, mb_s,
        model.covariates, model.timevary, mmlen,
        s
      )
    end

  else
    model.grandbalances = GrandDictNoStrat();
    __grandbalance!(
      model.grandbalances, model.meanbalances,
      model.covariates, model.timevary, mmlen
    )
  end

  return model
end

function __grandbalance!(
  grandbalances, meanbalances, covariates, timevary, mmlen, s
)
  for covar in covariates
    if timevary[covar]
      gbc = disallowmissing(
        _grandbalance(
          meanbalances[!, covar],
          mmlen
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

"""
    balancecheck(
      m::Dict{Symbol, Union{Float64, Vector{Float64}}};
      threshold = 0.1, stratareduce = true
    )

Simply check whether the grand means are above the std. balance threshold. Returns a Bool for each covariate.
If Stratareduce is true, then the strata balances will be agggregated to the covariate level, such that a violation in any caliper triggers a violation in the aggregated output.
"""
function balancecheck(
  model::AbstractCICModel;
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

function balance!(model::cicmodel, dat::DataFrame)

  fullbalance!(model, dat);
  meanbalance!(model);
  grandbalance!(model);
  
  if model.stratifier == Symbol("")
    model.treatednum = nrow(unique(model.meanbalances[!, [:treattime, :treatunit]]));
  else
    model.treatednum = Dict{Int64, Int64}();
    for s in unique(model.meanbalances.stratum)
      c1 = model.meanbalances.stratum .== s
      model.treatednum[s] = nrow(unique(@views(model.meanbalances[c1, [:treattime, :treatunit]])));
    end
  end

  return model
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
  caliper = Dict{Symbol, Float64}();
  for c in model.covariates
    caliper[c] = 1.0
  end
end

  cal = make_caliper(model, caliper);
  calr = make_refined(cal; refinementnum = refinementnum);

  # check calr only
  bc = balancecheck(calr; threshold = threshold)
  insight = inspectcaliper(calr);

  while any(values(bc)) & (nrow(insight) >= min_treated_obs) & all(values(caliper) .>= calmin)

    covset = [k for k in keys(bc)];
    for covar in tscsmethods.sample(covset, length(covset); replace = false)
      if bc[covar] & (caliper[covar] > calmin)
        caliper[covar] = caliper[covar] - step
      end
    end
    cal = make_caliper(model, caliper);
    calr = make_refined(cal; refinementnum = refinementnum);
    bc = balancecheck(calr; threshold = threshold)
    insight = inspectcaliper(calr);
  end

  return cal, calr
end
