# new balancing functions

"""
gives (1 / (std of treated units)) for each l in the matching period

outputs dict[t-l, covariate] where t-l is l days prior to the treatment

(adjust this for f, in sliding window F-defined match period case, which is the case...)
"""
function std_treated(cc::AbstractCICModel, dat::DataFrame)

  trtobs = unique(dat[dat[!, cc.treatment] .== 1, [cc.t, cc.id]])
  sort!(trtobs, [cc.t, cc.id])
  idx = Int[];
  L = Int[];

  for i in 1:nrow(trtobs)
    to = trtobs[i, :]

    c1 = dat[!, cc.id] .== to[2];
    ct = (dat[!, cc.t] .>= (to[1] - 1 + cc.mmin)) .& (dat[!, cc.t] .<= to[1] + cc.fmax);

    append!(idx, findall((c1 .& ct)))
    append!(L, dat[c1 .& ct, cc.t] .- to[1]) # L relative to tt, not the actual time
  end

  allvals = @view dat[idx, cc.covariates];

  Lset = unique(L);
  Lstd = zeros(Float64, length(Lset));
  Lstd = Dict{Tuple{Int64, Symbol}, Float64}()
  for l in Lset
    lvals = @view allvals[L .== l, :]
    for covar in cc.covariates
      Lstd[(l, covar)] = inv(std(lvals[!, covar]; corrected = true))
    end
  end
  return Lstd 
end

# sgi = gi[r, sv]
# Ls =  Lσ[cc.reference, sv]
function svcalc(iusv1, jusv1, Ls)
  return (iusv1 - jusv1) * Ls
end

function tvcalc(iul, jul, Ls)
  return (iul - jul) * Ls;
end

"""
    fullbalance!(cc::AbstractCICModel, dat::DataFrame)

calculate the full set of standardized balance scores, for each treated observation - matched unit pair.

cf. meanbalance for the averages by treated observation
"""
function fullbalance!(cc::AbstractCICModel, dat::DataFrame)

  varset = vcat([cc.t, cc.id], cc.covariates)
  
  sdat = dat[!, varset]
  tmin = minimum(dat[!, cc.t]);
  tmax = maximum(dat[!, cc.t]);

  Lσ = std_treated(cc, dat);

  cc.balances = unique(cc.matches[!, [:treattime, :treatunit, :matchunit]])
  
  Lrnge = length((cc.fmin + cc.mmin):(cc.fmax + cc.mmax))
  
  timevar = Vector{Symbol}();
  staticvar = Vector{Symbol}();
  for covar in cc.covariates
    if cc.timevary[covar]
      push!(timevar, covar)
      cc.balances[!, covar] = [Vector{Union{Missing, Float64}}(missing, Lrnge) for i in 1:nrow(cc.balances)]
    else
      push!(staticvar, covar)
      cc.balances[!, covar] = zeros(Float64, nrow(cc.balances));
    end
  end

  gdf = groupby(cc.balances, [:treattime, :treatunit]);

  _balance!(
    cc, sdat, Lσ,
    staticvar, timevar,
    tmin, tmax,
    gdf
  );

  return cc
end

function _balance!(
  cc::AbstractCICModel, sdat::DataFrame, Lσ,
  staticvar, timevar,
  tmin, tmax,
  gdf
)
  
  Threads.@threads for gi in gdf
    
    iu = gi[1, :treatunit]
    tt = gi[1, :treattime]

    ttl = (tt + cc.fmin + cc.mmin);
    ttu = (tt + cc.fmax - 1);
    sdf = @view sdat[(sdat[!, cc.t] .>= ttl) .& (sdat[!, cc.t] .<= ttu), :];

    iudat = @view sdf[sdf[!, cc.id] .== iu, :];

    _row_balance!(
      gi, sdf, cc, staticvar, timevar, iudat, tt, ttl, ttu, tmin, tmax, Lσ
    )
  end
  return cc
end

function _row_balance!(
  gi, sdf, cc, staticvar, timevar, iudat, tt, ttl, ttu, tmin, tmax, Lσ
)
  # @eachrow! cc.balances[parentindices(gi)[1], :] begin
  # cgi = @view cc.balances[parentindices(gi)[1], :];
  for r in eachrow(gi)
    judat = @view sdf[sdf[!, cc.id] .== r[:matchunit], :];

    # if unit series aren't all of same length, this would be an issue
    for sv in staticvar
      r[sv] = svcalc(iudat[1, sv], judat[1, sv], Lσ[cc.reference, sv])
      # cc.reference (tt - 1): a value that must exist, doesn't matter since it's static...
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
    meanbalance(cc!)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in cc.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(ccr::AbstractCICModel)
  ccr.meanbalances, groupedbalances = setup_meanbalance(ccr);
  mmlen = length(ccr.mmin:ccr.mmax);
  _meanbalance!(ccr.meanbalances, groupedbalances, mmlen, ccr.fmin)
  return ccr
end

### setup meanbalances
function setup_meanbalance(ccr::AbstractCICModel)
  
  MB = @chain ccr.matches[!, [:treattime, :treatunit, :matchunit, :f]] begin
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
  
  for covar in ccr.covariates
    if cc.timevary[covar]
      MB[!, covar] = Vector{Vector{Vector{Union{Float64, Missing}}}}(undef, nrow(MB))
    else 
      MB[!, covar] = Vector{Vector{Union{Float64, Missing}}}(undef, nrow(MB))
    end
  end
  
  GB = groupby(ccr.balances, [:treattime, :treatunit, :matchunit])

  return MB, GB
end

### meanbalances loop
function _meanbalance!(meanbalances, groupedbalances, mmlen, fmin)
  for r in eachrow(meanbalances)
    # r = eachrow(MB)[1];
    fs = r[:fs];
    muset = r[:matchunitsets];

    for covar in ccr.covariates
      # this is a [row, covar] for MB:
      #   the means across fs for a covariate (for a treated observation)
      if ccr.timevary[covar]
        r[covar] = Vector{Vector{Union{Float64, Missing}}}(undef, length(fs));
        row_covar_meanbalance!(
          r[covar], r[:treattime], r[:treatunit],
          fs, muset, groupedbalances, mmlen, fmin, covar
        );
      else
        r[covar] = Vector{Union{Float64, Missing}}(undef, length(fs));
        row_covar_meanbalance!(
          r[covar], r[:treattime], r[:treatunit],
          fs, muset, groupedbalances, covar
        );
      end
    end
  end
  return meanbalances
end

### inner functions
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
      Holding[i][l] = mean(skipmissing([vec[l] for vec in holding])) # assign to part of MB row, variable (for an f)
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

  mmlen = length(ccr.mmin:ccr.mmax)

  if model.stratifier != Symbol("")
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
  return model
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
  return model
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
  cc::AbstractCICModel;
  threshold = 0.1,
  stratareduce = true
)

  if cc.stratifier == Symbol("")
    chk = Dict{Symbol, Bool}();
    for (k, v) in cc.grandbalances
      _balancecheck!(chk, v, k, threshold)
    end
  else
    chk = Dict{Int, Dict{Symbol, Bool}}();
    [chk[s] = Dict{Symbol, Bool}() for s in keys(cc.grandbalances)]
    for (k, v) in cc.grandbalances
      for (ki, vi) in v
      _balancecheck!(chk[k], vi, ki, threshold)
      end
    end
  end

  if stratareduce & (cc.stratifier != Symbol(""))
    bc2 = Dict{Symbol, Bool}()
    for covar in cc.covariates
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

function balance!(cc::cicmodel, dat::DataFrame)

  fullbalance!(cc, dat);
  meanbalance!(cc);
  grandbalance!(cc);
  
  if cc.stratifier == Symbol("")
    cc.treatednum = nrow(unique(cc.meanbalances[!, [:treattime, :treatunit]]));
  else
    cc.treatednum = Dict{Int64, Int64}();
    for s in unique(cc.meanbalances.stratum)
      c1 = cc.meanbalances.stratum .== s
      cc.treatednum[s] = nrow(unique(@views(cc.meanbalances[c1, [:treattime, :treatunit]])));
    end
  end

  return cc
end

"""
    autobalance(
      cc;
      refinementnum = 5,
      calmin = 0.08, step = 0.05, initial_bals = false
    )

Automatically balance via a simple algorithm. Start with initial caliper of 1.0, and subtract `step` whenever the grand mean balance threshold (0.1) is not met.
"""
function autobalance(
  cc;
  threshold = 0.1,
  min_treated_obs = 10,
  refinementnum = 5, calmin = 0.1, step = 0.05, initial_bals = false
)


if !initial_bals
  caliper = Dict{Symbol, Float64}();
  for c in cc.covariates
    caliper[c] = 1.0
  end
end

  cal = make_caliper(cc, caliper);
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
    cal = make_caliper(cc, caliper);
    calr = make_refined(cal; refinementnum = refinementnum);
    bc = balancecheck(calr; threshold = threshold)
    insight = inspectcaliper(calr);
  end

  return cal, calr
end
