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
    meanbalance(cc)

Calculate the mean balances, for each treated observation from the full set of balances. This will limit the calculations to include those present in cc.matches, e.g. in case a caliper has been applied.
"""
function meanbalance!(cc::AbstractCICModel)
  mlen = length(cc.mmin:cc.mmax);
  flen = length(cc.fmin:cc.fmax);
  G, cc.meanbalances, ix, tts, tus = _meanbalance_setup(cc, mlen);
  gloop!(cc.meanbalances, ix, cc, G, mlen, flen, tts, tus);
  return cc.meanbalances
end

function _meanbalance_setup(cc::AbstractCICModel, mlen::Int)

  # do we need to sort this here, or can it be done locally?
  # it is not used again directly
  # sort!(cc.matches, [:treattime, :treatunit, :matchunit, :f, :mdist]);
  
  G = @chain cc.matches[!, [:treattime, :treatunit, :matchunit, :f]] begin
    groupby([:treattime, :treatunit, :matchunit])
    combine(:f => Ref∘unique => :fs)
    @orderby(:treattime, :treatunit, :matchunit)
    leftjoin(
      cc.balances,
      on = [:treattime, :treatunit, :matchunit]
    )
    groupby([:treattime, :treatunit])
  end;

  # indexing problem solvable by including the number of fs for each treated unit, then have rolling sum
  
  # G and M will be in the same order
  # G = groupby(C, [:treattime, :treatunit]);
  
  M = @chain cc.matches[!, [:treattime, :treatunit, :f]] begin
    unique([:treattime, :treatunit, :f])
    sort([:treattime, :treatunit, :f])
  end
    #= must be sorted to match G
       there will be an f missing for some (tt, tu) iff there are no exemplars
       in matches
       in which case, that row will not be called upon anyway
    =#
  
  for covar in cc.covariates
    if cc.timevary[covar]
      M[!, covar] = [zeros(Float64, mlen) for i in 1:nrow(M)];
      # if we want to store the numbers
      # M[!, Symbol(string(covar) * "N")] = [zeros(Int64, mlen) for i in 1:nrow(M)];
    else
      M[!, covar] = zeros(Float64, nrow(M));
      # M[!, Symbol(string(covar) * "N")] = zeros(Int64, nrow(M));
    end
  end

  # make index vector
  ix = gindices(M, G);

  Mg = unique(M[!, [:treattime, :treatunit]]);

  return G, M, ix, Mg[!, :treattime], Mg[!, :treatunit]
end
  
function gindices(M, G)

  rns = (@chain M begin
    groupby([:treattime, :treatunit])
    combine(nrow => :rn)
  end).rn

  ix = Vector{Int64}(undef, length(G) + 1);
  ix[1] = 1 
  nr = 0
  for (i, rn) in enumerate(rns)
    nr += rn
    ix[i + 1] = nr + 1
  end
  return ix
end

function make_Ns(covariates, timevary, mlen, flen)
  Ns = Dict{Symbol, Union{Vector{Int64}, Vector{Vector{Int64}}}}();
  for covar in covariates
    if timevary[covar]
      Ns[covar] = [zeros(Int64, mlen) for i in 1:flen]
    else
      Ns[covar] = zeros(Int64, flen)
    end
  end
  return Ns
end

function gloop!(M, ix, cc, G, mlen, flen, tts, tus)
  
  # for (i, g) in enumerate(G)
  @inbounds Threads.@threads for i = eachindex(tts)
    g = get(
      G,
      (treattime = tts[i], treatunit = tus[i]),
      false
    );
    Ns = make_Ns(cc.covariates, cc.timevary, mlen, flen);
    
    #=
    gfsmin strategy will not work if fs for some match to a treated unit is not a single segment. if it is not in some setting, will need to search.
    # the best way to do that would be, probably, enumerate over gfs and check if each element in gfs is in rfs (as it happens in the loop)

    but, it should form a single segment in all possible cases here, since we
    do not allow missingness between the start and end of a unit's time-series
    =#
    gfsmin = minimum(@views M[ix[i] : ix[i+1] - 1, :f]);
    dog!(M, ix[i], g, cc, mlen, Ns, gfsmin);
    Ndiv!(@view(M[ix[i]:ix[i + 1] - 1, :]), Ns, cc.covariates, cc.fmin);
  end
  return M
end

# @view(M[ix[i]:ix[i + 1] - 1, :])[1, covar] / Ns[covar][1][1]

# function ρdiv!(ρc, Nsc, f, fmin)
#   ρc = ρc ./ Nsc[f - fmin + 1]
#   return ρc
# end

function Ndiv!(Mn, Ns, covariates, fmin)
  for ρ in eachrow(Mn)
  # for φm1 in 0:(ix[i + 1] - ix[i])
    for covar in covariates
      # ρdiv!(ρ[covar], Ns[covar], ρ[:f], fmin);
      ρ[covar] = ρ[covar] ./ Ns[covar][ρ[:f] - fmin + 1]
    end
  end
  return Mn
end

function dog!(M, ixi, g, cc, mlen, Ns, gfsmin)
  # H = deepcopy(M);
  # Ns = make_Ns(cc.covariates, cc.timevary, mlen, flen);
  for r in eachrow(g) # [1:1726]
    # H = deepcopy(cc.meanbalances);
    addmean!(M, ixi, r, r[:fs], cc, mlen, Ns, gfsmin) # one row of M
  end
  #hcat(H[ixi, covar], Ns[covar][1], H[ixi, covar])
  return M
end
  
# #24,25
# rg = eachrow(g)[1725]
# rb = eachrow(g)[1726]

# the Ns are getting stored property (fixed # rows)
# unique f list
# need this to find M row
# this is important when there is f variation across matches to a treated observation
# this is per g
# gfsmin = minimum(@views M[ixi:ix[i+1] - 1, :f]);

function addmean!(M, ixi, r, rfs, cc, mlen, Ns, gfsmin)
  for covar in cc.covariates
    for f in rfs

      #=
      N.B. indexing:

      M over φ, Ns over f - fmin + 1
      M is restricted to elem of r[:fs], Ns is not
      Ns cannot really be, since it is a patchwork of fs present across matches to some treated observation
      =#

      if cc.timevary[covar]
        
        # PROBLEM M+φ is adding to 10, but we actually want 37
        # because other units have that range, we need the total range present
        # consider adding a check (maybe a test later?) that the things match...

        _innersum!(
          M[ixi + f - gfsmin, covar],
          @views(r[covar][fidx(f, mlen, cc.fmin)]),
          Ns[covar][f - cc.fmin + 1]
        )
          # alternatively, store in M
          # M[rownumber(r) + φ - 1, Symbol(string(covar) * "N")],

      else 
        M[ixi + f - gfsmin, covar] += r[covar]
        Ns[covar][f - cc.fmin + 1] += 1
        
        # M[rownumber(r) + φ - 1, Symbol(string(covar) * "N")] += 1
      end
    end
  end

  return M, Ns
end

function _innersum!(mvec, rvec, nvec)
  for l in eachindex(mvec)
    if !ismissing(rvec[l])
      mvec[l] += rvec[l]
      nvec[l] += 1
    end
  end
  return rvec, nvec
end

# rg = eachrow(cc.meanbalances)[(nrow(cc.meanbalances)-36363)]
# rb = eachrow(cc.meanbalances)[(nrow(cc.meanbalances)-36362)]
# rb[covar]
# rbs = eachrow(cc.meanbalances)[(nrow(cc.meanbalances)-36362:nrow(cc.meanbalances)-36350)]

function grandbalance!(cc::AbstractCICModel)

  if cc.stratifier != Symbol("")
    cc.grandbalances = GrandDictStrat();
    us = sort(unique(cc.meanbalances.stratum));
    [cc.grandbalances[s] = GrandDictNoStrat() for s in us];
    _grandbalance!(cc, us)
  else
    cc.grandbalances = GrandDictNoStrat();
    _grandbalance!(cc)
  end

  return cc
end

function _grandbalance!(cc)
  
  for covar in cc.covariates
    if cc.timevary[covar]
      cc.grandbalances[covar] = zeros(Float64, length(cc.mmin:cc.mmax))
    else
      cc.grandbalances[covar] = 0.0
    end
  end

  Ns = make_Ns_m(cc.covariates, cc.timevary, length(cc.mmin:cc.mmax));

  for r in eachrow(cc.meanbalances) # [1:(nrow(cc.meanbalances)-36363)]
    for covar in cc.covariates
      if cc.timevary[covar]
        _inner_sum_nan!(cc.grandbalances[covar], r[covar], Ns[covar])
      else
        cc.grandbalances[covar] += r[covar]
        Ns[covar] += 1
      end
    end  
  end

  # cc.grandbalances[covar]

  for covar in cc.covariates
    cc.grandbalances[covar] = cc.grandbalances[covar] ./ Ns[covar]
  end
  
  return cc
end

function _grandbalance!(cc, us)
  
  for s in us
    for covar in cc.covariates
      if cc.timevary[covar]
        cc.grandbalances[s][covar] = zeros(Float64, length(cc.mmin:cc.mmax))
      else
        cc.grandbalances[s][covar] = 0.0
      end
    end

    Ns = make_Ns_m(cc.covariates, cc.timevary, length(cc.mmin:cc.mmax));

    mbs = @views cc.meanbalances[cc.meanbalances.stratum .== s, :];

    for r in eachrow(mbs) # [1:(nrow(cc.meanbalances)-36363)]
      for covar in cc.covariates
        if cc.timevary[covar]
          _inner_sum_nan!(cc.grandbalances[s][covar], r[covar], Ns[covar])
        else
          cc.grandbalances[s][covar] += r[covar]
          Ns[covar] += 1
        end
      end  
    end

    # cc.grandbalances[covar]

    for covar in cc.covariates
      cc.grandbalances[s][covar] = cc.grandbalances[s][covar] ./ Ns[covar]
    end
  end
  
  return cc
end

function _inner_sum_nan!(mvec, rvec, nvec)
  for i in eachindex(mvec)
    if !isnan(rvec[i])
      mvec[i] += rvec[i]
      nvec[i] += 1
    end
  end
  return rvec, nvec
end

function make_Ns_m(covariates, timevary, mlen)
  Ns = Dict{Symbol, Union{Int64, Vector{Int64}}}();
  for covar in covariates
    if timevary[covar]
      Ns[covar] = zeros(Int64, mlen)
    else
      Ns[covar] = zero(Int64)
    end
  end
  return Ns
end

"""
    balancecheck(m::Dict{Symbol, Union{Float64, Vector{Float64}}}; threshold = 0.1)

Simply check whether the grand means are above the std. balance threshold. Returns a Bool for each covariate.
"""
function balancecheck(
  cc::AbstractCICModel;
  threshold = 0.1
)

  if cc.stratifier == Symbol("")
    chk = Dict{Symbol, Bool}();
    _balancecheck!(chk, v, k, threshold)
  else
    chk = Dict{Int, Dict{Symbol, Bool}}();
    [chk[s] = Dict{Symbol, Bool}() for s in keys(cc.grandbalances)]
    for (k, v) in cc.grandbalances
      for (ki, vi) in v
      _balancecheck!(chk[k], vi, ki, threshold)
      end
    end
  end

  return chk
end

function _balancecheck!(chki, v, k, threshold)
  if any(v .> threshold)
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

#=
## auto balancing

# do we need engage?
function chk_check(chk)
  for (k, v) in chk
    if v
      return x
    end
  end
end

# chk_check(chk)
# stepsize = 1

function docut(c::Vector{Int}, stepsize::Int)
  mx = maximum(c)
  if (mx < 6) | (mx + 1 <= stepsize)
    return mx - stepsize;
  else
    return mx - 1
  end
end

"""
remove highest ranked matches from cc.matches
"""
function threshold!(cc::AbstractCICModel, stepsize)
  
  G = groupby(cc.matches, [:treattime, :treatunit, :f]);
  idx = Int64[];
  for g in G
    rcut = docut(g[!, :rank], stepsize)
    append!(idx, parentindices(g)[1][(findall(g[!, :rank] .> rcut))])
  end

  delete!(cc.matches, idx)

  return cc
end
=#
