function sdtreated(dt, did, cmat, tmin, mlen, dtr)
  #= need to collect all match period indices, across all treated obs
  only for the treated units
  should this only draw from matches?
  what happens when we apply a caliper?

  per <https://cran.r-project.org/web/packages/MatchIt/vignettes/assessing-balance.html>
  we should use the full dataset, and not the matches based on a caliper
  =#

  treatment_points = get_all_treatment_points(dtr);

  minimum(matches.ttime)

  tobst = @view(dt[treatment_points]);
  tobsid = @view(did[treatment_points]);

  # should have 20 elems for each unit in tobs (then x # covar)
  indstore = zeros(Int64, mlen * length(tobst));
  matchtimes = zeros(Int64, mlen * length(tobst));

  indstore = Vector{Union{Missing, Int64}}(missing, mlen * length(tobst));
  matchtimes = similar(indstore);
  iloc = 1

  for i = eachindex(tobst)
    ttime = @views(tobst[i])
    unit = @views(tobsid[i])
    K = (tmin + (tpoint + ttime) + 1) : (tpoint + ttime)
    iind = getvarind(did, dt, unit, K);
    
    ulim = (iloc + length(iind) - 1);

    indstore[iloc : ulim] = iind;
    matchtimes[iloc : ulim] = dt[iind] .- tobst[i];

    iloc += length(iind)
  end

  arethere = findall(ismissing.(indstore) .== false);
  indstore = indstore[arethere];
  matchtimes = matchtimes[arethere];

  # we want a mlen by covariate matrix as output
  tsdmat = zeros(Float64, size(cmat)[2], mlen);

  J = unique(matchtimes);
  cmatn = @view(cmat[indstore, :]);
  
  for (i, j) in enumerate(J)
    mtind = findall(matchtimes .== j);
    tsdmat[:, i] = std(@view(cmatn[mtind, :]), dims = 1)'
  end

  return tsdmat'
end

# this version calculates and stores the distance for each unit
function getbalance(matches, covariates, dat, tmin, tpoint, id, t, treatment)
  # this probably excludes units without matches (even using mm)

  mlen = tpoint - tmin + 1;
  clen = length(covariates);

  cmat = convert(Matrix{Float64}, dat[!, covariates]);
  did = dat[!, id];
  dt = dat[!, t];

  uid = matches[!, :munit];
  utrtid = matches[!, :tunit];
  ut = matches[!, :ttime];

  # -20:-1
  sdcovariates = sdtreated(dt, did, cmat, tmin, mlen, dat[!, treatment]);

  # sold = sdtreatedold(dat, covariates, tmin, id, t, treatment);
  # sdcovariates = Matrix(sold[!, covariates]);

  trtind = findall(uid .== utrtid);

  covariaterepeat = Vector{String}(undef, 0)
  for c in covariates
    append!(covariaterepeat, fill(String(c), length(trtind) * mlen))
  end

  balances = DataFrame(
    score = 0,
    covariate = covariaterepeat);

  balval = zeros(length(trtind) * mlen, length(covariates));
  matchtimes = zeros(Int64, length(trtind) * mlen);
  treatedunits = zeros(Int64, length(trtind) * mlen);
  treatedtimes = zeros(Int64, length(trtind) * mlen);

  balance!(
    balval, matchtimes, treatedunits, treatedtimes,
    did, dt, cmat,
    trtind, utrtid, uid, ut, tmin,
    mlen, tpoint,
    sdcovariates);
    
  balances.score = reshape(balval, length(trtind) * mlen * clen);
  balances.matchtime = repeat(matchtimes, clen)
  balances.treatedunit = repeat(treatedunits, clen)
  balances.treatedtime = repeat(treatedtimes, clen)

  meanbalances = getmeanscore(balances);

  return balances, meanbalances
end

function balance!(
  balval, matchtimes, treatedunits, treatedtimes,
  did::Vector{Int64}, dt::Vector{Int64}, cmat::Matrix{Float64},
  trtind, utrtid, uid, ut, tmin,
  mlen, tpoint,
  sdcovariates)

  # pre-refinement is very slow without parallel
  # for i = eachindex(trtind)
  @inbounds Threads.@threads for i = eachindex(trtind)
    e = trtind[i]
    tunit = utrtid[e]
    tt = ut[e]
    units = uid[findall((utrtid .== tunit) .& (ut .== tt))]

    ws = (units .!= tunit) .* - 1 / (length(units) - 1)
    ws[units .== tunit] .= 1;

    # oc = Vector{Union{Float64, Missing}}(missing, 20);
    oc = zeros(Float64, mlen, size(cmat)[2]);
    
    K = (tmin : tpoint) .+ tt;
    for (ii, u) in enumerate(units)

      iix = getvarind(
        did, dt, u, K)
      if length(iix) != mlen
        println(tunit)
        println(u)
        println(tt)
        println(length(ix))
        error("treated, match, time is not correct length")
      end
      # these indexes are present
      oc .+= @view(cmat[iix, :]) .* ws[ii] # this is time consuming
    end

    lidx = (i - 1) * mlen + 1
    uidx = mlen * i
    balval[lidx:uidx, :] = oc ./ sdcovariates # balance at each time point for treated unit
    matchtimes[lidx:uidx] = tmin:tpoint # possible issue
    treatedunits[lidx:uidx] .= tunit
    treatedtimes[lidx:uidx] .= tt
  end

  return balval, matchtimes, treatedunits, treatedtimes
end

function getmeanscore(balances)
  return @linq meanbalance = balances |>
    groupby([:covariate, :matchtime]) |>
    combine(meanscore = mean(:score));
end

function getbalance_restricted(
  matches_pd,
  stratvar::Symbol,
  covariates, dat, tmin, tpoint, id, t, treatment)

  sv = Symbol(String(stratvar) * "_stratum");
  S = unique(matches_pd[!, sv]);

  balances_post_strat = DataFrame(
    [Float64, String, Int64, Int64, Int64, Int64],
    [:score, :covariate, :matchtime, :treatedunit, :treatedtime, sv]
    );

  meanbalances_post_strat = DataFrame(
    [String, Int64, Float64, Int64],
    [:covariate, :matchtime, :meanscore, sv]
    );

  for s in S
    idx = findall(matches_pd[!, sv] .== s);
    sub = @view(matches_pd[idx, :]);

    balances_post_s, meanbalances_post_s = getbalance(
      sub, covariates, dat, tmin, tpoint, id, t, treatment);

    balances_post_s[!, sv] = ones(Int64, nrow(balances_post_s)) * s
    meanbalances_post_s[!, sv] = ones(Int64, nrow(meanbalances_post_s)) * s

    append!(balances_post_strat, balances_post_s)
    append!(meanbalances_post_strat, meanbalances_post_s)
  end
  return balances_post_strat, meanbalances_post_strat
end
