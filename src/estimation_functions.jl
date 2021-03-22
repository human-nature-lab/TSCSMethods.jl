# estimation_functions.jl

"""
    outputprocess(bootests, ests, Fset)

process the bootstrap and estimation output into a dataframe
for easy access and plotting
"""
function outputprocess(booteststr, ests, Fset; ptiles = [0.025, 0.5, 0.975])
  lwr = zeros(Float64, length(ests));
  med = zeros(Float64, length(ests));
  upr = zeros(Float64, length(ests));

  for i = eachindex(ests)
    lwr[i], med[i], upr[i] = quantile(@view(booteststr[:, i]), ptiles);
  end

  results = DataFrame(
    f = Fset,
    att = ests,
    lwer = lwr,
    med = med,
    uper = upr
  );

  return results
end

DDict = Dict{Tuple{Int64, Int64}, Float64};

function outcomedict(did, dt, dout)
  nd = DDict();
  for i = eachindex(1:length(dt))
    nd[(did[i], dt[i])] = dout[i]
  end;
  return nd
end

function handlemats(
  uid, utrtid, ut, mcnts, nd, tpoint, tl
)
  outcomemat = Matrix{Union{Float64, Missing}}(missing, length(uid), tl);
  dwitmat = similar(outcomemat);
  makemats!(outcomemat, dwitmat, uid, utrtid, ut, mcnts, nd, tpoint);

  # type handling for att() (in bootstrap_estimation.jl)
  try
    outcomemat = convert(Array{Float64, 2}, outcomemat);
    dwitmat = convert(Array{Float64, 2}, dwitmat);
    omiss = false;
  catch
    omiss = true
    println("missing post-treatment outcome observations present")
  end

  if (typeof(outcomemat) == Array{Union{Missing, Float64},2}) & (typeof(dwitmat) == Array{Union{Missing, Float64},2})
    # om = deepcopy(outcomemat); wm = deepcopy(dwitmat);
    outcomemat, dwitmat, bl = missingmats(outcomemat, dwitmat, uid, utrtid);
  else bl = DataFrame();
  end

  omsm, wmsm, uidsm, utsm, trtbool, blmat = summarizemats(
    uid, ut, outcomemat, dwitmat, bl, keys(mcnts));
  return omsm, wmsm, uidsm, utsm, trtbool, blmat
end

#=
estimation without restricted averaging
=#
function standard_estimation(
  iternum::Int64, matches5::DataFrame,
  Fset, tpoint,
  dat::DataFrame,
  id::Symbol, t::Symbol, outcome::Symbol)
  # order makes a HUGE time difference
  sort!(matches5, [:ttime, :tunit, :possible, :munit]);

  utrtid = matches5[!, :tunit];
  uid = matches5[!, :munit];
  ut = matches5[!, :ttime];

  did = dat[!, id];
  dt = dat[!, t];
  doc = dat[!, outcome];

  mcnts = matchcounts(matches5);
  nd = outcomedict(did, dt, doc);
  
  tl = length(Fset) + 1;

  omsm, wmsm, uidsm, utsm, trtbool, blmat = handlemats(
    uid, utrtid, ut, mcnts, nd, tpoint, tl
  )
  
  println("mats have been made")

  bootests = attboot(
    iternum, Fset,
    utrtid, uid,
    did,
    omsm, wmsm,
    uidsm, utsm, trtbool, blmat,
    omiss
  );

  println("boots have been strapped")

  ests = att(omsm, wmsm, sum(trtbool));
    
  bootests = bootests'

  results = outputprocess(bootests, ests, Fset);

  return results
end

# input pre-constructed integer vect for the group splits
function restricted_estimation(
  iternum, m, Fset, tpoint,
  dat, stratvar::Symbol,
  id, t, outcome)

  stratname = Symbol(String(stratvar) * "_stratum");

  # must not include missing
  if eltype(m[!, stratname]) != Int64
    try
      m[!, stratname] = convert.(Int64, m[!, stratname]);
    catch e
      println(
        "there are missing or non-integer values in stratification variable"
      )
    end
  end

  mm = sort(m, [stratname, :ttime, :tunit, :munit]);
  
  did = dat[!, id];
  dt = dat[!, t];
  doc = dat[!, outcome];

  utrtid = mm[!, :tunit];
  uid = mm[!, :munit];
  ut = mm[!, :ttime];
  
  us = mm[!, stratname];
  S = unique(us);

  mcnts = matchcounts(mm);
  nd = outcomedict(did, dt, doc);
  
  tl = length(Fset) + 1;
  outcomemat = Matrix{Union{Float64, Missing}}(missing, length(uid), tl);
  dwitmat = similar(outcomemat);
  makemats!(outcomemat, dwitmat, uid, utrtid, ut, mcnts, nd, tpoint);

  println("mats have been made")

  Results = DataFrame(
    [Int64, Int64, Float64, Float64, Float64, Float64],
    [:stratum, :f, :att, :lwer, :med, :uper]
  )

  restricted_estimation_inner!(
    Results, outcomemat, dwitmat,
    utrtid, uid, ut,
    Fset, did, iternum,
    S, us
  )
  
  return Results
end

function restricted_estimation_inner!(
  Results, outcomemat, dwitmat,
  utrtid, uid, ut, Fset, did, iternum,
  S, us
  )

  @inbounds Threads.@threads for i = eachindex(S)
  # for s in S
    s = S[i];

    sind = findall(us .== s);

    utrtid_s = @view(utrtid[sind]);
    uid_s = @view(uid[sind]);
    ut_s = @view(ut[sind]);

    outcomemat_s = @view(outcomemat[sind, :])
    dwitmat_s = @view(dwitmat[sind, :])

    bootests_s = attboot(
      iternum, Fset,
      utrtid_s, uid_s,
      did,
      outcomemat_s,
      dwitmat_s
      );

    println("boots have been strapped")

    # otrtnums = gettrtnums(
    #   uid_s,
    #   utrtid_s,
    #   outcomemat_s);
    # ests_s = newattcalc(
    #   dwitmat_s,
    #   outcomemat_s,
    #   Fset, otrtnums);

    tn_s = sum(utrtid_s .== uid_s);
    ests_s = att(outcomemat_s, dwitmat_s, tn_s);
      
    bootests_s = bootests_s'

    results = outputprocess(bootests_s, ests_s, Fset);

    results.stratum = s * ones(Int64, nrow(results));

    append!(Results, results)
  end
  return Results
end

"""
    assignquantile(stratvar, id, fulldat)

takes a dataset with the relevant variables and outputs a df
of quantile assignments for the chosen variable, unit-by-unit

this is used to create a variable with strata for restricted averaging in
att estimation
"""
function assignquantile(
  stratvar::Symbol, id::Symbol, fulldat::DataFrame)
  ref = unique(
    fulldat[!, [id, stratvar]]);
  sq = quantile(ref[!, stratvar]);

  stratname = Symbol(String(stratvar) * "_stratum");
  ref[!, stratname] = zeros(Int64, nrow(ref));

  for i = eachindex(ref[!, stratname])
    sv = @view(ref[!, stratvar][i])[1]
    if sv >= sq[4]
      ref[!, stratname][i] = 4
    elseif (sv >= sq[3]) & (sv < sq[4])
      ref[!, stratname][i] = 3
    elseif (sv >= sq[2]) & (sv < sq[3])
      ref[!, stratname][i] = 2
    elseif (sv >= sq[1]) & (sv < sq[2])
      ref[!, stratname][i] = 1
    end
  end
  return select(ref, [id, stratname]);
end

#= old stuff
=#

# this simply tacks on the treated observations themselves as rows
# maybe want a better way to handle this
function add_trtobsrows(mm5)
  M = vcat(
    mm5,
    DataFrame(
      trt_t = unique(mm5[[:trt_t, :trt_fips]]).trt_t,
      trt_fips = unique(mm5[[:trt_t, :trt_fips]]).trt_fips,
      match_fips = unique(mm5[[:trt_t, :trt_fips]]).trt_fips,
      mdist = 0.0)
      );
  sort!(M, [:trt_t, :trt_fips]);
  return M
end

