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
  uid, utrtid, ut, mcnts, nd, tpoint, tl,
  fmin, fmax;
  summarize = true
)
  outcomemat = Matrix{Union{Float64, Missing}}(missing, length(uid), tl);
  dwitmat = similar(outcomemat);
  makemats!(
    outcomemat, dwitmat, uid, utrtid, ut,
    mcnts, nd, tpoint,
    fmin, fmax
  );

  # type handling for att() (in bootstrap_estimation.jl)
  omiss = false;
  try
    outcomemat = convert(Array{Float64, 2}, outcomemat);
    dwitmat = convert(Array{Float64, 2}, dwitmat);
  catch
    omiss = true
    println("missing post-treatment outcome observations present")
  end

  if (typeof(outcomemat) == Array{Union{Missing, Float64},2}) & (typeof(dwitmat) == Array{Union{Missing, Float64},2})
    # om = deepcopy(outcomemat); wm = deepcopy(dwitmat);
    outcomemat, dwitmat, bl = missingmats(outcomemat, dwitmat, uid, utrtid, ut);
  else bl = DataFrame();
  end

  if summarize == true
    omsm, wmsm, uidsm, utsm, trtbool, blmat = summarizemats(
      uid, ut, outcomemat, dwitmat, bl, keys(mcnts));
    return omsm, wmsm, uidsm, utsm, trtbool, blmat, omiss
  else
    return outcomemat, dwitmat, bl, omiss
  end
end

#=
estimation without restricted averaging
=#
function standard_estimation(
  iternum::Int64, matches5::DataFrame,
  Fset, tpoint,
  dat::DataFrame,
  id::Symbol, t::Symbol, outcome::Symbol,
  ATT::Bool,
  outboot::Bool
)
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

  omsm, wmsm, uidsm, utsm, trtbool, blmat, omiss = handlemats(
    uid, utrtid, ut, mcnts, nd, tpoint, tl,
    minimum(Fset), maximum(Fset)
  )
  
  println("mats have been made")

  bootests = attboot(
    iternum, Fset,
    utrtid, uid,
    did,
    omsm, wmsm,
    uidsm, utsm, trtbool, blmat,
    omiss,
    ATT
  );

  println("boots have been strapped")

  if ATT == true
    tn = att_trtnums(blmat, trtbool, omiss);
  else
    tn = 1
  end

  ests = att(omsm, wmsm, tn);
    
  bootests = permutedims(bootests)

  results = outputprocess(bootests, ests, Fset);

  if outboot == true
    return results, bootests
  elseif outboot == false
    return results
  end
end

# for use outside the bootstrap
function att_trtnums(blmat, trtbool, omiss)
  if omiss == false
    tn = sum(trtbool)
  else
    tcnt = sum(@views(trtbool));
    tn = tcnt .- sum(@views(blmat), dims = 1)
    tn = reshape(tn, length(tn))
  end
  return tn
end

# input pre-constructed integer vect for the group splits
function restricted_estimation(
  iternum, m, Fset, tpoint,
  dat, stratvar::Symbol,
  id, t, outcome,
  ATT::Bool,
  outboot::Bool)

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
  
  om, wm, bl, omiss = handlemats(
    uid, utrtid, ut, mcnts, nd, tpoint, tl,
    minimum(Fset), maximum(Fset);
    summarize = false
  )

  println("mats have been made")

  Results = DataFrame(
    [Int64, Int64, Float64, Float64, Float64, Float64],
    [:stratum, :f, :att, :lwer, :med, :uper]
  )

  Sbootests = Dict{Int64, Matrix{Float64}}();

  restricted_estimation_inner!(
    Results, om, wm,
    utrtid, uid, ut,
    bl, mcnts,
    Fset, did, iternum,
    S, us,
    omiss,
    ATT,
    outboot,
    Sbootests
  )
  
  if outboot == false
    return Results
  elseif outboot == true
    return Results, Sbootests
  end

end

function restricted_estimation_inner!(
  Results, om, wm,
  utrtid, uid, ut,
  bl, mcnts,
  Fset, did, iternum,
  S, us,
  omiss,
  ATT,
  outboot,
  Sbootests,
)

  @inbounds Threads.@threads for i = eachindex(S)
  # for s in S
    s = S[i];

    sind = findall(us .== s);

    utrtid_s = @view(utrtid[sind]);
    uid_s = @view(uid[sind]);
    ut_s = @view(ut[sind]);

    om_s = @view(om[sind, :]);
    wm_s = @view(wm[sind, :]);

    # want to select bl vals that exist in s
    if size(bl)[1] > 0
      blr = Int64[]
      for i = eachindex(utrtid_s)
        utrti = utrtid_s[i]
        uti = utrtid_s[i]

        for j = eachindex(1:nrow(bl))
          if (bl.id == utrti) & (bl.t == uti)
            push!(blr, j)
          end
        end
      end
      bl_s = @views(bl[blr, :])
    else
      bl_s = zeros(Int64, 0, 0);
    end

    omsm, wmsm, uidsm, utsm, trtbool_s, blmat = summarizemats(
      uid_s, ut_s, om_s, wm_s, bl_s, keys(mcnts));

    bootests_s = attboot(
      iternum, Fset,
      utrtid_s, uid_s,
      did,
      omsm,
      wmsm,
      uidsm, utsm, trtbool_s, blmat,
      omiss,
      ATT
    );

    println("boots have been strapped")

    tn_s = sum(utrtid_s .== uid_s);
    ests_s = att(omsm, wmsm, tn_s);
      
    bootests_s = permutedims(bootests_s)

    if outboot == true
      Sbootests[s] = bootests_s;
    end

    results = outputprocess(bootests_s, ests_s, Fset);

    results.stratum = s * ones(Int64, nrow(results));

    append!(Results, results)
  end
  
  if outboot == true
    return Results, Sbootests
  elseif outboot == false
    return Results
  end

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