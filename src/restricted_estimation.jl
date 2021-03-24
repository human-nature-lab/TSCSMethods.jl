# restricted estimation rework

# currently not working

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
  
  om, wm, bl, omiss = handlemats(
    uid, utrtid, ut, mcnts, nd, tpoint, tl;
    summarize = false
  )

  println("mats have been made")

  Results = DataFrame(
    [Int64, Int64, Float64, Float64, Float64, Float64],
    [:stratum, :f, :att, :lwer, :med, :uper]
  )

  restricted_estimation_inner!(
    Results, om, wm,
    utrtid, uid, ut,
    bl, mcnts,
    Fset, did, iternum,
    S, us,
  )
  
  return Results
end

function restricted_estimation_inner!(
  Results, om, wm,
  utrtid, uid, ut,
  bl, mcnts,
  Fset, did, iternum,
  S, us
)

  @inbounds Threads.@threads for i = eachindex(S)

  Sind = Dict{Int64, Vector{Int64}}();
  Soms = Dict{Int64, Matrix{Float64}}();
  Swms = Dict{Int64, Matrix{Float64}}();
  Suidsm = Dict{Int64, Vector{Int64}}();
  Sutsm = Dict{Int64, Vector{Int64}}();
  Strtbool = Dict{Int64, Vector{Int64}}();
  Sblmat = Dict{Int64, Matrix{Float64}}();

  Suid = Dict{Int64, Vector{Int64}}();
  Sutrtid = Dict{Int64, Vector{Int64}}();
  
  for s in S
    Sind[s] = findall(us .== s);
    sind = Sind[s];

    Sutrtid[s] = @view(utrtid[sind]);
    Suid[s] = @view(uid[sind]);
    ut_s = @view(ut[sind]);

    om_s = @view(om[sind, :]);
    wm_s = @view(wm[sind, :]);

    # want to select bl vals that exist in s
    if size(bl)[1] > 0
      blr = Int64[]
      for i = eachindex(utrtid_s)
        utrti = Sutrtid[s][i]
        uti = Sutrtid[s][i]

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

    Soms[s], Swms[s], Suidsm[s], Sutsm[s], Strtbool[s], Sblmat[s] = summarizemats(
      Suid[s], ut_s, om_s, wm_s, bl_s, keys(mcnts));

  end

  bootests_s = attboot(
    iternum, Fset,
    utrtid_s, uid_s,
    did,
    omsm,
    wmsm,
    uidsm, utsm, trtbool_s, blmat,
    omiss
    );

  println("boots have been strapped")

  tn_s = sum(utrtid_s .== uid_s);
  ests_s = att(omsm, wmsm, tn_s);
    
  bootests_s = bootests_s'

  results = outputprocess(bootests_s, ests_s, Fset);

  results.stratum = s * ones(Int64, nrow(results));

  append!(Results, results)

  return Results
end

## bootstrap estimation

# restricted-specific

function attboot_res(
  iternum, Fset,
  utrtid, uid,
  did,
  omsm, wmsm,
  uidsm, utsm, trtbool, blmat,
  omiss
)

  clusters = unique(did); # clusters comes from the full dataset
  
  Strtidunique = Dict{Int64, Vector{Int64}}();
  [Strtidunique[s] = unique(Sutrtid[s]) for s in S];
  
  Smunique = Dict{Int64, Vector{Int64}}();
  [Smunique[s] = unique(Suid[s]) for s in S];

  Sbootests = Dict{Int64, Matrix{Float64}}();
  [Sbootests[s] = zeros(length(Fset), iternum) for s in S];

  Suidsmrep = Dict{Int64, Dict{Int64, Int64}}();
  [Suidsmrep[s] = countmemb(Suidsm[s]) for s in S];

  Suidsmtoind = Dict{Int64, Dict{Int64, Vector{Int64}}}()
  uidsmtoind = Dict{Int64, Vector{Int64}}();
  [Suidsmtoind[s] = makeuiddict!(uidsmtoind, Suidsm[s]) for s in S];

  # 1 iter: 2.1 sec, 44 MiB
  attboot_inner!(
    bootests,
    clusters,
    omsm, wmsm,
    uidsm, utsm, trtbool, blmat,
    trtidunique, munique,
    uidsmrep, uidsmtoind,
    iternum,
    omiss
  )

  return bootests
end

"""
ensure that treated unit from each subgroup exists in bootstrap sample
"""
function make_check_sample_S!(units, clusters, Streatedunits, S)
  x = zeros(Int64, length(S))
  while any(x .== 0)
    units = sample(clusters, length(clusters), replace = true);
    
    [x[i] = sum([q in units for q in Streatedunits[s]]) for (s, i) in enumerate(S)]
  end;
  return units
end

function attboot_inner_res!(
  bootests,
  clusters,
  omsm, wmsm,
  uidsm, utsm, trtbool, blmat,
  trtidunique, munique,
  uidsmrep, uidsmtoind,
  iternum,
  omiss
)

  @inbounds Threads.@threads for k = eachindex(1:iternum)
  # for k = eachindex(1:iternum)
    
    units = sample(clusters, length(clusters), replace = true); # 29 μs
    make_check_sample_S!(units, clusters, Strtidunique, S)
    

    for s in S

      tn, inx = treatednum(
        Suidsm[s], Sutsm[s],
        units, clusters,
        Suidsmrep[s],
        Smunique[s],
        Suidsmtoind[s],
        Sblmat[s], Strtbool[s],
        omiss
      );

      Sbootests[s][:, k] = att(
        @views(Soms[s][inx, :]),
        @views(Swms[s][inx, :]), tn
      );
    end
  end
  
  return Sbootests
end