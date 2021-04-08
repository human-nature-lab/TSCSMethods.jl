
# union type to handle @views
Mtype = Union{
  SubArray{
    Float64,2,
    Array{Float64,2},
    Tuple{Array{Int64,1},
    Base.Slice{Base.OneTo{Int64}}},false},
  Matrix{Float64}}

function attinner!(
  jvals::Vector{Float64},
  om::Mtype,
  wm::Mtype,
  tn::Vector{Int64}
)
  @inbounds for j in 2:size(om)[2]
    @inbounds @simd for i in 1:size(om)[1]
      jvals[j - 1] += @views(om[i, 1]) * @views(wm[i, 1]) + @views(om[i, j]) * @views(wm[i, j])
    end
    jvals[j - 1] = @views(jvals[j - 1]) / tn[j - 1]
  end
  return jvals
end

#= 
  if some treated unit is missing at an f
  we must remove it, and its matches, from the sum (jval[i-1])

  if a match is missing at that f, then we just remove the match,
  although we will need to remove the treated if it is the only
  match left for that treated unit

  want to report the number of treated units at each f when there is
  missingness

  we cannot include a match without its treated unit, and we cannot include matches to a treated unit that doesn't exist

  for an f, skip certain i values (in, say uid[inx]) affected by missingness

  cannot just make the observations zero, b/c of matching-treated comments above
  but, we could make them all 
    
=#
function attinner!(
  jvals::Vector{Float64},
  om::Mtype,
  wm::Mtype,
  tn::Int64
)
  @inbounds for j in 2:size(om)[2]
    @inbounds @simd for i in 1:size(om)[1]
      jvals[j - 1] += @views(om[i, 1]) * @views(wm[i, 1]) + @views(om[i, j]) * @views(wm[i, j])
    end
    jvals[j - 1] = @views(jvals[j - 1]) / tn
  end
  return jvals
end

"""
    att(om, wm, trtnums)

this will need to be adapted for missingness, and we will need a new 
  way to get the treated unit number in the presence of missingness, which will vary across f

"""
function att(om, wm, tn)
  jvals = zeros(Float64, size(om)[2] - 1)
  attinner!(jvals, om, wm, tn)
  return jvals
end

"""
    make_check_sample!(units, clusters, treatedunits, x)

ensure that there are treated units in the resample
"""
function make_check_sample!(units, clusters, treatedunits, x)
  while x < 1
    units = sample(clusters, length(clusters), replace = true);
    x = sum([q in units for q in treatedunits]);
  end;
  return units
end

function attboot(
  iternum,
  fmin, fmax,
  utrtid, uid,
  did,
  omsm, wmsm,
  uidsm, utsm, trtbool, blmat,
  omiss,
  ATT
)

  clusters = unique(did); # clusters comes from the full dataset
  trtidunique = unique(utrtid);
  munique = unique(uid);

  bootests = zeros(length(fmin:fmax), iternum);

  uidsmrep = countmemb(uidsm);

  uidsmtoind = Dict{Int64, Vector{Int64}}();
  makeuiddict!(uidsmtoind, uidsm);

  # treated units left after missing handling, caliper
  trtleft = intersect(utrtid, uidsm);

  # 1 iter: 2.1 sec, 44 MiB
  attboot_inner!(
    bootests,
    clusters,
    omsm, wmsm,
    uidsm, utsm, trtbool, blmat,
    trtidunique, munique,
    uidsmrep, uidsmtoind,
    iternum,
    omiss,
    trtleft,
    ATT
  )

  return bootests
end

function makeuiddict!(uidtoind, uid)
  for u in unique(uid)
    uidtoind[u] = findall(uid .== u)
  end
  return uidtoind
end

function attboot_inner!(
  bootests,
  clusters,
  omsm, wmsm,
  uidsm, utsm, trtbool, blmat,
  trtidunique, munique,
  uidsmrep, uidsmtoind,
  iternum,
  omiss,
  trtleft,
  ATT
)

  @inbounds Threads.@threads for k = eachindex(1:iternum)
  # for k = eachindex(1:iternum)
    units = sample(clusters, length(clusters), replace = true); # 29 μs

    make_check_sample!(units, clusters, trtleft, 0); # 540.9

    tn, inx = treatednum(
      uidsm, utsm,
      units, clusters,
      uidsmrep,
      munique,
      uidsmtoind,
      blmat, trtbool,
      omiss,
      ATT
    );

    while any(tn < 1)
      # should work for missing and vector tn
      # may cause endless loop if there are zero observations for some f in F?
      # maybe use blmat to skip entirely if an f has no treated units
      units = sample(clusters, length(clusters), replace = true); # 29 μs

      tn, inx = treatednum(
        uidsm, utsm,
        units, clusters,
        uidsmrep,
        munique,
        uidsmtoind,
        blmat, trtbool,
        omiss,
        ATT
      );
    end

    bootests[:, k] = att(
      @views(omsm[inx, :]),
      @views(wmsm[inx, :]), tn
    );

  end
  return bootests
end

function getinxlen(uidsmrep, munique, units)
  ur = countmemb(units[(in.(units, Ref(munique)))]);

  inxlen = 0
  for (l, v) in ur
    inxlen += uidsmrep[l] * v
  end
  return inxlen
end

function addidx!(inx, ix, uidtoind, unitreps)
  for (k, v) in unitreps
    
    x = get(uidtoind, k, nothing);
    if !isnothing(x)
      xl = length(x);
      lst = ix + (v * xl) - 1;
      inx[ix : lst] .= xl > 1 ? repeat(x, v) : x
      ix = lst + 1
    end
  end
  return inx
end

function addidx(uidtoind, unitreps, inxlen)
  inx = Vector{Int64}(undef, inxlen); # 18 μs, 11.18 MiB
  ix = 1
  addidx!(inx, ix, uidtoind, unitreps);
  return inx
end

function treatednum(
  uid, ut,
  units, clusters, uidsmrep, munique, uidsmtoind,
  blmat, trtbool,
  omiss::Bool,
  ATT::Bool
)

  unitreps = countmemb(units, length(clusters));

  inxlen = getinxlen(uidsmrep, munique, units);
  inx = addidx(uidsmtoind, unitreps, inxlen);

  uinx = @views(uid[inx]); utinx = @views(ut[inx]);
  
  if ATT == true
    if omiss == false
      tn = sum(@views(trtbool[inx]))
    elseif omiss == true
      tcnt = sum(@views(trtbool[inx]));
      tn = tcnt .- sum(@views(blmat[inx, :]), dims = 1)
      tn = reshape(tn, length(tn))
    end
  else
    tn = 1
  end
  
  return tn, inx
end

function getresamplelen(unitreps, ulens)
  relen = Int64(0)
  for (key, val) in unitreps
    relen += ulens[key] * val
  end
  return relen
end