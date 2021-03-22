
# union type of handle @views
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
  iternum, Fset,
  utrtid, uid,
  did,
  omsm, wmsm,
  uidsm, utsm, trtbool, blmat,
  omiss
  )

  clusters = unique(did); # clusters comes from the full dataset
  trtidunique = unique(utrtid);
  munique = unique(uid);

  bootests = zeros(length(Fset), iternum);

  # uidrep = countmemb(uid);
  uidsmrep = countmemb(uidsm);
  
  # uidtoind = Dict{Int64, Vector{Int64}}();
  # makeuiddict!(uidtoind, uid);

  uidsmtoind = Dict{Int64, Vector{Int64}}();
  makeuiddict!(uidsmtoind, uidsm);

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
  omiss
  )

  @inbounds Threads.@threads for k = eachindex(1:iternum)
  # for k = eachindex(1:iternum)
    units = sample(clusters, length(clusters), replace = true); # 29 μs

    make_check_sample!(units, clusters, trtidunique, 0); # 540.9
      # kinda slow

    #=
      this gives the length of the resampled indices
      product of the number of times the unit appears in the data and 
      number of repeats in the bootstrapped sample
      think about it this way
      we are really trying to resample the data,
      and resample units in the data
      those then get assigned as matches with dwit values
      if a unit is resampled > 1x, then it should be reassigned that many times
    =#

    tn, inx = treatednum(
      uidsm, utsm,
      units, clusters,
      uidsmrep,
      munique,
      uidsmtoind,
      blmat, trtbool,
      omiss); # omiss

    bootests[:, k] = att( # 2.304 sec, 464 bytes
      @views(omsm[inx, :]),
      @views(wmsm[inx, :]), tn
    );

  end
  return bootests
end

# """
#     addindices!(indices, witsf_id, unitunique, unitreps, cnt)
# function barrier
# add the actual index values to grab from wits, with reptitions
# """
# function addindices!(indices, witsf_id, unitunique, unitreps, cnt)
#   for i = eachindex(unitunique)
#     ui = unitunique[i]
#     toadd = repeat(findall(witsf_id .== ui), unitreps[ui])
#     tal = length(toadd)
#     indices[(cnt + 1) : (cnt + tal)] = toadd
#     cnt += tal
#   end
#   return indices
# end

function getinxlen(uidsmrep, munique, units)
  ur = countmemb(units[(in.(units, Ref(munique)))]);

  inxlen = 0
  cnt = 0
  for (l, v) in ur
    cnt += 1
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
  omiss::Bool # used as type check for method
)

  unitreps = countmemb(units, length(clusters)); # 51.8 μs, 68.9 KiB

  inxlen = getinxlen(uidsmrep, munique, units); # 1.5 ms; 120.3 KiB
  inx = addidx(uidsmtoind, unitreps, inxlen); # 2.7 ms, 22.79 MiB -- ridiculous

  uinx = @views(uid[inx]); utinx = @views(ut[inx]);
  
  if omiss == false
    tn = sum(@views(trtbool[inx]))
  elseif omiss == true
    tcnt = sum(@views(trtbool[inx]));
    tn = tcnt .- sum(@views(blmat[inx, :]), dims = 1)
    tn = reshape(tn, length(tn)), inx
  end
  
  return tn, inx
end

# """
#     treatednum(
#       uid,
#       utrtid,
#       units, clusters, uidrep, munique,
#       outcomemat::Array{Union{Missing, Float64}, 2}
#     )

# return the number of treated units (or at each f in F, in case of missingness)
# and the indices needed to construct the bootstrap sample.
# """
# function treatednum(
#   uid,
#   utrtid,
#   units, clusters, uidrep, munique, uidtoind,
#   outcomemat::Array{Union{Missing, Float64}, 2}
# )

#   unitreps = countmemb(units, length(clusters)); # 51.8 μs, 68.9 KiB

#   inxlen = getinxlen(uidrep, munique, units); # 1.5 ms; 120.3 KiB
#   inx = addidx(uidtoind, unitreps, inxlen); # 2.7 ms, 22.79 MiB -- ridiculous

#   tl = uid[inx] .== utrtid[inx];
#   tn = sum(tl);

#   tns = Vector{Int64}(undef, size(outcomemat)[2] - 1);
#   treatednum!(tn, tns, tl, inx, om)
  
#   return tns, inx
# end

# """
#     treatednum!(
#       tn,
#       tns,
#       tl,
#       inx,
#       om::Array{Union{Missing, Float64}, 2}
#     )

# In case of missingness in outcomes, get the number of treated units
# in the bootstrap (or original) sample present at each point in the F range.
# """
# function treatednum!(
#   tn,
#   tns,
#   tl,
#   inx,
#   om::Array{Union{Missing, Float64}, 2}
# )

#   @inbounds for j in 2:size(outcomemat)[2]
#     tns[j - 1] = tn - sum(ismissing.(@view(om[inx[tl], j]))) # this should still be correct, provided that we remove treated units without matches left from missingness. 
#     # this won't be correct if we use summarized values, then we need to know which ones to remove from list

#   end
#   return tns
# end

# """
#     addindices!(indices, witsf_id, unitunique, unitreps, cnt)

# function barrier
# add the actual index values to grab from wits, with reptitions
# """
# function addindices!(indices, witsf_id, unitreps, cnt)
#   for (k, v) in unitreps
#     toadd = repeat(findall(witsf_id .== k), v)
#     tal = length(toadd)
#     indices[(cnt + 1) : (cnt + tal)] = toadd
#     cnt += tal
#   end
#   return indices
# end


function countmemb(itr)
  d = Dict{eltype(itr), Int}()
  for val in itr
      # we don't need these checks
      #if isa(val, Number) && isnan(val)
      #    continue
      #end
      d[val] = get(d, val, 0) + 1
  end
  return d
end

function countmemb(itr, len::Int64)
  d = Dict{eltype(itr), Int64}()
  sizehint!(d, len)
  for val in itr
      # we don't need these checks
      #if isa(val, Number) && isnan(val)
      #    continue
      #end
      d[val] = get(d, val, 0) + 1
  end
  return d
end

function getresamplelen(unitreps, ulens)
  relen = Int64(0)
  for (key, val) in unitreps
    relen += ulens[key] * val
  end
  return relen
end