
function attinner!(jvals::Vector{Float64}, om, wm, trtnums)
  @inbounds for j in 2:size(om)[2]
    @inbounds @simd for i in 1:size(om)[1]
      jvals[j - 1] += @views(om[i, 1]) * @views(wm[i, 1]) + @views(om[i, j]) * @views(wm[i, j])
    end
    jvals[j - 1] = @views(jvals[j - 1]) / trtnums
  end
  return jvals
end

"""
    att(om, wm, trtnums)

this will need to be adapted for missingness, and we will need a new 
  way to get the treated unit number in the presence of missingness, which will vary across f

"""
function att(om, wm, trtnums)
  jvals = zeros(Float64, size(om)[2] - 1)
  attinner!(jvals, om, wm, trtnums)
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
  outcomemat, dwitmat,
  )

  clusters = unique(did); # clusters comes from the full dataset
  trtidunique = unique(utrtid);
  munique = unique(uid);

  bootests = zeros(length(Fset), iternum);

  uidrep = countmemb(uid);
  
  uidtoind = Dict{Int64, Vector{Int64}}();
  makeuiddict!(uidtoind, uid);

  # 1 iter: 2.1 sec, 44 MiB
  attboot_inner!(
    bootests,
    clusters, outcomemat, dwitmat,
    uid, utrtid, trtidunique, munique,
    uidrep, uidtoind,
    iternum
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
  outcomemat, dwitmat,
  uid, utrtid, trtidunique, munique,
  uidrep, uidtoind,
  iternum
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

    unitreps = countmemb(units, length(clusters)); # 51.8 μs, 68.9 KiB
    
    inxlen = getinxlen(uidrep, munique, units); # 1.5 ms; 120.3 KiB
    inx = addidx(uidtoind, unitreps, inxlen); # 2.7 ms, 22.79 MiB -- ridiculous

    tn = sum(uid[inx] .== utrtid[inx]);

    #= will need way to handle missingness in treated units at some f points
      do this by getting the total number of appearances of a unit, and 
      subtracting that many times for a missing in outcomemat
      sum(ismissing.(outcomemat), dims = 2)
    =#

    #= with att(): ~210 ms (down from 1 sec), and ~500 bytes (from 181 GiB!!)
      this is faster, but still need missingness solution, and a
      way to deal with treatnums in case of missingness
    =#

    bootests[:, k] = att( # 2.304 sec, 464 bytes
      @views(outcomemat[inx, :]),
      @views(dwitmat[inx, :]), tn
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

function getinxlen(uidrep, munique, units)
  ur = countmemb(units[(in.(units, Ref(munique)))]);

  inxlen = 0
  for (l, v) in ur
    inxlen += uidrep[l] * v
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

"""
    addindices!(indices, witsf_id, unitunique, unitreps, cnt)

function barrier
add the actual index values to grab from wits, with reptitions
"""
function addindices!(indices, witsf_id, unitreps, cnt)
  for (k, v) in unitreps
    toadd = repeat(findall(witsf_id .== k), v)
    tal = length(toadd)
    indices[(cnt + 1) : (cnt + tal)] = toadd
    cnt += tal
  end
  return indices
end


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