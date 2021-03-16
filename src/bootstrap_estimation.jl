# maybe add line to remove entries that have pretreatment missingness

function newattcalc(dwitmat, outcomemat, Fset, trtnums::Vector{Int64})
  newatts = zeros(Float64, length(Fset)); # consider missingness
  for i = 2 : (length(Fset) + 1)
    dwit1i = dwitmat[[1, i], :];
      #= extend missingness to pretreatment dwits so that they are
         excluded from the sum =#
    dwit1i[[1], findall(ismissing.(@view(dwit1i[2, :])))] .= missing;
    out1i = @view(outcomemat[[1,i], :]);
    newatts[i - 1] = sum(skipmissing(dwit1i .* out1i)) /
      @view(trtnums[i - 1])[1]
  end
  return newatts
end

function newattcalc(dwitmat, outcomemat, Fset, trtnums::Int64)
  newatts = zeros(Float64, length(Fset)); # consider missingness
  for i = 2 : (length(Fset) + 1)
    dwit1i = dwitmat[[1, i], :];
      #= extend missingness to pretreatment dwits so that they are
         excluded from the sum =#
    dwit1i[[1], findall(ismissing.(@view(dwit1i[2, :])))] .= missing;
    out1i = @view(outcomemat[[1, i], :]);
    newatts[i - 1] = sum(skipmissing(dwit1i .* out1i)) / trtnums
  end
  return newatts
end

function attinner!(jvals, om, wm, trtnums::Vector{Int64})
  @inbounds for j in 2:size(om)[2]
    @inbounds @simd for i in 1:size(om)[1]
      jvals[j - 1] += @views(om[i, 1]) * @views(wm[i, 1]) + @views(om[i, j]) * @views(wm[i, j])
    end
    jvals[j - 1] = @views(jvals[j - 1]) / trtnums[j - 1]
  end
  return jvals
end

function attinner!(jvals, om, wm, trtnums::Float64)
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

** currently differs in result from newattcalc **
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
    x = sum([q in units for q in treatedunits])
  end
  return units
end

"""
  gettrtnum()

gets the number of treated observations in the data or bootstrap sample

treated units without a single matched unit should be removed during
the matching procedure (check this again)

this needs to output a vector for each post treatment timepoint,
in case of missingness
"""
function gettrtnums(indices, uid, utrtid, outcomemat)
  to = findall(@view(uid[indices]) .== @view(utrtid[indices]));
  tolen = length(to);

  trtnums = repeat([tolen], size(outcomemat)[1]);

  for i = eachindex(trtnums)
    missinglocsi = sum(ismissing.(@view(outcomemat[i, indices])));
    trtnums[i] -= length(intersect(to, missinglocsi));
      # number of missing treated units
  end

  return trtnums
end

function gettrtnums(uid, utrtid, outcomemat)
  to = findall(uid .== utrtid);
  tolen = length(to);

  trtnums = repeat([tolen], size(outcomemat)[1]);

  for i = eachindex(trtnums)
    missinglocsi = sum(ismissing.(@view(outcomemat[i, :])));
    trtnums[i] -= length(intersect(to, missinglocsi));
      # number of missing treated units
  end

  return trtnums
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

  # no huge gain from parallelizing this function
  # 131 to 45 sec -- why? instabilities?

  uidrep = countmemb(uid);
  
  uidtoind = Dict{Int64, Vector{Int64}}();
  makeuiddict!(uidtoind, uid);

  # @btime yields 38.078 s (53038237 allocations: 7.18 GiB)
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
    # clusters = dat_id_un;
    units = sample(clusters, length(clusters), replace = true); # 29 μs

    make_check_sample!(units, clusters, trtidunique, 0); # 521.1 μs

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

    # begin old way
    # unitreps = countmap(units, alg = :auto); # 74 μs
    # unitunique = unique(units); # 65 μs

    # # # preallocate index list (just make 2x as easy solution)
    # indices = zeros(Int64, length(uid) * 2); # bad for possible estimate
    # addindices!(indices, uid, unitunique, unitreps, 0); # 5.9 ms -- !!
    # indices = indices[indices .> 0]; # 9.3 μs (up to 11.2 with easy sol)

    # trtnumsold = gettrtnums(indices, uid, utrtid, outcomemat');

    
    # end old way
    
    # new mod way
    
    unitreps = countmap(units, alg = :auto); # 74 μs
    
    inxlen = getinxlen(uidrep, munique, units); # 1.5 ms
    inx = addidx(uidtoind, unitreps, inxlen); # 153.7 μs
    
    # end new mod way
    
    # trtnums = gettrtnums(inx, uid, utrtid, outcomemat'); # 173 ms -> 30 minutes
    
    # sort(inx) == sort(indices)

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

    # bt = newattcalc(dwitmat[inx, :]', outcomemat[inx, :]', Fset, trtnums)

    # cmp = hcat(bt, bootests[:, k])

    bootests[:, k] = att(
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
  inx = Vector{Int64}(undef, inxlen); # 312.8 μs
  #inx = zeros(Int64, inxlen); # 312.8 μs
  ix = 1
  addidx!(inx, ix, uidtoind, unitreps)
  return inx
end

function repeatnew!(tofill, v, x)
  st = 1
  @inbounds for i = eachindex(x)
    en = st + v - 1
    tofill[st : en] .= x[i]
    st = en + 1
  end
  return tofill
end

function repeatnew(x, v)
  tofill = Vector{Int64}(undef, length(x) * v)
  repeatnew!(tofill, v, x)
  return tofill
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

function getresamplelen(unitreps, ulens)
  relen = Int64(0)
  for (key, val) in unitreps
    relen += ulens[key] * val
  end
  return relen
end