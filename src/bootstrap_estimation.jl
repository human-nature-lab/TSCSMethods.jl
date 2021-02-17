# maybe add line to remove entries that have pretreatment missingness

function newattcalc(dwitmat, outcomemat, Fset, trtnums::Vector{Int64})
  newatts = zeros(Float64, length(Fset)); # consider missingness
  for i = 2:(length(Fset)) + 1
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
  for i = 2:(length(Fset)) + 1
    dwit1i = dwitmat[[1, i], :];
      #= extend missingness to pretreatment dwits so that they are
         excluded from the sum =#
    dwit1i[[1], findall(ismissing.(@view(dwit1i[2, :])))] .= missing;
    out1i = @view(outcomemat[[1,i], :]);
    newatts[i - 1] = sum(skipmissing(dwit1i .* out1i)) / trtnums
  end
  return newatts
end

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

  bootests = zeros(length(Fset), iternum);

  # no huge gain from parallelizing this function
  # 131 to 45 sec -- why? instabilities?
  attboot_inner!(
    bootests,
    clusters, outcomemat, dwitmat,
    uid, utrtid, trtidunique, Fset, iternum
  )

  return bootests
end

"""
    addindices!(indices, witsf_id, unitunique, unitreps, cnt)

function barrier
add the actual index values to grab from wits, with reptitions
"""
function addindices!(indices, witsf_id, unitunique, unitreps, cnt)
  for i = eachindex(unitunique)
    ui = unitunique[i]
    toadd = repeat(findall(witsf_id .== ui), unitreps[ui])
    tal = length(toadd)
    indices[(cnt + 1) : (cnt + tal)] = toadd
    cnt += tal
  end
  return indices
end

function attboot_inner!(
  bootests,
  clusters,
  outcomemat, dwitmat,
  uid, utrtid, trtidunique,
  Fset, iternum)
  #@inbounds Threads.@threads for k = eachindex(1:iternum)
  for k = eachindex(1:iternum)
    # clusters = dat_id_un;
    units = sample(clusters, length(clusters), replace = true);

    make_check_sample!(units, clusters, trtidunique, 0); # 521.1 μs

    unitreps = countmap(units, alg = :auto); # 74 μs
    unitunique = unique(units); # 65 μs
    
    # preallocate index list (just make 2x as easy solution)
    indices = zeros(Int64, length(uid) * 2); # bad for possible estimate
    addindices!(indices, uid, unitunique, unitreps, 0); # 5.9 ms -- !!
    indices = indices[indices .> 0]; # 9.3 μs (up to 11.2 with easy sol)

    trtnums = gettrtnums(indices, uid, utrtid, outcomemat);

    bootests[:, k] = newattcalc(
      @view(dwitmat[:, indices]),
      @view(outcomemat[:,indices]),
      Fset,
      trtnums) # this needs a solution
  end
  return bootests
end