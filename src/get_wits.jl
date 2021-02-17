#=
functions to calculate the sub observation level weights
that is, we calculate a weight for each appearance of an observation
as a non-zero weight

the last function makes those weights observation level
=#

"""
    checkmatchnums(mm)

get the number of matches for each treated observation
"""
function checkmatchnums(mm)
  uid = unique(@view(mm[[:tunit, :ttime]]).tunit);
  ut = unique(@view(mm[[:tunit, :ttime]]).ttime);

  nms = similar(uid);

  for i = eachindex(uid)
    nms[i] = length(
      findall((mm.trt_t .== @view(ut[i])) .& (mm.tunit .== @view(uid[i])))
      );
  end
  return uid, ut, nms
end

function getoutcomes!(
  outcol, Fset, tpoint, uti, unit, dat, t, id, outcome)

  K = vcat(uti + tpoint, (uti .+ (Fset)));

  ind = findall(
    (dat[id] .== unit) .& (in.(dat[t], Ref(K))));

  if dat[ind, t][1] != uti + tpoint
    nothing
  elseif length(ind) < length(Fset) + 1;
    q = dat[ind, t] .- uti;
    q[1] += 2
    q[2:end] .-= (minimum(Fset) - 2);
    outcol[q] = dat[outcome][ind];
  else
    outcol[:] = dat[ind, outcome]
  end

  return outcol
end

function makedwits!(
  wcol, Fset, tpoint, uti, munit, trtunit, mnumi)

  fmin = minimum(Fset);
  K = vcat(uti + tpoint, (uti .+ (Fset)));
  for k in K
    
    if (k == K[1]) & (munit == trtunit)
      wcol[1] = -1
    elseif (k == K[1]) & (munit != trtunit)
      wcol[1] = (1 / mnumi)
    elseif (k > K[1]) & (munit == trtunit)
      # matchnum
      wcol[k - uti - fmin + 1 + 1] = 1
    elseif (k > K[1]) & (munit != trtunit)
      # matchnum
      wcol[k - uti - fmin + 1 + 1] = -(1 / mnumi)
    end
  end
  return wcol
end

function getmatchnums!(mnum, uid, utrtid, ut)
  for i = eachindex(uid)
    mnum[i] = length(
      findall(
        (uid .!= utrtid[i]) .& (utrtid .== utrtid[i]) .& (ut .== ut[i])
        )
    )
  end
  return mnum
end

function findneglocs(mat)
  negones = zeros(size(mat)[2]);
  for i = 1:(size(mat)[2])
    negones[i] = 1 * any(mat[:,i] .== -999.0)
  end
  neglocs = findall(negones .== 1);
  return neglocs
end

# maybe better to input the df and extract the cols with @view()
# make this the outermost function

# this is faster than vectorized / concat operation
function makeK!(K, Fset, uti, tpoint)
  for k = eachindex(Fset)
    K[k + 1] = uti + (Fset[k])
  end
  K[1] = uti + tpoint
  return K
end

function populate_outcome_mat!(
  outcomemat, uid, ut, did, dt, doc, tpoint, Fset, tlen, fmin)

  @inbounds Threads.@threads for i = eachindex(uid)
  # for i = eachindex(uid)
    unit = @view(uid[i])[1]
    uti = @view(ut[i])[1]

    # getoutcomes!(
    #   @view(outcomemat[:, i]), Fset, tpoint, uti, unit, dat, t, id, outcome)

    ocmatind, dqq = makeoutcomematrow(
      tlen, Fset, uti, tpoint, did, dt, doc, unit, fmin)

    outcomemat[ocmatind, i] = dqq;

  end
  return outcomemat
end

function makeoutcomematrow(tlen, Fset, uti, tpoint, did, dt, doc, unit, fmin)
  K = zeros(Int64, tlen);
  makeK!(K, Fset, uti, tpoint);

  qq = getvarindset(did, dt, unit, K);

  ocmatind = @view(dt[qq]) .- (uti + fmin - 2)
  ocmatind[1] = 1 # this point must exist
  return ocmatind, doc[qq]
end

function make_dwit_mat!(
  dwitmat, outcomemat, uid, utrtid, ut, mnum, Fset, tpoint)
  
  @inbounds Threads.@threads for i = eachindex(uid)
  # for i = eachindex(uid)
    munit = @view(uid[i])[1]
    trtunit = @view(utrtid[i])[1]
    uti = @view(ut[i])[1]
    mnumi = @view(mnum[i])[1]

    # enfore same missingness from outcomes onto dwits
    # after getting the weights -- easiest solution
    makedwits!(
      @view(dwitmat[:, i]),
      Fset, tpoint, uti, munit, trtunit, mnumi)
    dwitmat[findall(ismissing.(outcomemat))] .= missing;
  end
  return dwitmat
end

function get_dwits_outcomes(
  utrtid, uid, ut,
  Fset,
  did, dt, doc,
  tpoint)

  tlen = length(Fset) + 1 # +1 for tpoint
  lrow = length(uid);

  fmin = minimum(Fset);

  outcomemat = Matrix{Union{Float64, Missing}}(
    missing, (tlen, lrow));
  dwitmat = similar(outcomemat);

  populate_outcome_mat!(
    outcomemat, uid, ut, did, dt, doc, tpoint, Fset, tlen, fmin);

  mnum = zeros(Int64, length(uid));
  getmatchnums!(mnum, uid, utrtid, ut);

  make_dwit_mat!(dwitmat, outcomemat, uid, utrtid, ut, mnum, Fset, tpoint);

  return outcomemat, dwitmat
end