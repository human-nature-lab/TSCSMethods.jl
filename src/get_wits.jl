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
  checkmatchnums_inner!(nms, uid, ut, mm)

  return uid, ut, nms
end

function checkmatchnums_inner!(nms, uid, ut, mm)
  for i = eachindex(uid)
    nms[i] = length(
      findall((mm.trt_t .== @view(ut[i])) .& (mm.tunit .== @view(uid[i])))
      );
  end
  return nms
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

"""
this is the new, faster function to make the outcome and dwit mats
we cannot simply add the wits to the dataframe, since we are keeping them
disaggregated
"""
function makemats!(om, wm, uid, utid, ut, mcnts, nd, tpoint, fmin, fmax)
  for (j, φ) in enumerate(vcat(tpoint , collect(fmin:fmax))) # over f
    @inbounds Threads.@threads for i = eachindex(uid) # over matches
      mu = @views(uid[i])
      tu = @views(utid[i])
      tt = @views(ut[i])
      mc = @views(mcnts[(tt, tu)])

      om[i, j] = get(nd, (mu, tt + φ), missing)

      if ismissing(om[i, j])
        wm[i, j] = missing
      else
        wm[i, j] = witcal(mu, tu, j, mc)
      end
    end
  end
  return om, wm
end

"""
    missingmats(om, wm)

output:
- adjusted om
- adjusted wm
- bl: black list of treated observations with treatment period (f)
"""
function missingmats(
  om::Array{Union{Missing, Float64}, 2},
  wm::Array{Union{Missing, Float64}, 2},
  uid, utid, ut, utrtid
)

  # check bl
  tulx = findall((uid .== utrtid));
  sum(ismissing.(om[tulx, :])) # this difference could be the problem

  hcat(
    sum(ismissing.(om[tulx, :]), dims = 2),
    utrtid[tulx]
  )
  
  # MM = DataFrame(ismissing.(outcomemat[tulx, :]));
  # MM.tunit = utrtid[tulx];
  # MM.tt = ut[tulx];

  # MM = stack(MM, Not([:tunit, :tt]))
  # MM = MM[MM.value .!= 0, :]
  # MM.variable = string.(MM.variable);
  # [MM.variable[i] = split(MM.variable[i], "x")[2] for i in eachindex(MM.variable)]
  # MM.variable = parse.(Int64, MM.variable)
  # bl.count = 1
  # MC = leftjoin(MM, bl, on = [:tunit => :id, :tt => :t, :variable => :f])
  # mc = @where(MC, ismissing.(:count))

  # upper bound: number of missing entries
  # possible missing num treated x f
  poss = sum(uid .== utrtid) * size(om)[2]

  # idbl = Vector{Int64}(undef, poss)
  # tbl = Vector{Int64}(undef, poss)
  # fbl = Vector{Int64}(undef, poss)

  om = deepcopy(om);
  wm = deepcopy(wm);

  idbl = zeros(Int64, poss, Threads.nthreads());
  tbl = zeros(Int64, poss, Threads.nthreads());
  fbl = zeros(Int64, poss, Threads.nthreads());

  cnts = zeros(Int64, Threads.nthreads());

  om, wm, idbl, tbl, fbl = innermissing(
    idbl, tbl, fbl,
    uid, utid, ut,
    om, wm, cnts
  )
  
  idbl = reshape(idbl, poss * Threads.nthreads())
  tbl = reshape(tbl, poss * Threads.nthreads())
  fbl = reshape(fbl, poss * Threads.nthreads())

  # these treated obs are missing at 
  bl = unique(DataFrame(id = idbl, t = tbl, f = fbl))
  bl = bl[bl.id .!= 0, :]
  return om, wm, bl
end

function innermissing(
  idbl, tbl, fbl,
  uid, utid, ut,
  om, wm, cnt
)
  @inbounds Threads.@threads for j in 2:size(om)[2]
    q = Threads.threadid()
    for i = eachindex(uid)
      # cannot use Threads.@threads with push!
      # cannot use @inbounds -- need safe

      ui = uid[i]
      uti = utid[i]
      tt = ut[i]

      # if outcome is missing, make all other uti x ut values missing
      if (ui == uti) & ismissing(@views(om[i, j]))
        c1 = ut .== tt;
        c2 = utid .== uti;
        wm[c1 .& c2, j] .= 0.0 # give all matches to missing treated outcome weight zero
        om[c1 .& c2, j] .= 0.0

        # add treated x t x f to blacklist
        # append!(idbl, uti)
        # append!(tbl, tt)
        # append!(fbl, j)
        cnt[q] += 1
        idbl[cnt[q], q] = uti
        tbl[cnt[q], q] = tt
        fbl[cnt[q], q] = j
      elseif (ui != uti) & ismissing(om[i, j])
        c1 = ut .== tt
        c2 = utid .== uti # c2 and c3 problem?
        c3 = utid .!= uid
        chk = ismissing.(@views(om[c1 .& c2 .& c3, j]))
        if sum(chk) == length(chk) # if all matches are missing
          wm[c1 .& c2, j] .= 0.0 # give all matches to missing treated outcome zero weight & outcome, also for treated unit itself
          om[c1 .& c2, j] .= 0.0
          # could the above just be [i, j]?
          
          # add treated x t x f to blacklist
          # append!(idbl, uti)
          # append!(tbl, tt)
          # append!(fbl, j)
          cnt[q] += 1
          idbl[cnt[q], q] = uti
          tbl[cnt[q], q] = tt
          fbl[cnt[q], q] = j
        elseif sum(chk) < length(chk)
          # change weights of remaining
          c4 = (!).(ismissing.(@views(om[:, j])))
          chk2 = @views(om[c1 .& c2 .& c3 .& c4, j])
          wm[c1 .& c2 .& c3 .& c4, j] .= -1.0 / length(chk2) # matches so negative
          om[i, j] = 0.0
          wm[i, j] = 0.0
        end
      end
    end
  end
  return om, wm, idbl, tbl, fbl
end

# slow because operating over rows...
# cannot summarize to above unit
function summarizemats(uid, ut, om, wm, bl, mk)
  uid = convert(Vector{Int64}, uid);
  # ut
  # om
  # wm

  obs = unique(hcat(uid, ut), dims = 1);

  # tobs = unique(matches5[[:ttime, :tunit]])

  omsm = Matrix{Float64}(undef, size(obs)[1], size(om)[2]); # depends on tranpose situation
  wmsm = similar(omsm);

  trtsm = zeros(Int64, size(obs)[1]);

  if size(bl)[1] > 0
    blmat = zeros(Int64, size(obs)[1], size(omsm)[2] - 1);
  else blmat = zeros(Int64, 0, 0);
  end

  summarizemats_inner!(obs, mk, trtsm, bl, blmat, uid, ut, omsm, wmsm, om, wm)

  return omsm, wmsm, obs[:, 1], obs[:, 2], trtsm, blmat
end

function summarizemats_inner!(
  obs, mk, trtsm, bl, blmat, uid, ut, omsm, wmsm, om, wm
)

  for i in 1:size(obs)[1]
    lid = @views(obs[i , 1]) # id
    lt = @views(obs[i , 2]) # time

    summarizemats_inner_inner!(trtsm, i, mk, lt, lid)

    if size(bl)[1] > 0
      blf = @views(bl[(bl[!, :id] .== lid) .& (bl[!, :t] .== lt), :f]) .- 1;
      if length(blf) > 0
        blmat[i, blf] .= 1
      end
    end

    found = findall((uid .== lid) .& (ut .== lt))

    omsm[i, :] = om[@view(found[1]), :]
    wmsm[i, :] = sum(@view(wm[found, :]), dims = 1)
  end
  return omsm, wmsm, blmat, trtsm
end

function summarizemats_inner_inner!(trtsm, i, mk, lt, lid)
  # for j in 1:length(tobs.ttime)
  for (kt, kid) in mk
    # if (lt == tobs.ttime[j]) & (lid == tobs.tunit[j])
    if (lt == kt) & (lid == kid)
      trtsm[i] += 1 # can only ever add a total of one for each i anyway
      continue # may save some time
    end
  end
  return trtsm
end

"""
witcal(mu, tu, j, mc)

get the disaggregated weight value for an observation. this function is only applied to observations that are non-zero. (although else statement implies that it could.)

Arguments
≡≡≡≡≡≡≡≡≡≡≡

- mu: match unit
- tu: treated unit for mu
- j: vector position for relative timepoint, j == 1 => tpoint, j > 1 => f in F
- mc: number of matches to the tu (of which mu is one)
"""
function witcal(mu, tu, j, mc)
  cif = mu == tu
  c1st = (j == 1)

  if cif & c1st
    ς = -1.0
  elseif cif & !c1st
    ς = 1.0
  elseif !cif & c1st
    ς = 1.0 / mc
  elseif !cif & !c1st
    ς = -1.0 / mc
  else
    ς = 0.0
  end
  return ς
end