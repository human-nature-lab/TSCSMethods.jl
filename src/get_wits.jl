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

#= this would need to be written for summarized dwits,
  adding the value to each entry of dw, each time 
=#
# function assign_dwits_to_dat!(dw, utrtid, uid, ut, mcnt)
#   # dw = zeros(Float64, length(did));

#   # utrtid # treated
#   # uid # matches and treated
#   # ut # times

#   for i = eachindex(1:length(ut))
#     tu = utrtid[i]
#     mu = uid[i]
#     tt = ut[i]
#     mc = mcnt[(tt, tu)]

#     ctp = dt .== tt + tpoint;
#     cfpl = dt .>= tt + fmin;
#     cfpu = dt .<= tt + fmax;
#     c2 = did .== mu;
    
#     cif = mu == tu;
#     if cif
#       dw[ctp .& c2][] = -1.0
#       dw[cfpl .& cfpu .& c2] .= 1.0
#     else
#       dw[ctp .& c2][] = 1.0/mc
#       dw[cfpl .& cfpu .& c2] .= -1.0/mc # try countmap
#     end

#   end
#   return dw
# end

"""
this is the new, faster function to make the outcome and dwit mats
we cannot simply add the wits to the dataframe, since we are keeping them
disaggregated
"""
function makemats!(om, wm, uid, utid, ut, mcnts, nd, tpoint)
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
  uid, utid
  )

  idbl = Vector{Int64}()
  tbl = Vector{Int64}()
  fbl = Vector{Int64}()
  for j in 2:size(om)[2]
    @inbounds for i = eachindex(uid)

      ui = uid[i]
      uti = utid[i]
      tt = ut[i]

      # if outcome is missing, make all other uti x ut values missing
      if (ui == uti) & ismissing(om[i, j])
        c1 = ut .== tt;
        c2 = utid .== uti;
        wm[c1 .& c2, j] .= 0 # give all matches to missing treated outcome weight zero
        om[c1 .& c2, j] .= 0

        # add treated x t x f to blacklist
        push!(idbl, uti)
        push!(tbl, tt)
        push!(fbl, j)
      elseif (ui != uti) & ismissing(om[i, j])
        c1 = ut .== tt
        c2 = utid .== uti
        c3 = utid .!= uid
        chk = ismissing.(@views(om[c1 .& c2 .& c3, j]))
        if sum(chk) == length(chk) # if all matches are missing
          wm[c1 .& c2, j] .= 0 # give all matches to missing treated outcome zero weight & outcome, also for treated unit itself
          om[c1 .& c2, j] .= 0
          
          # add treated x t x f to blacklist
          push!(idbl, uti)
          push!(tbl, tt)
          push!(fbl, j)
        elseif sum(chk) < length(chk)
          # change weights of remaining
          c4 = (!).(ismissing.(om[:, j]))
          chk2 = @views(om[c1 .& c2 .& c3 .& c4, j])
          wm[c1 .& c2 .& c3 .& c4, j] .= -1.0 / length(chk2) # matches so negative
        end
      end
    end
  end   
  # these treated obs are missing at 
  bl = unique(DataFrame(id = idbl, t = tbl, f = fbl))
  return om, wm, bl
end

# slow because operating over rows...
# cannot summarize to above unit
function summarizemats(uid, ut, om, wm, bl, mk)
  uid2 = convert(Vector{Int64}, uid);
  # ut
  # om
  # wm

  obs = unique(hcat(uid2, ut), dims = 1);

  # tobs = unique(matches5[[:ttime, :tunit]])

  omsm = Matrix{Float64}(undef, size(obs)[1], size(om)[2]); # depends on tranpose situation
  wmsm = similar(omsm);

  trtsm = zeros(Int64, size(obs)[1]);

  if size(bl)[1] > 0
    blmat = zeros(Int64, size(obs)[1], size(omsm)[2] - 1);
  else blmat = zeros(Int64, 0, 0);
  end
  for i in 1:size(obs)[1]
    lid = @views(obs[i , 1]) # id
    lt = @views(obs[i , 2]) # time

    # for j in 1:length(tobs.ttime)
    for (kt, kid) in mk
      # if (lt == tobs.ttime[j]) & (lid == tobs.tunit[j])
      if (lt == kt) & (lid == kid)
        trtsm[i] += 1 # can only ever add a total of one for each i anyway
        continue # may save some time
      end
    end

    if size(bl)[1] > 0
      blf = @views(bl[(bl[!, :id] .== lid) .& (bl[!, :t] .== lt), :f]) .- 1;
      if length(blf) > 0
        blmat[i, blf] .= 1
      end
    end

    found = findall((uid2 .== lid) .& (ut .== lt))

    omsm[i, :] = om[@views(found[1]), :]
    wmsm[i, :] = sum(@views(wm[found, :]), dims = 1)
  end

  return omsm, wmsm, obs[:, 1], obs[:, 2], trtsm, blmat
end

        


# function makemats_old!(om, wm, uid, utid, ut, mcnts, nd, tpoint)
#   # om = similar(outcomemat);

#   @inbounds Threads.@threads for j = eachindex(uid)
#   # for j = eachindex(uid)
#     mu = @views(uid[j])
#     tu = @views(utid[j])
#     tt = @views(ut[j])
#     mc = @views(mcnts[(tt, tu)])

#     for (i, φ) in enumerate(vcat(tpoint , collect(fmin:fmax)))
#       # outcomemat[i, j]
#       om[i, j] = get(nd, (mu, tt + φ), missing)

#       if ismissing(om[i, j])
#         wm[i, j] = missing
#       else
#         wm[i, j] = witcal(mu, tu, i, mc)
#       end
#     end
#   end
#   return om, wm
# end

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