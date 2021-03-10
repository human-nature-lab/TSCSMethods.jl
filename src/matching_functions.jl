using StatsBase, Statistics

import LinearAlgebra
import IterTools.product

#=
new matching functions
=#

# 402 sec with maha dist calculation
# 408 sec with caliper (on two variables)
# TODO
# - clean up the variables
# - make mlen an arg or calculate from inputs, etc.
# possible options: keep/drop impossible matches
# methods: caliper or not, distance or not (maybe options instead?)
function matching(
  dat, covariates::Vector{Symbol},
  id::Symbol, t::Symbol,
  fmin, fmax,
  treatment::Symbol,
  calvars::Vector{Symbol},
  variancesonly::Bool)

  sort!(dat, [id, t]);

  treatment_points = findall(dat[!, treatment] .== 1)

  cvkey = [findfirst(covariates .== e) for e in calvars];

  allvars = vcat([t, id, treatment], covariates);

  did = dat[!, id];
  dt = dat[!, t];

  cmat = convert(Matrix{Float64}, dat[!, covariates]);

  # treatment point reference
  rid = @view(did[treatment_points]);
  rt = @view(dt[treatment_points]);
  ky = 1:length(rid);
  
  unid = unique(did);

  mlen = 20;

  pp = collect(product(ky, unid)); # collect since required by @threads
  ppl = length(pp);

  idI = zeros(Int64, ppl);
  tI = zeros(Int64, ppl);
  idJ = zeros(Int64, ppl);
  possible = zeros(Int64, ppl);
  mdist = Inf .* zeros(Float64, ppl);
  caldists = Inf .* zeros(Float64, length(cvkey), ppl);

  ΣΣ, ΣΣtimes = get_Σs(cmat, dt, unique(rt), mlen, tpoint);

  matching_inner!(
    pp, rid, rt,
    idI, tI, idJ, possible, mdist, caldists,
    did, dt, cmat,
    cvkey,
    fmin, fmax,
    ΣΣ, ΣΣtimes, mlen,
    variancesonly)

  matches = DataFrame(
    ttime = tI, tunit = idI, munit = idJ,
    possible = possible,
    mdist = mdist);

  caldf = DataFrame(caldists')
  calvarsnmes = [Symbol(String(calvar) * "_mdist") for calvar in calvars]
  rename!(caldf, calvarsnmes)

  matches = hcat(matches, caldf)

  keepind = findall((isnan.(matches.mdist) .- 1) .* - 1 .== 1);
  matches = matches[keepind, :]

  # make sure there is a match for every treated unit
  matches1 = groupby(matches, [:ttime, :tunit]);
  combine(
    matches1,
    nrow
  )

  #= this is an inefficient solution to an apparent problem. make better.
  it should happen, most efficiently, in the inner function =#
  tobst = unique(matches[!, [:tunit, :ttime]])[!, :ttime];
  tobsu = unique(matches[!, [:tunit, :ttime]])[!, :tunit];

  rmlocs = Vector{Int64}(undef, 0);
  
  for i = eachindex(tobst)
    locs = findall(
      (matches[:ttime] .== tobst[i]) .& (matches[:tunit] .== tobsu[i])
    )

    if length(locs) == 1 # if unit has no matches
      append!(rmlocs, locs)
    end
  end

  delete!(matches, rmlocs)

  sort!(matches, [:ttime, :tunit, :possible, :mdist])
  
  return matches
end

function matching_inner!(
  pp, rid, rt,
  idI, tI, idJ, possible, mdist, caldists,
  did, dt, cmat, # dat,
  cvkey,
  fmin, fmax,
  ΣΣ, ΣΣtimes, mlen,
  variancesonly)

  @inbounds Threads.@threads for i = eachindex(pp)
  # for (i, p) in enumerate(pp)

    idI[i] = rid[pp[i][1]];
    tI[i] = rt[pp[i][1]];
    idJ[i] = pp[i][2];

    idi = @views(idI[i]);
    ti = @views(tI[i]);
    idj = @views(idJ[i]);

    tmn = max(ti - mlen - fmax, 0);
    tmx = ti + fmax - fmin;

    if (idi == idj)
      mdist[i] = 0.0
      caldists[:, i] .= 0.0

    # the match (idj) must exist at or before tmn
    # save time by relying on [id, t] order (28 seconds w/ 2mill len loop)
    # elseif (minimum(dt[did .== idj]) <= tmn)
    elseif dt[findfirst(did .== idj)] <= tmn;
    
      idpl = findall(rid .== idj)
      tref = @view(rt[idpl])

      v = 1 * !any((tref .>= tmn) .& (tref .<= tmx))

      possible[i] = v

      if (v == 1)

        ix = getvarind(
          did, dt, idi, (ti - 20) : (ti - 1)); # t-unit
        jx = getvarind(
          did, dt, idj, (ti - 20) : (ti - 1)); # m-unit
          
        if (length(ix) == length(jx) == mlen)

          Σs = indexΣΣ(ΣΣ, ΣΣtimes, ti) # slow for index

          # length check here? somewhere?
          y = @view(cmat[ix, :]) # just sel and arrange columns up front
          x = @view(cmat[jx, :])

          dissum = [0.0] # Float64 is immutable
          calsums = zeros(Float64, length(calvars))
          
          #try 
          matching_distancecalc_barrier!(
            dissum, calsums, y, x, Σs, cvkey, covariates, mlen,
            variancesonly)
          #catch
          #  @warn "error at i = " * string(i)
          #end

          # avg dist between treated observation and possible match
          # over match period
          mdist[i] = @views(dissum[1]) / mlen
          caldists[:, i] = calsums / mlen
        end
      end
    end
  end
  return idI, tI, idJ, possible, mdist, caldists
end

function matching_distancecalc_barrier!(
  dissum, calsums, y, x, Σs, cvkey, covariates, mlen,
  variancesonly)
  for l = eachindex(1:mlen)
    
    Σ_t = @view(Σs[:, :, l]); # slow?
    
    yy = @views(y[l, :]);
    xx = @views(x[l, :]);

    diffs = yy - xx; # trt - cntrl
    # maha dist at t-l:
    if variancesonly == true
      dissum[1] += sqrt(
        transpose(diffs) * LinearAlgebra.pinv(
          LinearAlgebra.Diagonal(Σ_t)
        ) * diffs
        # changed to pinv()
      );
    else
      dissum[1] += sqrt(
        transpose(diffs) * LinearAlgebra.pinv(Σ_t) * diffs);
        # changed to pinv()
    end

    # caliper var calc here
    # calsums[m,i] rows on inner loop
    # NEW METHOD, NOT CONDITION
    # diag as 1d index for a given covariate
    dval = (cvkey .- 1) .* length(covariates) + cvkey
    σ_t = @view(Σ_t[dval]);

    matching_distancecalc_caliper_barrier!(calsums, y, x, cvkey, σ_t, l)
  end
  return dissum, calsums
end

function matching_distancecalc_caliper_barrier!(calsums, y, x, cvkey, σ_t, l)
  for (h, m) in enumerate(cvkey) # key from names to cols in dat
    yyc = @views(y[l, m]);
    xxc = @views(x[l, m]);
    caldiff = yyc - xxc
    calsums[h] += sqrt(caldiff * 1 / @views(σ_t[h]) * caldiff)
  end
  return calsums
end

function refine(matches::DataFrame, setsize::Int64)
    sort!(matches, [:ttime, :tunit, :possible, :mdist])

    tobs = unique(matches[!, [:ttime, :tunit]]);

    indls = Vector{Int64}(undef, 0);
    ss = setsize + 1

    for i = eachindex(1:nrow(tobs))
      ind = findall(
        (matches[!, :ttime] .== tobs[i, :ttime]) .&
        (matches[!, :tunit] .== tobs[i, :tunit])
        )
      if length(ind) >= (ss)
        append!(indls, ind[1 : ss])
      end
    end
    return matches[indls, :]
end

"""
    caliper(cdict::Dict, matches::DataFrame)

Enforce a caliper on desired set of covariates, or overall Mahalanabois
distance.

Examples
≡≡≡≡≡≡≡≡≡≡

cdict = Dict(:cum_death_rte => 0.5, :first_case => 0.5);
matchescal = caliper(cdict, matches);
"""
function caliper(cdict::Dict, matches5::DataFrame)
  ntobsp = nrow(unique(matches5[!, [:ttime, :tunit]]));
  # need to add mdist

  keep = ones(Int64, nrow(matches5));
  
  for (k, v) in cdict
    if k != :mdist
      k = Symbol(String(k) * "_mdist")
      keep = keep .& (matches5[!, k] .< v)
    end
  end

  keeploc = findall(keep .== 1);

  matches5 = matches5[keeploc, :];

  matches5 = @linq matches5 |>
    groupby([:ttime, :tunit]) |>
    combine(remove = sum(:mdist)) |>
    leftjoin(matches5, on = [:tunit, :ttime])

  lost = @where(matches5, :remove .== 0);
  lost = unique(lost[!, [:ttime, :tunit]]);

  matches5 = @where(matches5, :remove .> 0);

  matches5 = select(matches5, Not(:remove));
    
  ntobs = nrow(unique(matches5[!, [:ttime, :tunit]]))
  println(
    string(ntobs) *
    " out of " *
    string(ntobsp) *
    " treated observations remain."
    )

  return matches5, lost
end

#=
needed old code, needed, likely to be modified
=#

# calculate the standard deviations
# need all c_list covariates in a dataframe, with
# a col to standardize the running vals (with trt as day 1)
# then, work with groupby for each time point

function get_trt_sd(dat, c_list, mm, tmin, id, t)
  c_types = get_c_types(dat, covariates)

  # trt_obs = dat[treatment_points, [id, t]];
  trt_obs = unique(mm[!, [:trt_fips, :trt_t]]);

  # need to get treatment points in dat from
  # trt_obs
  trt_obs.trt_pts = 0
  for i = 1:nrow(trt_obs)
    trt_obs.trt_pts[i] = findall(
      (dat.fips .== trt_obs.trt_fips[i]) .&
      (dat.running .== trt_obs.trt_t[i])
    )[1]
  end
  rename!(trt_obs, [id, t, :trt_pts])

  trt_covs_nrow = nrow(trt_obs) * -tmin; # only works with tmin before trt
  trt_covs = DataFrame(
    vcat(c_types, [Int64, Int64, Int64]),
    vcat(c_list, [t, :mtch_t, id]),
    trt_covs_nrow);
  
  # trt_covs[c_list[1]] = Int64[]
  for i in 1:nrow(trt_obs)
    fps = trt_obs[id][i]
    trt_run = trt_obs[t][i]
    run_min = trt_run + tmin

    i_obs = get_trt_covariates(
      trt_run, fps, dat, c_list, tmin, id, t, delay = 0);
    i_obs = @transform(
      i_obs,
      mtch_t = :running .- (minimum(:running) - 1)
      )
    i_obs[id] = fps;
    # earliest matching period is 1
    # highest is 20, just before treatment onset

    # put into dataframe -- we a priori know its length this time
    # only work with tmin before trt (negative in value)
    idex_s = (-tmin) * (i - 1) + 1
    idex_e = (-tmin) * i
    trt_covs[idex_s:idex_e, :] = i_obs;
  end
  return trt_covs
end

function get_Σs(cmat, dt, TT, mlen, tpoint)

  sort!(TT)

  xy_dim = size(cmat)[2];
  ΣΣ = Array{Union{Float64, Missing}}(
    missing, xy_dim, xy_dim, mlen, length(TT)
  );
  ΣΣtimes = Vector{Int64}(undef, length(TT));
  ΣΣ, ΣΣtimes = get_Σs_inner(
      cmat, dt,
      ΣΣ, ΣΣtimes,
      TT, tmin, tpoint, xy_dim)
  return ΣΣ, ΣΣtimes
end

function get_Σs_inner(
  cmat, dt,
  ΣΣ, ΣΣtimes,
  TT, tmin, tpoint, xy_dim)
  for i = eachindex(TT)
    tt = TT[i]
    match_times = [(tt + tmin) : 1 : (tt + tpoint);]
    # storage? 3x3x20 x length(uniq_treat_times)
    ΣΣ[:, :, :, i] = get_sample_covars(
      match_times, dt, cmat, xy_dim)
    ΣΣtimes[i] = tt
  end
  return ΣΣ, ΣΣtimes
end

function get_sample_covars(
  match_times, dt, cmat, xy_dim)
  
  Σs = Array{Union{Float64, Missing}}(
    undef, xy_dim, xy_dim, length(match_times)
  );
  t_adjust = minimum(match_times) - 1;

  get_sample_covar_t!(Σs, match_times, dt, cmat, t_adjust, xy_dim);
  
  return Σs
end

function get_sample_covar_t!(Σs, match_times, dt, cmat, t_adjust, xy_dim)
  
  for mt in match_times
    mtind = findall(dt .== mt);
    #=
    essentially, the number of units that exist at some match point
    there should maybe be a lower bound above 1 here for calculating the sample
    covariance mat at that time
    =#
    if length(mtind) > 0
      MCov_t = StatsBase.cov(@view(cmat[mtind, :]));
    else
      MCov_t = Matrix{Union{Float64, Missing}}(missing, xy_dim, xy_dim);
    end

    Σs[:, :, mt - t_adjust] = MCov_t
  end

  return Σs
end

"""
      indexΣΣ(ΣΣ, ΣΣtimes, trt_t)

helper function to grab the correct set of sample covars:
there is a unique one for each treatment time in the data
"""
function indexΣΣ(ΣΣ, ΣΣtimes, trt_t)
  idex = findall(ΣΣtimes .== trt_t)
  Σs = ΣΣ[:, :, :, idex];
  return reshape(Σs, size(Σs)[1], size(Σs)[2], size(Σs)[3])
end

# functions required by functions to be modified

"""
    get_trt_covariates(trt_t, trt_fips, dat, c_list, tmin)

this method grabs the covariates in the matching period
"""
function get_trt_covariates(trt_t::Int,
                            trt_fips::Int,
                            dat::DataFrame,
                            c_list::Array{Symbol,1},
                            tmin,
                            id::Symbol, t::Symbol;
                            delay::Int = 0)
  t_min = trt_t + tmin

  y = findall(
    (dat[!, t] .>= t_min) .&
      (dat[!, t] .< (trt_t + delay)) .&
      (dat[!, id] .== trt_fips))

  return dat[y, vcat(c_list, t)]
end

"""
    get_trt_covariates(trt_t, trt_fips, dat, c_list, tmin, tmax)

this should be made more general, to get covariates
either for the match period, or some post-primary period

tmin, days before/after treatment, tmax = days before/after treatment
constraint: tmin < tmax
"""
function get_trt_covariates(trt_t::Int,
                            trt_fips::Int,
                            dat::DataFrame,
                            c_list::Array{Symbol,1},
                            tmin::Int,
                            tmax::Int,
                            id::Symbol, t::Symbol)

  if (tmin >= tmax)
    error("tmax must be strictly larger than tmin")
  end

  t_min = trt_t + tmin
  t_max = trt_t + tmax

  y = findall(
    (dat[!, t] .>= t_min) .&
    (dat[!, t] .<= t_max) .&
    (dat[!, id] .== trt_fips))

  return dat[y, vcat(c_list, :running)]
end

# # this should get the treated points directly from the data
# function sdtreatedold(dat, covariates, tmin, id, t, treatment)
  
#   c_types = get_c_types(dat, covariates)
  
#   treatment_points = get_all_treatment_points(dat, treatment);

#   trt_obs = dat[treatment_points, [id, t]];

#   # need to get treatment points in dat from
#   # trt_obs
#   trt_obs.trt_pts = zeros(Int64, nrow(trt_obs))
#   for i = eachindex(1:nrow(trt_obs))
#     a = dat.fips .== trt_obs.fips[i]
#     b = dat[!, t] .== trt_obs[!, t][i]
#     trt_obs.trt_pts[i] = findall(a .& b)[1]
#   end

#   trt_covs_nrow = nrow(trt_obs) * -tmin; # only works with tmin before trt
#   trt_covs = DataFrame(
#     vcat(c_types, [Int64, Int64, Int64]),
#     vcat(covariates, [t, :mtch_t, id]),
#     trt_covs_nrow);
  
#   # trt_covs[c_list[1]] = Int64[]
#   for i = eachindex(1:nrow(trt_obs))
#     fps = trt_obs[!, id][i]
#     trt_run = trt_obs[!, t][i]
#     run_min = trt_run + tmin

#     i_obs = get_trt_covariates(
#       trt_run, fps, dat, covariates, tmin, id, t, delay = 0);
#     i_obs = @transform(
#       i_obs,
#       mtch_t = :running .- (minimum(:running) - 1)
#       )
#     i_obs[:, id] = fps * ones(Int64, nrow(i_obs));
#     # earliest matching period is 1
#     # highest is 20, just before treatment onset

#     # put into dataframe -- we a priori know its length this time
#     # only work with tmin before trt (negative in value)
#     idex_s = (-tmin) * (i - 1) + 1
#     idex_e = (-tmin) * i
#     trt_covs[idex_s:idex_e, :] = i_obs;
#   end

#   trt_covs = groupby(trt_covs, :mtch_t);
#   trt_std = combine(
#     trt_covs,
#     covariates .=> std,
#     );
#   new_names_std = replace.(names(trt_std), "_std" => "")
#   rename!(trt_std, new_names_std)
#   return trt_std
# end