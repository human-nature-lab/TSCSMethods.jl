#=
new matching functions
=#

function get_Σs!(Σdict, cmat, covnum, dt, uniqt, mmin, mmax)

  for tt in unique(uniqt)
    for mt in mmin:mmax
      inx = findall(dt .== tt + mt);
      if length(inx) > 0
        Σdict[(tt, mt)] = StatsBase.cov(@views(cmat[inx, :]));
      else
        Σdict[(tt, mt)] = Matrix{Union{Float64, Missing}}(
          missing, covnum, covnum
        );
      end
    end
  end

  return Σdict
end

# 402 sec with maha dist calculation
# 408 sec with caliper (on two variables)
# TODO
# - clean up the variables
# - make mlen an arg or calculate from inputs, etc.
# possible options: keep/drop impossible matches
# methods: caliper or not, distance or not (maybe options instead?)
function matching(
  dat::DataFrame,
  covariates::Vector{Symbol},
  id::Symbol, t::Symbol,
  fmin, fmax,
  mmin, mmax,
  treatment::Symbol,
  calvars::Vector{Symbol},
  variancesonly::Bool
)

  sort!(dat, [id, t]);

  treatment_points = findall(dat[!, treatment] .== 1);

  cvkey = [findfirst(covariates .== e) for e in calvars];

  allvars = vcat([t, id, treatment], covariates);

  did = dat[!, id];
  dt = dat[!, t];
  dtrt = dat[!, treatment];

  # this is kind of wasteful
  # this needs to be sorted a bit better
  cmat = convert(Matrix{Float64}, Matrix(dat[!, covariates]));

  # treatment point reference
  rid = @view(did[treatment_points]);
  rt = @view(dt[treatment_points]);
  #  ky = 1:length(rid);
  
  unid = unique(did);

  # tpoint: reference day (usually t - 1)
  # mmax: first day of match period w.r.t treatment day (NOT tpoint)
  # mmin: last day of match period w.r.t treatment day (NOT tpoint)

  mlen = length(mmin:mmax)

  # pp = collect(product(ky, unid)); # collect since required by @threads
  # ppl = length(pp);

  L = length(unid) * length(treatment_points);

  idtreat = zeros(Int64, L);
  ttimes = zeros(Int64, L);
  idmatch = zeros(Int64, L);
  possible = [false for i in 1:L];
  mdist = Inf .* zeros(Float64, L);
  caldists = Inf .* zeros(Float64, length(cvkey), L);

  Σdict = Dict{Tuple{Int64, Int64}, Array{Union{Float64, Missing}, 2}}();
  get_Σs!(Σdict, cmat, size(cmat)[2], dt, unique(rt), mmin, mmax);

  # matching_inner!(
  #   pp, rid, rt,
  #   idI, tI, idJ, possible, mdist, caldists,
  #   did, dt, cmat,
  #   calvars,
  #   cvkey,
  #   fmin, fmax,
  #   Σdict,
  #   mmin, mmax, mlen,
  #   variancesonly
  # )

  matching_inner_alt!(
    # preall
    idtreat, idmatch, ttimes,
    possible, mdist, caldists,
    # data
    rid, rt,
    did, dt, dtrt, cmat,
    unid,
    Σdict,
    # params
    fmin, fmax,
    mmin, mmax, mlen,
    calvars, cvkey,
    # options
    variancesonly;
    nopriortrt = true
  );

  matches = DataFrame(
    ttime = ttimes,
    tunit = idtreat,
    munit = idmatch,
    possible = possible,
    mdist = mdist
  );

  caldf = DataFrame(permutedims(caldists), :auto)
  calvarsnmes = [Symbol(String(calvar) * "_mdist") for calvar in calvars]
  rename!(caldf, calvarsnmes);

  matches = hcat(matches, caldf);

  # keepind = findall((isnan.(matches.mdist) .- 1) .* - 1 .== 1);
  # matches = matches[keepind, :]

  matches = matches[isnan.(matches.mdist) .== false, :];

  #= this is an inefficient solution to an apparent problem. make better.
  it should happen, most efficiently, in the inner function =#
  tobst = unique(matches[!, [:tunit, :ttime]])[!, :ttime];
  tobsu = unique(matches[!, [:tunit, :ttime]])[!, :tunit];

  rmlocs = Vector{Int64}(undef, 0);
  
  for i = eachindex(tobst)
    locs = findall(
      (matches[!, :ttime] .== tobst[i]) .& (matches[!, :tunit] .== tobsu[i])
    )

    if length(locs) == 1 # if unit has no matches
      append!(rmlocs, locs)
    end
  end

  delete!(matches, rmlocs)

  sort!(matches, [:ttime, :tunit, :possible, :mdist])
  
  return matches
end

function matching_inner_alt!(
  # preall
  idtreat, idmatch, ttimes,
  possible, mdist, caldists,
  # data
  rid, rt,
  did, dt, dtrt, cmat,
  unid,
  Σdict,
  # params
  fmin, fmax,
  mmin, mmax, mlen,
  calvars, cvkey,
  # options
  variancesonly;
  nopriortrt = true
)

  # for i = eachindex(rid)
  @inbounds Threads.@threads for i = eachindex(rid)

    tu = @views(rid[i])
    tt = @views(rt[i])

    tmn = tt + mmin - fmax;

    ftu = dt[findfirst(did .== tu)]; # first time that treated unit exists
    # only proceed if data for treated exists during whole match period
    # not back to tmn, so: far enough back or first point in data
    # "or first point" assumes no treatments prior to data...

    if (ftu > tt + mmin) .& nopriortrt
      # set possible values to zero by default for the combinations
      continue # if treatment, for tu, is too early skip over the tobs and its possible matches
    elseif (ftu > tmn) .& !nopriortrt
      continue # if we don't assume no prior treatment, the data needs to go further back
    else 
      tmx = tt + fmax - fmin;

      # select only data in time range that we need to check for possible matches
      # maybe faster than repeatedly searching the whole dataset?
      ct = ((dt .>= tmn) .& (dt .<= tmx));
      dtsub = @views(dt[ct]);
      didsub = @views(did[ct]);
      dtrtsub = @views(dtrt[ct]);
      cmatsub = @views(cmat[ct, :]);

      tmnadj = max(tmn, minimum(dt));

      # t-unit
      tux = getvarind(didsub, dtsub, tu, (tt + mmin), (tt + mmax));

      # section to alter for a thread
      jsect = ((i - 1) * length(unid) + 1) : (length(unid) * i);

      idtreatr = @views(idtreat[jsect]);
      idmatchr = @views(idmatch[jsect]);
      ttimesr = @views(ttimes[jsect]);
      possibler = @views(possible[jsect]);
      mdistr = @views(mdist[jsect]);
      caldistsr = @views(caldists[:, jsect]);

      idtreat[jsect],
      idmatch[jsect],
      ttimes[jsect],
      possible[jsect],
      mdist[jsect],
      caldists[:, jsect] = 
      matchdistances(
        idtreatr,
        idmatchr,
        ttimesr,
        possibler,
        mdistr,
        caldistsr,
        didsub, dtsub, dtrtsub, cmatsub, # data
        tmn, tmnadj,
        unid,
        tt, tu,
        Σdict,
        calvars, cvkey,
        mmin, mmax, mlen,
        tux,
        nopriortrt,
        variancesonly
      );

    end
  # i loop level
  end
  return idtreat, idmatch, ttimes, possible, mdist, caldists
end
 
"""
will operate on
(i - 1) * length(unid) + 1 : length(unid) * i
"""
function matchdistances(
  idtreatr, idmatchr, ttimesr, # preall
  possibler, mdistr, caldistsr,
  didsub, dtsub, dtrtsub, cmatsub, # data
  tmn, tmnadj,
  unid,
  tt, tu,
  Σdict,
  calvars, cvkey,
  mmin, mmax, mlen,
  tux,
  nopriortrt,
  variancesonly
)

  # reference / checking
  # idtreatr = @views(idtreat[jsect]);
  # idmatchr = @views(idmatch[jsect]);
  # ttimesr = @views(ttimes[jsect]);
  # possibler = @views(possible[jsect]);
  # mdistr = @views(mdist[jsect]);
  # caldistsr = @views(caldists[:, jsect]);
  # sort(
  #   @where(v_standard.matches, :tunit .== 4001),
  #   [:ttime, :munit, :possible]
  # )[1:10, :]
  # sort(
  #   @where(vnew.matches, :tunit .== 4001),
  #   [:ttime, :munit, :possible]
  # )[1:10, :]
  # hcat(ttimesr, idtreatr, idmatchr, mdistr, caldistsr')[224, :]'

  # setdiff(unid, unique(didsub))

  @inbounds for (j, mu) in enumerate(unid)

    # there are I (nrow robs) * J (length(unid)) combinations
    # l = (i - 1) * length(unid) + j

    idtreatr[j] = tu
    idmatchr[j] = mu
    ttimesr[j] = tt
    
    if (tu == mu)
      mdistr[j] = 0.0
      caldistsr[:, j] .= 0.0
    
    # the match (idj) must exist at or before tmn
    # save time by relying on [id, t] order (28 seconds w/ 2mill len loop)
    # elseif (minimum(dt[did .== idj]) <= tmn)
    else
      # data must exist far enough back
      # first time point for potential match must be less than tmnadj
      # that is, less-than-or-equal to (tt + mmin - fmax) or (the first time point in the data)
      # tmnadj becomes tmn if not assuming prior treatment
      if nopriortrt
        cpt = findfirst(didsub .== mu)
        if isnothing(cpt)
          possibler[j] = false
          continue
        else cp1 = dtsub[cpt] <= tmnadj;
        end
      elseif !nopriortrt
        cpt = findfirst(didsub .== mu)
        if isnothing(cpt)
          possibler[j] = false
          continue
        else
          cp1 = (dtsub[cpt] <= tmn);
        end
      end
      
      # may want check on length of data for a potential match, or even the treated, not just the range

      if !cp1
        possibler[j] = false
        continue
      else
        poss = (sum(dtrtsub[didsub .== mu]) == 0);
        possibler[j] = poss # no treatments in range
      end

      if poss # if possible do dist calcs

        # get unit indices for the match period
        # data must exist for whole match period

        # m-unit
        mux = getvarind(didsub, dtsub, mu, (tt + mmin), (tt + mmax));
          
        if (length(tux) == length(mux) == mlen)

          y = @views(cmatsub[tux, :]) # just sel and arrange columns up front
          x = @views(cmatsub[mux, :])

          dissum = [0.0] # Float64 is immutable
          calsums = zeros(Float64, length(calvars))
          
          matching_distancecalc_barrier!(
            dissum, calsums, y, x, tt,
            Σdict, cvkey, calvars,
            mmin, mlen,
            variancesonly
          )

          # avg dist between treated observation and possible match
          # over match period
          # must be mlen since we don't allow missingness in match period
          mdistr[j] = @views(dissum[1]) / mlen
          caldistsr[:, j] = calsums / mlen

        end
      end
    end
  end
  return idtreatr, idmatchr, ttimesr, possibler, mdistr, caldistsr
end

function matching_distancecalc_barrier!(
  dissum, calsums, y, x, ti,
  Σdict,
  cvkey, covariates,
  mmin, mlen,
  variancesonly
)

  @inbounds for l in 1:mlen # for each l in match period

    mp = (l - 1) + mmin

    if !variancesonly
      Σinv = LinearAlgebra.pinv(Σdict[(ti, mp)]);
    elseif variancesonly
      Σinv = LinearAlgebra.pinv(LinearAlgebra.Diagonal(Σdict[(ti, mp)]));
    end

    # slower about the same to get y and yy; but we do half that time and then half it again 20 times rather than full time 20 times
    # yy = @views(cmat[(dt .== ti + l) .& (did .== idi), :]);
    # xx = @views(cmat[(dt .== ti + l) .& (did .== idj), :]);
    
    yy = @views(y[l, :]);
    xx = @views(x[l, :]);

    diffs = yy - xx; # trt - cntrl
    # maha dist at t-l:
    if variancesonly
      dissum[1] += sqrt(transpose(diffs) * Σinv * diffs);
    end

    # this is not faster b/c of reshape?

    # xxt = @views(cmat2[:, (dt .== ti + l) .& (did .== idj)]);
    # yyt = @views(cmat2[:, (dt .== ti + l) .& (did .== idi)]);
    # @btime r = pairwise(Mahalanobis(Σinv), yyt, xxt);

    # caliper var calc here
    # calsums[m,i] rows on inner loop
    # NEW METHOD, NOT CONDITION
    # diag as 1d index for a given covariate
    dval = (cvkey .- 1) .* length(covariates) + cvkey
    σ_t = @view(Σinv[dval]);

    matching_distancecalc_caliper_barrier!(calsums, y, x, cvkey, σ_t, l)

  end
  return dissum, calsums
end

function matching_distancecalc_caliper_barrier!(calsums, y, x, cvkey, σ_t, l)

  for (h, m) in enumerate(cvkey) # key from names to cols in dat
    yyc = @views(y[l, m]);
    xxc = @views(x[l, m]);
    caldiff = yyc - xxc
    calsums[h] += sqrt(caldiff * @views(σ_t[h]) * caldiff)
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
      if length(ind) >= 1
        en = min(length(ind), ss) # include tobs with < setsize matches
        append!(indls, ind[1 : en])
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
function caliper(cdict::Dict{Symbol, Float64}, matches5::DataFrame)

  ntobsp = nrow(unique(matches5[!, [:ttime, :tunit]]));
  # need to add mdist

  keep = ones(Int64, nrow(matches5));
  
  for (k, v) in cdict
    if k != :mdist
      k = Symbol(String(k) * "_mdist")
      keep = keep .& (abs.(matches5[!, k]) .< v) # wasn't abs
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

  return matches5, lost, ntobs, ntobsp
end
