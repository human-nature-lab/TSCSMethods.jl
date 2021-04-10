using CSV, DataFrames, DataFramesMeta

"""
  getvarind()

get the indices of a unit over time time periode
"""
function getvarind(
    dat::DataFrame, unit, id, t, K)
  
    return findall(
      (dat[!, t] .>= minimum(K)) .& (dat[!, t] .<= maximum(K)) .&
      (dat[!, id] .== unit)
      )
  end
  
function getvarind(
  did::Vector{Int64},
  dt::Vector{Int64},
  unit,
  kmin, kmax
)

return findall(
  (dt .>= kmin) .& (dt .<= kmax) .&
  (did .== unit)
  )
end

function getvarind(
  didsub::SubArray,
  dtsub::SubArray,
  unit,
  kmin, kmax
)

return findall(
  (dtsub .>= kmin) .& (dtsub .<= kmax) .&
  (didsub .== unit)
  )
end

function getvarindset(
  did::Vector{Int64},
  dt::Vector{Int64},
  unit, K)

  a = (did .== unit);
  b1 = (dt .== K[1]);
  b2 = (dt .>= K[2]) .& (dt .<= maximum(K));

  return findall(a .& (b1 .| b2))
end

# this should replace getoutcomes!() above, just use enumerate and input K
# for a treated observation
function getvarvals!(
  outcol, tpoint, uti,
  dat, munit, id, t, covariate::Symbol;
  K = (tmin : tpoint) .+ uti)

  for (i, k) in enumerate(K)
    # inner loop should be rows                      
    ll = findall((dat[!, t] .== k) .& (dat[!, id] .== munit))
    # outcomemat[k - ut[i] - fmin + 1, i] = dat[outcome][ll][1]

    # what are the reasons that would cause an outcome to not exist? only date cutoff?

    if (length(ll) > 0) # if there is an outcome that exists
      outcol[i] = dat[covariate][ll][1]
    end

  end
  return outcol
end

function getvarvals!(
  outcol, tpoint, uti,
  dat, munit, id, t, covariate::Symbol,
  w;
  K = (tmin : tpoint) .+ uti)

  for (i, k) in enumerate(K)
    # inner loop should be rows                      
    ll = findall((dat[!, t] .== k) .& (dat[!, id] .== munit))
    # outcomemat[k - ut[i] - fmin + 1, i] = dat[outcome][ll][1]

    # what are the reasons that would cause an outcome to not exist? only date cutoff?

    if (length(ll) > 0) # if there is an outcome that exists
      outcol[i, j] += dat[covariate][ll][1] * w
    end
  end
  return outcol
end

  # check this again
function get_the_data(
    loc, c_list, covariates,
    fmin, fmax,
    treatment, outcome, id, t)
  
    dat = CSV.read(loc, DataFrame);
    # should filter to c_list + necessary variables
  
    dat = get_cols(dat, c_list, treatment, outcome, id, t);
  
    # get all treatment points
    treatment_points = get_all_treatment_points(dat, treatment);
  
    dat = remove_incomplete(
      dat, covariates,
      treatment_points, fmin, fmax, tmin, id, t);
  
    treatment_points = get_all_treatment_points(dat, treatment);
    return dat, treatment_points
  end
  
  function get_cols(dat, c_list, treatment, outcome, id, t)
    nec = [t, id, treatment, outcome] # no matching vars
    sels = vcat(nec, c_list)
    dat = select(dat, sels)
    return dat
  end
  
  # rely on dataframe order
  # this will need to be altered for multi-length days
  function get_all_treatment_points(dat, treatment)
    treatment_points = findall(dat[!, treatment] .== 1)
    return treatment_points
  end
  
  function get_all_treatment_points(dat_treatment)
    treatment_points = findall(dat_treatment .== 1)
    return treatment_points
  end
  
  function get_c_types(dat, covariates)

    desc = describe(dat[!, covariates])

    return c_types = desc.eltype
  end
  
  """
      remove_incomplete(dat, c_list, treatment_points, tmin, tmax, id, t)
  
  Remove, from the dataframe, all treated observations without
  long-enough / incomplete match periods. Removes
  observations for that treated unit from the match period up to
  the duration of the treatment period, on the idea
  that they are strictly not usable.
  
  NOTE may not need to remove the before-treatment periods
  """
  function remove_incomplete(
    dat, covariates, fmin, fmax, tmin, id, t, treatment, tpoint)

    treatment_points = get_all_treatment_points(dat, treatment);
  
    trt_obs = dat[treatment_points, [id, t]];
    trt_obs.trt_pts = treatment_points;
    tobsc = zeros(Int64, nrow(trt_obs));
    # tobsc[:] .= true
    trt_obs[!, :complete] = tobsc;
  
    remind = Array{Int64,1}(undef, 0)
    for i = eachindex(1:nrow(trt_obs))
      fps = trt_obs[!, id][i]
      trt_run = trt_obs[!, t][i] # treatment-admin point
      mtcht_min = trt_run + tpoint + tmin + 1
      mtcht_max = trt_run + tpoint
  
      rmin = trt_run + fmin # last day of fwindow
      rmax = trt_run + fmax # last day of fwindow
  
      #=
      check only the match period, allow missingness
      in the F window
      =#
      iind = getvarind(dat, fps, id, t, mtcht_min : mtcht_max);
  
      # check whether whole match period is present
      mc1 = length(iind) == length(mtcht_min:mtcht_max)
      #=
      check whether any matching covariates are missing
      at any point in matching period
      =#
      if mc1 == true
        datf = @view(dat[iind, covariates]);
        anymiss = 0
        for co = eachcol(datf)
          anymiss += any(ismissing(co)) * 1
        end
      else anymiss = 1
      end
      trt_obs.complete[i] = (mc1 & (anymiss == 0)) * 1
  
      #=
      if mtch covars and data not wholly present,
      remove the treated observation period
      for that tunit x ttime,
      remove observations from fmin (NOW trt_run, the admin point)
      to fmax
      =#
  
      if trt_obs.complete[i] == 0
        a = dat[!, id] .== fps;
        b = (dat[!, t] .>= trt_run);
        c = (dat[!, t] .<= rmax);
  
        remindi = findall(a .& b .& c);
  
        append!(remind, remindi)
      end
  
    end
    remind = sort(unique(remind))
  
    delete!(dat, remind)
    return dat
  end
  
MCDict = Dict{Tuple{Int64, Int64}, Int64}
function matchcounts(matches5)
  mm = @linq matches5 |>
    groupby([:ttime, :tunit])
  mm = combine(mm, nrow)
  
  mcd = MCDict()

  @eachrow mm begin
    mcd[(:ttime, :tunit)] = (:nrow - 1)
  end
  return mcd
end

"""
faster than countmap()
"""
function countmemb(itr)
  d = Dict{eltype(itr), Int}()
  for val in itr
      d[val] = get(d, val, 0) + 1
  end
  return d
end

"""
faster version when length is known
"""
function countmemb(itr, len::Int64)
  d = Dict{eltype(itr), Int64}()
  sizehint!(d, len)
  for val in itr
      d[val] = get(d, val, 0) + 1
  end
  return d
end

function namemodel(model::cicmodel)
  ttl = string(model.title)
  outc = string(model.outcome)
  cal = length(model.caliper) > 0 ? "cal" : ""
  strt = string(model.stratvar)

  return ttl * "_" * outc * "_" * strt * "_" * cal
end

"""
    assignquantile(stratvar, id, fulldat)

takes a dataset with the relevant variables and outputs a df
of quantile assignments for the chosen variable, unit-by-unit

this is used to create a variable with strata for restricted averaging in
att estimation

- treated_only: (default = true) calculates the quantiles only over the unique set of treated observations (t x unit id).
"""
function assignquantile(
  stratvar::Symbol,
  id::Symbol,
  t::Symbol,
  tmin::Int64,
  treatment::Symbol,
  fulldat::DataFrame;
  treatedonly = true
)

  if treatedonly == false
    ref = unique(
      fulldat[!, [id, stratvar]]);
    sq = quantile(ref[!, stratvar]);
  else
    c1 = fulldat[!, treatment] .== 1;
    c2 = fulldat[!, t] .> -tmin;
    ref = unique(
      @views(fulldat[c1 .& c2, [id, t, stratvar]]));
    sq = quantile(ref[!, stratvar]);
  end

  stratname = Symbol(String(stratvar) * "_stratum");
  ref[!, stratname] = zeros(Int64, nrow(ref));

  for i = eachindex(ref[!, stratname])
    sv = @view(ref[!, stratvar][i])[1]
    if sv >= sq[4]
      ref[!, stratname][i] = 4
    elseif (sv >= sq[3]) & (sv < sq[4])
      ref[!, stratname][i] = 3
    elseif (sv >= sq[2]) & (sv < sq[3])
      ref[!, stratname][i] = 2
    elseif sv < sq[2]
      ref[!, stratname][i] = 1
    end
  end
  return select(ref, [id, stratname]), sq;
end

function labelquantile(sq)
  pdkey = Dict{Int64, String}()
  for i in 2:5
    pdkey[i-1] = string(round(sq[i-1], digits = 2)) *
    " to " *
    string(round(sq[i], digits = 2))
  end
  return pdkey
end

"""
helper to sum the bootstrap distributions over some f range
"""
function estsum(bs)
  Z = sum(bs, dims = 2);
  Z = reshape(Z, size(Z)[1])
  return mean(Z), quantile(Z, [0.025, 0.5, 0.975])
end

"""
gives separate sum for each week, respecting days of week, for f window
note that, depending on f window, weeks may have different numbers of days
"""
function weeklyattci(bs, atts, startday::Int64, fmin, fmax)
  # create vector denoting week assignments over cols of bs
  dow = fld.(collect(startday .+ (fmin:fmax)), 7) # from week of primary (0) on
  udow = unique(dow)

  bss = Matrix{Float64}(undef, size(bs)[1], length(udow));
  attsums = Vector{Float64}(undef, length(udow));

  for (i, e) in enumerate(udow)
    lx = findall(dow .== e)
    bss[:, i] = sum(@views(bs[:, lx]), dims = 2)
    attsums[i] = sum(@views(atts[lx]))
  end

  return bss, attsums, udow
end

function processattsum(bss, attsums, udow)
  qtm = Matrix{Float64}(undef, length(attsums), 3)
  for i = eachindex(attsums)
    qtm[i, :] = quantile(bss[:, i], [0.025, 0.5, 0.975])
  end

  sumres = DataFrame(
    week = udow,
    attsum = attsums,
    lwer = qtm[:, 1],
    med = qtm[:, 2],
    uper = qtm[:, 3]
  )
  return sumres
end

function weeklyatt(bs::Matrix{Float64}, res, startday)
  
  res = sort(res, [:f])

  bss, attsums, udow = weeklyattci(
    bs,
    res.att,
    startday,
    minimum(res.f),
    maximum(res.f)
  )
  sumres = processattsum(bss, attsums, udow)
  return sumres
end

function weeklyatt(
  bs::Dict{Int64, Matrix{Float64}},
  res,
  startday # day of the week as Int
)
  res = sort(res, [:stratum, :f])

  SumRes = DataFrame(
    [Int64, Int64, Float64, Float64, Float64, Float64],
    [:stratum, :week, :attsum, :lwer, :med, :uper],
    0
  )

  for (k, v) in bs
    resk = @views(res[res.stratum .== k, :])
    bss, attsums, udow = weeklyattci(
      v,
      resk[!, :att],
      startday,
      minimum(resk.f),
      maximum(resk.f)
    )
    sumres = processattsum(bss, attsums, udow)
    sumres[!, :stratum] .= k
    append!(SumRes, sumres)
  end
  return SumRes
end

function plot_att_sum(sumres; savename = "")
  wmin = minimum(sumres.week)

  att_plt = plot(
    sort(sumres, [:week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      Geom.point,
      Geom.vline(style = :dot, color = "black", size = [0.2mm]),
      Geom.errorbar,
      Guide.title("ATT (weekly sum)"),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Coord.Cartesian(ymin = wmin)
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end


function plot_att_sum(sumres, stratvar::Symbol; savename = "")
  wmin = minimum(sumres.week)

  att_plt = plot(
      sort(sumres, [:stratum, :week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      xgroup = :stratum,
      Guide.title("ATT (weekly sum)" * " by " * string(stratvar)),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Geom.subplot_grid(
        Geom.point,
        Geom.vline(style = :dot, color = "black", size = [0.2mm]),
        Geom.errorbar,
        Coord.Cartesian(ymin = wmin)
      )
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end

function plot_att_sum(sumres, stratvar::Symbol, stratdict; savename = "")
  wmin = minimum(sumres.week)

  sumres.svlabel = getindex.(Ref(stratdict), sumres[!, :stratum])  

  att_plt = plot(
      sort(sumres, [:stratum, :week]),
      xintercept = [0],
      y = :week,
      x = :attsum,
      xmin = :lwer, xmax = :uper,
      xgroup = :svlabel,
      Guide.title("ATT (weekly sum)" * " by " * string(stratvar)),
      Guide.ylabel("weeks after the primary"),
      Guide.xlabel("estimate"),
      Geom.subplot_grid(
        Geom.point,
        Geom.vline(style = :dot, color = "black", size = [0.2mm]),
        Geom.errorbar,
        Coord.Cartesian(ymin = wmin)
      )
  )

  if length(savename) > 0
    draw(PNG(savename, 9inch, 5inch), att_plt)
  end
  return att_plt
end

# SumRes = weeklyatt(
#   bs, res, Dates.Tuesday
# )

# plot_att_sum(
#   SumRes,
#   :running;
#   savename = pltl * namemodel(v_cc_cal) * "_att_wk" * ".png"
# )

"""
      datamake(covariates, varset, fmin, fmax, tmin,
      id, tpoint, treatment
      )

use this function when the original data is mutated in some model run...
"""
function datamake(
  covariates, varset,
  fmin, fmax, tmin,
  id, tpoint, treatment
)

  @load "../covid-19-data/data/cvd_dat.jld2"

  dat = @where(dat, :running .>= 0);

  dat = remove_incomplete(
    dat, covariates, fmin, fmax, tmin, id, t, treatment, tpoint);

  select!(dat, varset);
  return dat
end


import CSV.read

"""
      datamakecc(covariates, varset, fmin, fmax, tmin,
      id, tpoint, treatment
      )

use this function when the original data is mutated in some model run...
"""
function datamakecc(
  covariates, varset,
  fmin, fmax, mmin,
  id, tpoint, treatment
)

  @load "../covid-19-data/data/cvd_dat.jld2"

  dat = @where(dat, :running .>= 0);

  select!(dat, varset)
  
  evac = read("data/hurricane_evacuation_orders.csv", DataFrame);
  select!(evac, [:Hurricane, :Date, :Place, :fips, :Mandatory])

  evac = @where(evac, :fips .!= "NA");

  evac.fips = parse.(Int64, evac.fips);

  evac = @where(evac, ismissing.(:Date) .== false);

  # move date back
  evac.Date = evac.Date .- Dates.Day(5)

  evac.evac = ones(Int64, nrow(evac))

  dat = leftjoin(dat, evac, on = [:date => :Date, :fips]);
  dat.evac[ismissing.(dat.evac)] .= 0;

  dat = remove_incomplete(
    dat, covariates, fmin, fmax, mmin, id, t, treatment, tpoint);

  codist = CSV.read("/Users/emf/Library/Mobile Documents/com~apple~CloudDocs/Yale/yale research/COVID19/covid-19-data/data/sf12010countydistance500miles.csv", DataFrame)
  
  exc = @where( # exclude county2
    codist,
    :mi_to_county .< 250,
    :county1 .∈ Ref(evac.fips),
    :county2 .∉ Ref(evac.fips)
    );

  destfips = [48453, 48029];

  dat = @where(
    dat,
    :fips .∉ Ref(setdiff(exc.county2, destfips))
  );

  c1 = dat.fips .∈ Ref([48453, 48029]);
  c2 = dat.date .== (Dates.Date("2020-08-25") - Dates.Day(5));
  dat.evac[c1 .& c2, :] .= 1

  return dat
end