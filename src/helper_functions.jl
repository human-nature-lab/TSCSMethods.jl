using CSV, DataFrames, DataFramesMeta

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
    unit, K)
  
    return findall(
      (dt .>= minimum(K)) .& (dt .<= maximum(K)) .&
      (did .== unit)
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
    dat, covariates, fmin, fmax, tmin, id, t, treatment)

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
  