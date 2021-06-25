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
helper to sum the bootstrap distributions over some f range
"""
function estsum(bs)
  Z = sum(bs, dims = 2);
  Z = reshape(Z, size(Z)[1])
  return mean(Z), quantile(Z, [0.025, 0.5, 0.975])
end

"""
https://stackoverflow.com/questions/40694157/how-to-reverse-a-dictionary-in-julia
"""
function invert_dict(dict, warning::Bool = false)
  vals = collect(values(dict))
  dict_length = length(unique(vals))

  if dict_length < length(dict)
      if warning
          warn("Keys/Vals are not one-to-one")
      end 

      linked_list = Array[]

      for i in vals 
          push!(linked_list,[])
      end 

      new_dict = Dict(zip(vals, linked_list))

      for (key,val) in dict 
          push!(new_dict[val],key)
      end
  else
      key = collect(keys(dict))

      counter = 0
      for (k,v) in dict 
          counter += 1
          vals[counter] = v
          key[counter] = k
      end
      new_dict = Dict(zip(vals, key))
  end 

  return new_dict
end

function namemodel(model::cicmodel)
  ttl = string(model.title)
  outc = string(model.outcome)
  cal = length(model.caliper) > 0 ? "cal" : ""
  strt = string(model.stratvar)

  return ttl * "_" * outc * "_" * strt * "_" * cal
end

function mkDataFrame(cts)
  df = DataFrame()
  for (k, v) in cts
    df[!, k] = v
  end
  return df
end