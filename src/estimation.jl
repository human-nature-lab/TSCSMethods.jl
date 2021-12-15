# estimation.jl

function outcomedict(dat, t, id, outcome)
  outdict = Dict{Tuple{Int64, Int64}, Float64}();
  @eachrow dat begin
    outdict[($t, $id)] = $outcome
  end;
return outdict
end

"""
    att!(results, M)

calculate the ATT.
"""
function att!(results, M)

  if "stratum" ∈ names(M)
    mgroup = [:stratum, :f]
  else
    mgroup = :f
  end

  # combine to get one value for each f, via sums over W^*_(i',t) * Y_(i',t)
  atts = @chain M begin
    groupby(mgroup)
    @combine(
      :att = sum(:wstar),
      :treatnum = sum(:treatev)
    ) # MUST DIVIDE BY N. TREATED UNITS
    @transform(:att = :att ./ :treatnum) # same order
  end

  append!(results, atts)

  return results
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

"""
faster version when length is known
"""
function countmembinput!(itr, d)
  for val in itr
      d[val] = get(d, val, 0) + 1
  end
  return d
end

function bootstrap(W, uid, observations, strata; iter = 500)
  # ADD TREAT EVE RESTRICTION MINIMUM?
  # >= one unit suffices

  c1 = "stratum" ∈ names(W);
  wgroup = if c1
    [:stratum, :f];
  else
    [:f];
  end
  
  uf = sort(unique(W[!, wgroup]), wgroup);
  flen = nrow(uf);

  uinfo = Vector{Vector{Int64}}(undef, flen);
  wout = Vector{Vector{Float64}}(undef, flen);
  trt = Vector{Vector{Int64}}(undef, flen);
  
  # actually, don't need to do this
  # strt = Vector{Vector{Int64}}(undef, len);
  
  G = groupby(W, wgroup);
  for (i, k) in enumerate(sort(keys(G)))
    g = get(G, k, 0);
    uinfo[i] = g[!, :unit];
    wout[i] = g[!, :wstar];
    trt[i] = g[!, :treatev];
    # if c1
    #   strt[i] = g[!, :stratum]; # this is assigning for each unit, need to ignore not-event items
    # end
  end

  boots = zeros(Float64, flen, iter);
  # boots = zeros(Float64, iter, length(uf));

  luid = length(uid);
  checkunitsets = if c1
    getcheckset(W, observations, strata)
  else
    getcheckset(W, observations)
  end
  
  _boot!(boots, uid, luid, wout, uinfo, trt, checkunitsets, iter);

  return boots
end

function getcheckset(W, observations, strata)

  cgroup = [:stratum, :f, :unit]
  checkset = @chain W[W.treatev, cgroup] begin
    groupby(cgroup[1:end-1])
    combine(
      :unit => Ref ∘ unique => :units,
    )
  end

  # treated units left after caliper application
  # we need at least one of these in the set
  trtleft = Dict{Int, Vector{Int}}();
  for s in unique(strata)
    trtleft[s] = [ob[2] for ob in observations[strata .== s]]
  end
  
  @eachrow! checkset begin
    :units = intersect(:units, trtleft[:stratum])
  end;

  # unique since we don't need to check for unitsets that are redundant
  # for each row, at least one unit (of that row) must exist in the sample

  # this is the set of treated units, such that at least one from
  # each row  must occur in the sample
  return unique(checkset.units)
end

function getcheckset(W, observations)

  cgroup = [:f, :unit];

  checkset = @chain W[W.treatev, cgroup] begin
    groupby(cgroup[1:end-1])
    combine(
      :unit => Ref ∘ unique => :units,
    )
  end

  # treated units left after caliper application
  # we need at least one of these in the set
  trtleft = [ob[2] for ob in observations]
  
  @eachrow! checkset begin
    :units = intersect(:units, trtleft)
  end;

  # unique since we don't need to check for unitsets that are redundant
  # for each row, at least one unit (of that row) must exist in the sample

  # this is the set of treated units, such that at least one from
  # each row  must occur in the sample
  return unique(checkset.units)
end

function _checkunits(units, samp)
  for unit in units # at least one must be true
    if unit ∈ samp # is the unit from a row of checkunitsets in the sample?
      return true
      break # if > 0 units are present, that row is satisfied by the sample
    end
  end
  return false
end

function samplepass(checkunitsets, samp)
  for units in checkunitsets
    if !_checkunits(units, samp)
      break
      return false
    end
  end
  return true
end

"""
    sample_check(uid, luid, checkset, uslen)

Sample the units, and ensure that there is a treated unit in the resample. If the data is stratified, ensure that there is at least one treated unit in each stratum.
"""
function sample_check(uid, luid, unitsets)
  samp = sample(uid, luid);

  while !(samplepass(unitsets, samp))
    samp = sample(uid, luid);
  end

  return samp
end

function _check(unitsets, samp)
  for units in unitsets
    for unit in units
      if unit ∈ samp
        break
      else
        # if last entry of a row does not yield true
        # (and haven't already skipped out to next row), then exit
        return false
      end
    end
  end
  return true # all rows have a treated unit present
end

"""
    _boot!(boots, uid, luid, wout, uinfo, trt, checkset, uslen, iter)

Inner function to bootstrap(), which actually executes the bootstrapping. N.B. that the seed should be set globally.
"""
function _boot!(boots, uid, luid, wout, uinfo, trt, checkunitsets, iter)
  @inbounds Threads.@threads for i in 1:iter
    reuid = countmemb(
      sample_check(uid, luid, checkunitsets),
      length(uid)
    );
    for c in eachindex(wout) # fs
      wo = wout[c]; ufo = uinfo[c]; to = trt[c]
      # n = stratified ? zeros(Int64, uslen) : zero(Int64)
      n = zero(Int64)
      for r in eachindex(wo)
        boots[c, i] += get(reuid, ufo[r], 0) * wo[r] # adding up to get a single estimate for an f
        n += get(reuid, ufo[r], 0) * to[r]
      end
      if n == 0
        println("lack of treated unit failure at c = ", c)
        boots[c, i] = NaN
        break
      end
      boots[c, i] = boots[c, i] / n
    end
  end
  return boots
end

"""
    bootinfo!(atts, boots; qtiles = [0.025, 0.5, 0.975])

Format the bootstrap matrix into the results dataframe. Assumes that att() has been executed.
"""
function bootinfo!(atts, boots; qtiles = [0.025, 0.5, 0.975])
  qnmes = Vector{Symbol}();
  for q in qtiles
    atts[!, :mean] = Vector{Float64}(undef, nrow(atts))
    qn = Symbol(string(q * 100) * "%");
    push!(qnmes, qn)
    atts[!, qn] = Vector{Float64}(undef, nrow(atts))
  end
  
  for (c, r) in enumerate(eachrow(atts))
    r[qnmes] = quantile(boots[c, :], qtiles)
    r[:mean] = mean(boots[c, :])
  end
  return atts
end

# import tscsmethods:@unpack,observationweights,att!,DataFrame,nrow
# using Accessors
# @reset calmodel.results = DataFrame()

"""
    estimate!(ccr::AbstractCICModel, dat; iterations = 500)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
  model::AbstractCICModel,
  dat::DataFrame;
  iterations = nothing,
  qtiles::Union{Float64, Vector{Float64}} = [0.025, 0.5, 0.975]
)

  if !isnothing(iterations)
    @reset model.iterations = iterations;
  end

  @unpack observations, ids, results, iterations = model;
  
  c1 = (length(observations) == 0);
  c2 = sum([isassigned(observations, i) for i in 1:length(observations)]) == 0
  if c1 | c2
    return "There are no matches."
  end

  W = observationweights(model, dat);
  results = att!(results, W);
  if typeof(model) <: AbstractCICModelStratified
    boots = bootstrap(W, ids, observations, model.strata; iter = iterations);
  else
    boots = bootstrap(W, ids, observations, nothing; iter = iterations);
  end
  bootinfo!(results, boots; qtiles = qtiles)
  
  return model
end
