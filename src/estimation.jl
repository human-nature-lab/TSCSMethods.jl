# new_estimation.jl

function outcomedict(dat, t, id, outcome)
  outdict = Dict{Tuple{Int64, Int64}, Float64}();
  @eachrow dat begin
    outdict[($t, $id)] = $outcome
  end;
return outdict
end

function Mprep(matches, reference; mgroup = [:f, :t, :unit])

  c1 = :stratum ∈ mgroup

  if c1
    mgroup2 = [:stratum, :treattime, :treatunit, :matchunit, :f];
    mgroup3 = [:stratum, :treattime, :treatunit, :f];
  else
    mgroup2 = [:treattime, :treatunit, :matchunit, :f];
    mgroup3 = [:treattime, :treatunit, :f];
  end
  
  M = matches[:, mgroup2];
  
  Mt = begin
    Mt = unique(M[!, mgroup3]);
    Mt[!, :matchunit] .= Mt[!, :treatunit];
    Mt[!, :w2] .= 1
    Mt[!, :t2] = Mt[!, :treattime] + Mt[!, :f]
    Mt[!, :w1] .= -Mt.w2
    Mt[!, :t1] = Mt[!, :treattime] .+ reference
    Mt;
  end;

  M = @chain M begin
    groupby(mgroup3)
    transform(nrow => :w2)
    @transform(w2 = -1 ./ :w2, t2 = :treattime + :f)
    @transform(w1 = -:w2, t1 = :treattime .+ reference)
  end
  
  append!(M, Mt)

  # at (i,t) level
  if c1
    M = DataFrame(
      stratum = vcat(M[!, :stratum], M[!, :stratum]),
      f = vcat(M[!, :f], M[!, :f]),
      unit = vcat(M[!, :matchunit], M[!, :matchunit]),
      t = vcat(M[!, :t1], M[!, :t2]),
      w = vcat(M[!, :w1], M[!, :w2])
    );
  else
    M = DataFrame(
      f = vcat(M[!, :f], M[!, :f]),
      unit = vcat(M[!, :matchunit], M[!, :matchunit]),
      t = vcat(M[!, :t1], M[!, :t2]),
      w = vcat(M[!, :w1], M[!, :w2])
    );
  end

  M = @chain M begin
    groupby(mgroup) # group by t?
    @combine(wstar = sum(:w))
    @subset(:wstar .!= 0)
  end

  return M
end

# treated observation weights
function observationweights(ccr, dat)

  if ccr.stratifier == Symbol("")
    M = Mprep(ccr.matches, ccr.reference);
  else
    M = Mprep(ccr.matches, ccr.reference; mgroup = [:stratum, :f, :t, :unit]);
  end

  outdict = outcomedict(dat, ccr.t, ccr.id, ccr.outcome);

  M = @eachrow M begin
    :wstar = :wstar * outdict[(:t, :unit)]
  end

  trtdict = begin
    trtdict = Dict{Tuple{Int64, Int64}, Bool}();
    tobs = ccr.matches[!, [:treattime, :treatunit]];
    unique(tobs, [:treattime, :treatunit])
    @eachrow tobs begin
      trtdict[(:treattime, :treatunit)] = true
    end
    trtdict
  end;
  
  M[!, :treatev] .= false
  @eachrow! M begin
    :treatev = get(trtdict, (:t - ccr.reference, :unit), false)
  end

  return M
end

"""
    att(M, outdict, trtdict)

calculate the ATT.
"""
function att(M)

  if "stratum" ∈ names(M)
    mgroup = [:stratum, :f]
  else
    mgroup = :f
  end

  # combine to get one value for each f, via sums over W^*_(i',t) * Y_(i',t)
  atts = @chain M begin
    groupby(mgroup)
    @combine(
      att = sum(:wstar),
      treatnum = sum(:treatev)
    ) # MUST DIVIDE BY N. TREATED UNITS
    @transform(att = :att ./ :treatnum) # same order
  end

  return atts
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

function bootstrap(W, uid; iter = 500)
  # ADD TREAT EVE RESTRICTION MINIMUM?
  # >= one unit suffices

  c1 = "stratum" ∈ names(W);
  wgroup = if c1
    [:stratum, :f];
  else
    [:f];
  end
  
  uf = sort(unique(W[!, wgroup]), wgroup);
  len = nrow(uf);

  uinfo = Vector{Vector{Int64}}(undef, len);
  wout = Vector{Vector{Float64}}(undef, len);
  trt = Vector{Vector{Int64}}(undef, len);
  
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

  boots = zeros(Float64, len, iter);
  # boots = zeros(Float64, iter, length(uf));

  #=
  if !(us == 1:uslen)
    error("strata definition problem")
  end
  =#

  luid = length(uid);
  if c1
    cgroup = [:stratum, :f, :unit]
    us = unique(uf.stratum);
    uslen = length(us);
  else
    cgroup = [:f, :unit]
    uslen = 1
  end
  
  # new checkset
  checkunitsets = begin
    checkset = @chain W[W.treatev, cgroup] begin
      groupby(cgroup[1:end-1])
      combine(
        :unit => Ref ∘ unique => :units,
      )
    end
    # unique since we don't need to check for unitsets that are redundant
    # for each row, at least one unit (of that row) must exist in the sample
    unique(checkset.units)
  end
  
  _boot!(boots, uid, luid, wout, uinfo, trt, checkunitsets, iter);

  return boots
end

"""
    sample_check(uid, luid, checkset, uslen)

Sample the units, and ensure that there is a treated unit in the resample. If the data is stratified, ensure that there is at least one treated unit in each stratum.
"""
function sample_check(uid, luid, unitsets)
  samp = sample(uid, luid);

  while !(_check(unitsets, samp))
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
      for r in 1:length(wo)
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

"""
    estimate!(ccr::AbstractCICModel, dat; iter = 500)

Perform ATT estimation, with bootstrapped CIs.
"""
function estimate!(
  ccr::AbstractCICModel,
  dat::DataFrame;
  iter::Int = 500,
  qtiles::Union{Float64, Vector{Float64}} = [0.025, 0.5, 0.975]
)

  if nrow(ccr.matches) == 0
    return "There are no matches."
  end

  uid = unique(dat[!, ccr.id])
  W = observationweights(ccr, dat);
  ccr.results = att(W);
  boots = bootstrap(W, uid; iter = iter);
  bootinfo!(ccr.results, boots; qtiles = qtiles)
  
  return ccr.results
end
