# new_estimation.jl

function outcomedict(dat, t, id, outcome)
  outdict = Dict{Tuple{Int64, Int64}, Float64}();
  @eachrow dat begin
    outdict[($t, $id)] = $outcome
  end;
return outdict
end

# obs = matches[1]

# hcat(obs.mus[obs.mus], obs.fs[obs.mus])
# sum(obs.fs[obs.mus]) # treatnums

# fsout = convert(
#   Vector{Vector{Float64}}, similar(matches[1].fs[matches[1].mus])
# );
# # w for matches
# for i in 1:length(matches[1].mus[matches[1].mus])
#   fsout[i] = obs.fs[obs.mus][i] .* 1 ./ sum(obs.fs[obs.mus])
# end

function Mloop!(M, Mt, observations, matches)
  for i in eachindex(observations)
    # matches
    append!(
      M,
      DataFrame(
        treatob = fill(observations[i], sum(matches[i].mus)),
        id = ids[matches[i].mus],
        fbs = matches[i].fs[matches[i].mus]
      )
    );
    
    # treated obs.
    append!(
      Mt,
      DataFrame(
        treatob = observations[i],
        id = observations[i][2],
        fbs = [sum(matches[i].fs[matches[i].mus]) .> 0] # check whether every f is included
      )
    );
  end  
  return M, Mt
end

function Mprocess(Mt, F, Fmin)
  Mt[!, :f] = Vector{Vector{Int}}(undef, nrow(Mt));
  @inbounds Threads.@threads for i in 1:nrow(Mt)
    Mt.f[i] = Mt.fbs[i] .* collect(eachindex(F)) # index since F may be 0
  end

  select!(Mt, Not(:fbs));
  Mt = flatten(Mt, :f);

  Mt = Mt[Mt.f .!= 0, :]

  Mt[!, :f] = Mt[!, :f] .+ (Fmin - 1);
  return Mt
end

function Mprep(model)

  if typeof(model) <: AbstractCICModelStratified
    c1 = true
    mgroup = [:stratum, :f, :t, :unit]
  else
    c1 = false
    mgroup = [:f, :t, :unit]
  end

  # construct two dataframes
  # - matches
  # - treated units
  # in each case we have the unique set, at each f

  @unpack observations, matches, ids, F, reference = model;
  Fmin = minimum(F);

  M = DataFrame(
    treatob = Tuple{Int, Int}[],
    id = Int[],
    fbs = Vector{Bool}[]
  );

  Mt = similar(M, 0);

  Mloop!(M, Mt, observations, matches);

  Mt = Mprocess(Mt, F, Fmin);
  M = Mprocess(M, F, Fmin);

  Mt[!, :t2] = Vector{Int}(undef, nrow(Mt));
  Mt[!, :t1] = Vector{Int}(undef, nrow(Mt));
  Mt[!, :w2] = Vector{Float64}(undef, nrow(Mt));
  Mt[!, :w1] = Vector{Float64}(undef, nrow(Mt));
  @eachrow! Mt begin
    :t1 = :treatob[1] + :f;
    :t2 = :treatob[1] + reference;
    :w2 = 1.0;
    :w1 = -1.0;
  end;

  M = @chain M begin
    groupby([:treatob, :f]) # get stratum for free
    transform(nrow => :w2)
  end

  M[!, :t2] = Vector{Int}(undef, nrow(M));
  M[!, :t1] = Vector{Int}(undef, nrow(M));
  M[!, :w1] = Vector{Float64}(undef, nrow(M));
  M[!, :w2] = convert(Vector{Float64}, M[!, :w2]);
  @eachrow! M begin
    :t1 = :treatob[1] + :f;
    :t2 = :treatob[1] + reference;
    :w2 = -1.0 / :w2;
    :w1 = -:w2;
  end
  
  append!(M, Mt)

  # at (i,t) level
  # if c1
  #   M = DataFrame(
  #     stratum = vcat(M[!, :stratum], M[!, :stratum]),
  #     f = vcat(M[!, :f], M[!, :f]),
  #     unit = vcat(M[!, :matchunit], M[!, :matchunit]),
  #     t = vcat(M[!, :t1], M[!, :t2]),
  #     w = vcat(M[!, :w1], M[!, :w2])
  #   );
  # else
  #   M = DataFrame(
  #     f = vcat(M[!, :f], M[!, :f]),
  #     unit = vcat(M[!, :id], M[!, :id]),
  #     t = vcat(M[!, :t1], M[!, :t2]),
  #     w = vcat(M[!, :w1], M[!, :w2])
  #   );
  # end

  M = DataFrame(
    f = vcat(M[!, :f], M[!, :f]),
    unit = vcat(M[!, :id], M[!, :id]),
    t = vcat(M[!, :t1], M[!, :t2]),
    w = vcat(M[!, :w1], M[!, :w2])
  );

  ## add strata here
  if c1
    stratamap = Dict(observations .=> model.strata)
    M[!, :stratum] = Vector{Int}(undef, nrow(M))
    @eachrow! M begin
      :stratum = stratmap[:treatob]
    end
  end

  M = @chain M begin
    groupby(mgroup) # group by t?
    @combine(wstar = sum(:w))
    @subset(:wstar .!= 0)
  end

  return M
end

# treated observation weights
function observationweights(model, dat)

  @unpack id, t, outcome = model;

  M = Mprep(model);

  outdict = outcomedict(dat, t, id, outcome);
  
  trtdict = Dict{Int64, Bool}();
  for obs in observations
    trtdict[obs[2]] = true
  end
  
  M[!, :treatev] .= false
  
  @eachrow! M begin
    if :wstar == 1
      :treatev = get(trtdict, :unit, false)
    end
  end

  M = @eachrow M begin
    :wstar = :wstar * outdict[(:t, :unit)]
  end

  
  # M[!, :treatev] .= false
  # @eachrow! M begin
  #   :treatev = get(trtdict, (:t - ccr.reference, :unit), false)
  # end

  return M
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
      att = sum(:wstar),
      treatnum = sum(:treatev)
    ) # MUST DIVIDE BY N. TREATED UNITS
    @transform(att = :att ./ :treatnum) # same order
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
  model::AbstractCICModel,
  dat::DataFrame;
  iter::Int = 500,
  qtiles::Union{Float64, Vector{Float64}} = [0.025, 0.5, 0.975]
)

  @unpack observations, ids, results = model;
  
  c1 = (length(observations) == 0);
  c2 = sum([isassigned(observations, i) for i in 1:length(observations)]) == 0
  if c1 | c2
    return "There are no matches."
  end

  W = observationweights(model, dat);
  results = att!(results, W);
  boots = bootstrap(W, ids; iter = iter);
  bootinfo!(results, boots; qtiles = qtiles)
  
  return model
end
