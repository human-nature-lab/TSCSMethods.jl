# estimation_observationweights.jl

"""
treated observation weights
"""
function observationweights(model, dat)

  @unpack id, t, outcome, observations = model;

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

  Mloop!(M, Mt, observations, matches, ids, length(F));

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
    treatob = vcat(M[!, :treatob], M[!, :treatob]),
    f = vcat(M[!, :f], M[!, :f]),
    unit = vcat(M[!, :id], M[!, :id]),
    t = vcat(M[!, :t1], M[!, :t2]),
    w = vcat(M[!, :w1], M[!, :w2])
  );

  ## add strata here
  if c1
    stratmap = Dict(observations .=> model.strata)
    M[!, :stratum] = Vector{Int}(undef, nrow(M))
    @eachrow! M begin
      :stratum = stratmap[:treatob]
    end
  end

  M = @chain M begin
    groupby(mgroup) # group by t?
    @combine(:wstar = sum(:w))
    @subset(:wstar .!= 0)
  end

  return M
end

function get_anymus!(anymus, mus)
  for (i, r) in enumerate(eachrow(mus))
    if any(r)
      anymus[i] = true
    else 
      anymus[i] = false
    end
  end
  return anymus
end

function Mloop!(M, Mt, observations, matches, ids, Flen)
  for i in eachindex(observations)
    @unpack mus = matches[i]

    anymus = Vector{Bool}(undef, size(mus)[1]);
    get_anymus!(anymus, mus);

    fthere = Vector{Bool}(undef, Flen);
    getfunion!(fthere, mus);

    append!(
      M,
      DataFrame(
        treatob = fill(observations[i], sum(anymus)),
        id = ids[anymus],
        # fbs = matches[i].fs[matches[i].mus] # vector of vectors
        fbs = [r for r in eachrow(mus[anymus, :])]
      )
    );

    # treated obs.
    append!(
      Mt,
      DataFrame(
        treatob = observations[i],
        id = observations[i][2],
        fbs = [fthere]
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
