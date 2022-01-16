# information.jl
# various functions that give information about the models, etc.

"""
    eligibility(model)

"eligibility" for units, over all treated units => num. times a particular unit
is eligible to be a match, for each F in the outcome window.
"""
function eligibility(model)
  @unpack matches = model
  eligibles = similar(matches[1].mus) .* 0;
  
  for objet in matches
    eligibles += objet.mus
  end
  return eligibles
end

"""
    whentreated(fips, model)

Input a unit id, and see when (if) a treatment occured.
"""
function whentreated(fips, model)

  @unpack matches, observations = model;

  loc = Int[];
  cnt = 0
  found = false
  for ob in observations
    cnt += 1
    if fips == ob[2]
      found = true
      push!(loc, cnt)
    end
  end

  if !found
    return "treated unit not found"
  end

  return observations[loc]
end

"""
    showmatches(model, treatob)
  
Show the ranked (best to worst) eligible matches to a chosen
treated observation.
"""
function showmatches(model, treatob)
  @unpack matches, observations = model;

  cnt = 0
  found = false
  for ob in observations
    cnt += 1
    if treatob == ob
      found = true
      break
    end
  end

  if !found
    return "treated observation not found"
  end

  return sort(matches[cnt].ranks)
end

function matchinfo(model::Union{CIC, CICStratified}; maxrank = 5)
  @unpack observations, matches, ids, F = model;
  mf = DataFrame(
    timetreated = Int[],
    treatedunit = Int[],
    f = Int[],
    matchunits = Vector{Vector{Int}}(),
    ranks = Vector{Vector{Int}}()
  );

  for (i, treatob) in enumerate(observations)
    for f in F
      φ = f - minimum(F) + 1
      vlen = min(length(matches[i].ranks[φ]), maxrank);
      push!(
        mf,
        [
          treatob[1],
          treatob[2],
          f,
          matches[i].ranks[φ][1:vlen], # ids[matches[i].mus[:,φ]],
          (1:vlen)[1:vlen],
        ]
      );
    end
  end

  # remove fs where there are no matches
  keepidx = Int[];
  for (i, r) in enumerate(eachrow(mf))
    if length(r[:matchunits]) > 0
      append!(keepidx, i)
    end
  end

  mf = mf[keepidx, :];

  sort!(mf, [:timetreated, :treatedunit])

  return mf
end

function matchinfo(rc, model; maxrank = 5)

  mf = DataFrame(
    timetreated = Int[],
    treatedunit = Int[],
    f = Int[],
    matchunits = Vector{Vector{Int}}(),
    ranks = Vector{Vector{Int}}()
  );

  for (i, treatob) in enumerate(rc.observations)

    # find in original model
    # since indices may differ
    cnt = 0
    found = false
    for ob in model.observations
      cnt += 1
      if treatob == ob
        found = true
        break
      end
    end

    if !found
      return error("observation not found in original model")
    end
    
    for f in rc.F
      φ = f - minimum(rc.F) + 1      
      mus = model.ids[rc.matches[i].mus[:,φ]];
      mus = mus[1:min(length(mus), maxrank)];
      fullmus = model.ids[model.matches[cnt].ranks[φ]];
      rnks = [findfirst(x -> x == mu, fullmus) for mu in mus];
    
      rnks = rnks[sortperm(rnks)];
      mus = mus[sortperm(rnks)];
  
      push!(
        mf,
        [
          treatob[1],
          treatob[2],
          f,
          mus,
          rnks
        ]
      );
    end
  end
  
  # remove fs where there are no matches
  keepidx = Int[];
  for (i, r) in enumerate(eachrow(mf))
    if length(r[:matchunits]) > 0
      append!(keepidx, i)
    end
  end

  mf = mf[keepidx, :];

  sort!(mf, [:timetreated, :treatedunit])

  return mf
end

"""
    obsinfo(
      rcinfo, dat, vars;
      fullmodobs = nothing,
      t = :running, id = :fips
    )

Get the requested info for the treated units. If the var is timevarying,
it will get it from the time of the treatment event.
"""
function obsinfo(
  rcinfo, dat, vars;
  fullmodobs = nothing,
  t = :running, id = :fips
)

  obinfo = unique(rcinfo[!, [:timetreated, :treatedunit]]);

  if !isnothing(fullmodobs)
    obinfo[!, :removed] .= false
    calobs = [(x, y) for (x, y) in zip(obinfo[!, :timetreated], obinfo[!, :treatedunit])];
    
    for (x,y) in setdiff(fullmodobs, calobs)
      push!(
        obinfo,
        [x, y, true])
    end
  end


  for var in vars
    obinfo[!, var] = Vector{eltype(dat[!, var])}(undef, nrow(obinfo))
  end

  for r in eachrow(obinfo)
    for var in vars
      c1 = (dat[!, t] .== r[:timetreated]) .& (dat[!, id] .== r[:treatedunit]);
      
      r[var] = dat[c1, var][1];
    end
  end

  if !isnothing(fullmodobs)
    sort!(obinfo, [:removed, :timetreated, :treatedunit])
  else
    sort!(obinfo, [:timetreated, :treatedunit])
  end

  return obinfo
end