# information.jl
# various functions that give information about the models, etc.

"""
    unitcounts(model::VeryAbstractCICModel)

Count treated and matched units for each time period F.

# Arguments
- `model`: Fitted TSCSMethods model

# Returns
- For non-stratified models: `(treated_counts, matched_counts)`
  - `treated_counts`: Vector{Int} - number of treated units per F period
  - `matched_counts`: Vector{Int} - total matched units available per F period
- For stratified models: `(treated_by_stratum, matched_by_stratum)`
  - `treated_by_stratum`: Dict{Int, Vector{Int}} - treated counts by stratum
  - `matched_by_stratum`: Dict{Int, Vector{Int}} - matched counts by stratum

# Description
Counts how many treated units remain eligible and how many matched control units 
are available for each time period in the post-treatment window F. Essential for 
assessing statistical power and balance of the matching procedure.

For stratified models, counts are organized by stratum to enable stratum-specific
power analysis and balance assessment.

# Examples
```julia
# Non-stratified model
treated_counts, matched_counts = unitcounts(model)
println("F=1: \$(treated_counts[1]) treated, \$(matched_counts[1]) matched")

# Stratified model  
treated_strata, matched_strata = unitcounts(stratified_model)
for (stratum, counts) in treated_strata
    println("Stratum \$stratum: \$(counts[1]) treated units in F=1")
end
```
"""
function unitcounts(model::VeryAbstractCICModel)
    # Input validation
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    if length(model.matches) != length(model.observations)
        throw(ArgumentError("Mismatch between observations and matches"))
    end
    
    # Extract match availability for each treated observation
    # Each element is a vector showing matches available per F period
    match_availability_per_obs = [vec(sum(match.eligible_matches, dims=1)) 
                                 for match in model.matches]
    
    n_periods = length(model.F)
    is_stratified = typeof(model) <: AbstractCICModelStratified
    
    if !is_stratified
        # Non-stratified: aggregate across all treated observations
        treated_counts = zeros(Int, n_periods)
        _count_treated_units!(treated_counts, match_availability_per_obs)
        
        # Total matched units available per period
        matched_counts = sum(match_availability_per_obs)
        return treated_counts, matched_counts
    else
        # Stratified: organize by strata
        unique_strata = sort(unique(model.strata))
        treated_by_stratum = Dict{Int, Vector{Int}}()
        matched_by_stratum = Dict{Int, Vector{Int}}()
        
        for stratum in unique_strata
            stratum_mask = model.strata .== stratum
            stratum_match_availability = match_availability_per_obs[stratum_mask]
            
            # Count treated units in this stratum
            treated_by_stratum[stratum] = zeros(Int, n_periods)
            _count_treated_units!(treated_by_stratum[stratum], stratum_match_availability)
            
            # Count matched units available in this stratum (BUG FIX: was Us[i] = ...)
            matched_by_stratum[stratum] = sum(stratum_match_availability)
        end
        
        return treated_by_stratum, matched_by_stratum
    end
end

"""
    _count_treated_units!(treated_counts, match_availability)

Internal helper to count treated units with available matches for each F period.

Modifies `treated_counts` in place by incrementing the count for each period
where a treated unit has at least one available match.
"""
function _count_treated_units!(treated_counts::Vector{Int}, match_availability::Vector{Vector{Int}})
    for unit_matches in match_availability
        for (period_idx, available_matches) in enumerate(unit_matches)
            if available_matches > 0
                treated_counts[period_idx] += 1
            end
        end
    end
    return nothing
end

"""
    eligibility(model)

"eligibility" for units, over all treated units => num. times a particular unit
is eligible to be a match, for each F in the outcome window.
"""
function eligibility(model)
  (; matches) = model
  eligibles = similar(matches[1].eligible_matches) .* 0;
  
  for objet in matches
    eligibles += objet.eligible_matches
  end
  return eligibles
end

"""
    whentreated(fips, model)

Input a unit id, and see when (if) a treatment occured.
"""
function whentreated(fips, model)

  (; observations) = model;

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
  (; matches, observations) = model;

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

  return sort(matches[cnt].match_rankings)
end

function matchinfo(model::Union{CIC, CICStratified}; maxrank = 5)
  (; observations, matches, F) = model;

  mf = DataFrame(
    timetreated = Int[],
    treatedunit = Int[],
    f = Int[],
    matchunits = Vector{Vector{Int}}(),
    ranks = Vector{Vector{Int}}()
  );

  for (i, treatob) in enumerate(observations)
    for f in F
      outcome_period_index = f - minimum(F) + 1
      vlen = min(length(matches[i].match_rankings[outcome_period_index]), maxrank);
      push!(
        mf,
        [
          treatob[1],
          treatob[2],
          f,
          matches[i].match_rankings[outcome_period_index][1:vlen],
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

  for r in 1:length(mf.matchunits)
    for j in 1:length(mf.matchunits[r])
      mf.matchunits[r][j] = model.ids[mf.matchunits[r][j]]
    end
  end

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
      outcome_period_index = f - minimum(rc.F) + 1      
      mus = model.ids[rc.matches[i].eligible_matches[:,outcome_period_index]];
      mus = mus[1:min(length(mus), maxrank)];
      fullmus = model.ids[model.matches[cnt].match_rankings[outcome_period_index]];
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