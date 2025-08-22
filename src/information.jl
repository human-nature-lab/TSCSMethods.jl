# information.jl
# various functions that give information about the models, etc.

"""
    unitcounts(m)

count the total number of treated and the number of matched units for each F.
"""
function unitcounts(m)
    
    X = Vector{Vector{Int}}(undef, 0);

    Flen = length(m.F)
    
    for i in eachindex(m.matches)
        push!(
            X,
            vec(sum(m.matches[i].eligible_matches, dims = 1))
        )
    end

    strat = any(
        [
            typeof(m) == x for x in [
                CICStratified, RefinedCICStratified, RefinedCaliperCICStratified
            ]
        ]
    );

    if !strat
        # compute the number of treated units left (for each F)
        Y = zeros(Int, Flen)
        _treatedcount!(Y, X)
        
        # compute the total number of match units (for each F)
        U = sum(X)
        return Y, U
    else
        S = sort(unique(m.strata));
        # with stratification
        Ys = Dict{Int, Vector{Int}}();
        Us = Dict{Int, Vector{Int}}();

        # compute the number of treated units left (for each F)
        for (i, s) in enumerate(S)
            Ys[s] = zeros(Int, Flen)
            Us[s] = zeros(Int, Flen)
            c1 = m.strata .== s;
            _treatedcount!(Ys[s], X[c1])
        
            # compute the total number of match units (for each F)
            Us[i] = sum(X[c1])
        end
        return Ys, Us
    end
end

function _treatedcount!(Y, X)
    for x in X
        for (i, xi) in enumerate(x)
            if xi > 0
                Y[i] += 1
            end
        end
    end
end

"""
    eligibility(model::VeryAbstractCICModel)

Calculate unit eligibility across all treated observations.

# Arguments  
- `model`: Fitted TSCSMethods model

# Returns
- `Matrix{Int}`: Eligibility matrix (units × time periods F)

# Description
For each potential control unit, calculates how many times it is eligible
to serve as a match across all treated observations and time periods.
Higher eligibility indicates units that are consistently good potential matches
for multiple treated units.

Useful for identifying:
- Units that are consistently good matches across multiple treated observations
- Potential issues with match quality or overlap
- Whether certain units dominate the matching process
- Balance in the matching pool

# Examples
```julia
eligibility_matrix = eligibility(model)

# Find highly eligible units
total_eligibility = sum(eligibility_matrix, dims=2)
highly_eligible_indices = findall(x -> x > 5, vec(total_eligibility))
highly_eligible_ids = model.ids[highly_eligible_indices]

println("Units eligible for >5 matches: \$highly_eligible_ids")

# Check eligibility by time period
for (f_idx, f) in enumerate(model.F)
    eligible_count = sum(eligibility_matrix[:, f_idx])
    println("F=\$f: \$eligible_count total eligible matches")
end
```
"""
function eligibility(model::VeryAbstractCICModel)
    # Input validation
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    
    # Initialize eligibility matrix with same structure as match eligibility
    total_eligibility = zeros(Int, size(model.matches[1].eligible_matches))
    
    # Sum eligibility across all treated observations
    for match in model.matches
        total_eligibility += match.eligible_matches
    end
    
    return total_eligibility
end

"""
    whentreated(unit_id::Int, model::VeryAbstractCICModel)

Find when a specific unit received treatment.

# Arguments
- `unit_id`: ID of the unit to look up  
- `model`: Fitted TSCSMethods model

# Returns
- `Vector{Tuple{Int, Int}}`: Treatment events for this unit as (treatment_time, unit_id)
- `"Unit not found"`: If unit was never treated

# Description
Searches through all treatment observations to find when (if ever)
a specific unit received treatment. Returns all treatment events
for units that may have multiple treatments over time.

This function is useful for:
- Verifying treatment timing for specific units
- Checking if a unit received treatment multiple times
- Debugging treatment assignment issues
- Understanding the treatment history of particular units

# Examples
```julia
# Check when unit 1001 was treated
treatment_events = whentreated(1001, model)

if treatment_events isa String
    println("Unit 1001 was never treated")
else
    treatment_times = [event[1] for event in treatment_events]
    println("Unit 1001 treated at times: \$treatment_times")
    
    if length(treatment_events) > 1
        println("Note: Unit had multiple treatments!")
    end
end

# Check multiple units
for unit_id in [1001, 1002, 1003]
    events = whentreated(unit_id, model)
    if events isa String
        println("Unit \$unit_id: never treated")
    else
        println("Unit \$unit_id: treated \$(length(events)) times")
    end
end
```
"""
function whentreated(unit_id::Int, model::VeryAbstractCICModel)
    # Input validation
    if isempty(model.observations)
        throw(ArgumentError("Model has no treatment observations"))
    end
    
    # Find all treatment events for this unit using modern functional approach
    unit_events = filter(obs -> obs[2] == unit_id, model.observations)
    
    if isempty(unit_events)
        return "Unit not found"
    end
    
    return unit_events
end

"""
    showmatches(model::VeryAbstractCICModel, treatment_observation::Tuple{Int, Int})

Show ranked matches for a specific treatment observation.

# Arguments
- `model`: Fitted TSCSMethods model  
- `treatment_observation`: Tuple of (treatment_time, unit_id)

# Returns
- `Vector{Vector{Int}}`: Ranked matches for each time period F (sorted best to worst)
- `"Treatment observation not found"`: If observation doesn't exist in model

# Description
For a specific treatment event, shows the ranked list of matched control unit indices
for each time period in the post-treatment window F. Units are ranked from
best match (rank 1) to worst match based on the matching algorithm's distance
calculations and balancing constraints.

The returned structure is a vector where each element corresponds to a time period F,
containing the ranked list of control unit indices (internal to model.matches).
To get actual unit IDs, use `model.ids[index]`.

This function is useful for:
- Debugging matching quality for specific treatment events
- Understanding which units serve as matches across different time periods
- Assessing match consistency over time
- Manual inspection of matching results

# Examples
```julia
# Show matches for unit 1001 treated at time 15
matches = showmatches(model, (15, 1001))

if matches isa String
    println(matches)  # "Treatment observation not found"
else
    println("Matches across F periods:")
    for (f_idx, f) in enumerate(model.F)
        period_matches = matches[f_idx]
        if !isempty(period_matches)
            # Convert to actual unit IDs
            match_ids = [model.ids[idx] for idx in period_matches[1:min(3, length(period_matches))]]
            println("F=\$f: top 3 matches = \$match_ids")
        else
            println("F=\$f: no matches available")
        end
    end
end

# Check consistency of top match across periods
if matches isa Vector && all(length.(matches) .> 0)
    top_matches = [model.ids[period_matches[1]] for period_matches in matches]
    println("Top match by period: \$top_matches")
end
```
"""
function showmatches(model::VeryAbstractCICModel, treatment_observation::Tuple{Int, Int})
    # Input validation
    if isempty(model.observations)
        throw(ArgumentError("Model has no treatment observations"))
    end
    if isempty(model.matches)
        throw(ArgumentError("Model has no matches - run matching first"))
    end
    
    # Find the index of this treatment observation using modern approach
    obs_index = findfirst(obs -> obs == treatment_observation, model.observations)
    
    if isnothing(obs_index)
        return "Treatment observation not found"
    end
    
    # Return sorted rankings for this observation across all F periods
    return [sort(ranking) for ranking in model.matches[obs_index].match_rankings]
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