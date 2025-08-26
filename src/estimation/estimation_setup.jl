# estimation_setup.jl

"""
    getoutcomemap(outcome, data, t, id)

Create a fast lookup dictionary mapping (time, unit) tuples to outcome values.

# Purpose
This function creates an optimized lookup structure for outcome data that will be
accessed repeatedly during estimation. Instead of searching through the DataFrame
for each outcome lookup, we pre-build a hash table for O(1) access.

# Algorithm
1. **Dictionary Creation**: Build `Dict{Tuple{Int, Int}, Float64}` where keys are `(time, unit_id)`
2. **Data Population**: Iterate through all rows, storing non-missing outcomes
3. **Missing Data Handling**: Skip any rows where the outcome is missing

# Performance Benefits
- **Fast Lookups**: O(1) average case vs O(n) DataFrame search
- **Memory Efficient**: Only stores non-missing outcomes
- **Estimation Optimization**: Critical for bootstrap performance where outcomes are accessed thousands of times

# Arguments
- `outcome`: Symbol indicating the outcome column name
- `data`: Input DataFrame containing the panel data
- `t`: Symbol for the time variable column
- `id`: Symbol for the unit identifier column

# Returns
`Dict{Tuple{Int, Int}, Float64}`: Dictionary mapping `(time, unit_id)` → outcome value

# Usage in Estimation
During ATT estimation, we need outcome values for:
- Treated units at outcome periods: `outcomemap[(outcome_time, treated_unit)]`
- Matched units at outcome periods: `outcomemap[(outcome_time, matched_unit)]`
- All units at reference periods: `outcomemap[(reference_time, unit)]`

# Example
```julia
# If data has columns :time, :unit_id, :outcome
outcome_lookup = getoutcomemap(:outcome, data, :time, :unit_id)
# Later: quickly get outcome for unit 5 at time 10
value = outcome_lookup[(10, 5)]
```
"""
function getoutcomemap(outcome, dat, t, id)

    # Pre-allocate dictionary for (time, unit_id) → outcome mapping
    outcomemap = Dict{Tuple{Int, Int}, Float64}()

    # Note: Missing outcomes should be rare at this stage since matching
    # typically filters to units with complete data in the relevant periods
    # If needed, could use Dict{Tuple{Int, Int}, Union{Float64, Missing}}
    for r in eachrow(dat)
        if !ismissing(r[outcome])
            # Store as (time, unit_id) → outcome_value for fast lookup
            outcomemap[(r[t], r[id])] = r[outcome]
        end
    end

    return outcomemap
end

"""
    processunits(matches, observations, outcome, F, ids, reference, t, id, data)

Pre-process matched units and outcomes for estimation and bootstrapping.

# Algorithm Overview
This function transforms the match structure into arrays optimized for ATT estimation:

1. **Data Size Calculation**: For each treated unit, calculate total data points needed:
   - `sum(eligible_matches)`: Total number of matched control units across all F periods
   - `valid_outcome_periods_count`: Number of F periods with at least one match (treated unit appears once per valid F period)
   - Total = matches + treated unit observations

2. **Memory Pre-allocation**: Create vectors to store:
   - `Wos`: Weighted outcomes (positive for treated, negative for matches)
   - `Wrs`: Weighted reference period outcomes  
   - `Tus`: Treated unit IDs
   - `Mus`: Match unit IDs (includes treated unit ID for treated observations)
   - `Fs`: F period indicators

3. **Data Population**: Use `unitstore!` to populate arrays with outcome data and weights

# Mathematical Foundation
For ATT estimation, we need paired (treated, control) observations for each F period.
The weighting scheme is:
- Treated unit: weight = +1.0 for outcome, -1.0 for reference
- Matched units: weight = -1/n_matches for outcome, +1/n_matches for reference

# Arguments
- `matches`: Vector of match objects with `eligible_matches` matrices
- `observations`: Vector of (time, unit) tuples for treated units
- `outcome`: Outcome variable symbol
- `F`: Post-treatment periods for effect estimation
- `ids`: Unit identifier mapping
- `reference`: Reference period offset (typically negative)
- `t`, `id`: Time and unit variable symbols
- `data`: Input DataFrame

# Returns
Tuple of vectors: `(Tus, Mus, Wos, Wrs, Fs)` where each element contains
data for all treated units, structured for efficient estimation.

# Performance Notes
- Uses threaded allocation and population for scalability
- Pre-calculates data sizes to avoid dynamic resizing
- Structured for efficient bootstrapping operations
"""
function processunits(
    matches, observations, outcome, F, ids, reference, t, id,
    dat
)

    # for quick access to outcomes from unit & time

    outcomemap = getoutcomemap(outcome, dat, t, id);
    
    Fmin = minimum(F)
    obsnum = length(matches)

    # Calculate total data points needed for each treated unit
    # This includes both matched control units and treated unit observations
    mtchsums = Vector{Int}(undef, length(matches));
    for (k, mtch) in enumerate(matches)

        # Count F periods that have at least one valid match
        # For each such period, the treated unit contributes one observation
        # (This is why we add it to the match count - treated unit appears once per valid F period)
        valid_outcome_periods_count = 0
        for col in eachcol(mtch.eligible_matches); if any(col); valid_outcome_periods_count += 1 end end

        # Total data points = all matched units + treated unit observations (one per valid F period)
        mtchsums[k] = sum(mtch.eligible_matches) + valid_outcome_periods_count
    end

    # preallocate objects that hold the outcome-level
    # information, by treated obs
    # each element of these will be the set of relevant outcomes
    # for matches + treated unit for that treated observation
    Wos = Vector{Vector{Float64}}(undef, obsnum);
    Wrs = Vector{Vector{Float64}}(undef, obsnum);
    Tus = Vector{Vector{Int64}}(undef, obsnum);
    Mus = Vector{Vector{Int64}}(undef, obsnum);
    Fs = Vector{Vector{Int64}}(undef, obsnum);
    
    Threads.@threads :greedy for i in eachindex(mtchsums)
        ms = mtchsums[i]
        Wos[i] = fill(0.0, ms)
        Wrs[i] = fill(0.0, ms)
        Tus[i] = fill(0, ms)
        Mus[i] = fill(0, ms)
        Fs[i] = fill(0, ms)
    end

    # populate these
    # essentially, use the match.eligible_matches information
    # to find the relevant outcomes and restructure them
    # for bootstrapping & estimation
    Threads.@threads :greedy for i in eachindex(observations)
        (tt, tu) = observations[i]
        mtch = matches[i]
        
        #wo = Matrix{Union{Float64, Missing}}(missing, size(mtch.eligible_matches)) # outcomes
        #wo2 = similar(wo) # reference

        wos = Wos[i]
        wrs = Wrs[i]
        tux = Tus[i]
        mux = Mus[i]
        fux = Fs[i]

        matchnums = vec(sum(mtch.eligible_matches, dims = 1))

        unitstore!(
            wos, wrs, tux, mux, fux,
            tt, tu,
            mtch.eligible_matches, ids, Fmin, outcomemap, matchnums,
            reference
        )
    end
    return (Tus, Mus, Wos, Wrs, Fs)
end

"""

Use the mu matrix to generate vectors of weighted outcomes
and unit information, as a preprocessing step for estimation.
"""
function unitstore!(
    wos, wrs, tux, mux, fux,
    tt, tu,
    mus, ids, Fmin, outcomemap, matchnums,
    reference
)

    # k is sort of a linear index for the matrix
    # counts true elements in mu (matches to treated), and
    # an (additional) an element for every column (for the treated unit)
    k = 0
    for (j, col) in enumerate(eachcol(mus))
        if any(col)
            ow = tt + j + Fmin - 1
            k += 1
            # treated unit
            wos[k] = outcomemap[(ow, tu)] * 1.0 # treated unit weight, f in outcome window
            wrs[k] = outcomemap[(tt+reference, tu)] * -1.0 # treated unit weight, reference
            tux[k] = tu
            mux[k] = tu
            fux[k] = j + Fmin - 1
            for (i, e) in enumerate(col)
                if e
                    k += 1
                    # matches to treated
                    wos[k] = outcomemap[(ow, ids[i])] * (-1.0 / matchnums[j])
                    wrs[k] = outcomemap[(tt+reference, ids[i])] / matchnums[j]
                    tux[k] = tu
                    mux[k] = ids[i]
                    fux[k] = j + Fmin - 1
                end
            end
        end
    end

   return wos, wrs, tux, mux, fux
end

"""
        makefblocks(subTus, subMus, subWos, subWrs, subFs)

Populate a set of fblocks to estimate and bootstrap the ATT.
"""
function makefblocks(subTus, subMus, subWos, subWrs, subFs)

    red = DataFrame(
        treatedunit = reduce(vcat, subTus),
        matchunit = reduce(vcat, subMus),
        woutcome = reduce(vcat, subWos),
        wref = reduce(vcat, subWrs),
        f = reduce(vcat, subFs)
    )

    red[!, :treatment] = red[!, :treatedunit] .== red[!, :matchunit]
    g1 = groupby(red, [:matchunit, :f]);
    g2 = @combine(
        g1,
        :woutcome = sum(:woutcome),
        :wref = sum(:wref),
        :treatment = any(:treatment)
    )

    G = groupby(g2, [:f], sort= true);

    fblocks = Vector{Fblock}(undef, length(G))

    for (i, gix) in enumerate(eachindex(G))
        g = G[gix]
        fblocks[i] = Fblock(
            gix[1],
            g.matchunit,
            g.woutcome,
            g.wref,
            g.treatment,
        )
    end
    return fblocks
end

"""
        stratifyinputs(X, s, strata)

Stratify the output from processunits().
"""
function stratifyinputs(X, s, strata)

    stratdex = strata .== s;

    return [@views(x[stratdex]) for x in X]
end
