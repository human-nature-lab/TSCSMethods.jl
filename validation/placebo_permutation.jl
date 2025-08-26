#!/usr/bin/env julia

using TSCSMethods
using DataFrames
using Statistics
using Random

function parse_args()
    args = Dict{String,String}()
    i = 1
    while i <= length(ARGS)
        if startswith(ARGS[i], "--")
            key = ARGS[i][3:end]
            i += 1
            i <= length(ARGS) || error("Missing value for --$key")
            args[key] = ARGS[i]
        else
            error("Unexpected arg: " * ARGS[i])
        end
        i += 1
    end
    return args
end

function to_json(x)
    # Minimal JSON encoder for simple types used here
    if x === nothing
        return "null"
    elseif x isa AbstractString
        return '"' * replace(x, '"' => "\\\"") * '"'
    elseif x isa Bool
        return string(x)
    elseif x isa Real
        return string(x)
    elseif x isa AbstractVector
        return "[" * join(map(to_json, x), ",") * "]"
    elseif x isa Dict
        parts = String[]
        for (k,v) in x
            push!(parts, '"' * String(k) * '"' * ":" * to_json(v))
        end
        return "{" * join(parts, ",") * "}"
    else
        return to_json(string(x))
    end
end

function permute_treatment_times(df::DataFrame)
    """
    Randomly permute treatment event times within units, preserving the 
    overall structure but breaking any true causal relationship.
    """
    df_perm = copy(df)
    
    # Get all units and their original treatment times (example_data schema)
    treated_units = unique(df[df.gub .== 1, :fips])
    original_times = Dict{Int,Int}()
    
    for unit in treated_units
        treatment_rows = findall((df.fips .== unit) .& (df.gub .== 1))
        if length(treatment_rows) > 0
            original_times[unit] = df.day[treatment_rows[1]]
        end
    end
    
    if isempty(treated_units)
        @warn "No treated units found in data"
        return df_perm
    end
    
    # Reset all treatment indicators
    df_perm.gub .= 0
    
    # Set valid time windows (post F, pre L)
    F = 1:10
    L = -20:-1
    
    for unit in treated_units
        unit_times = df_perm[df_perm.fips .== unit, :day]
        min_time = minimum(unit_times)
        max_time = maximum(unit_times)
        
        # Valid treatment times that allow full F and L windows
        valid_t0_min = min_time - minimum(L)  # ensure t0 + Lmin >= min_time
        valid_t0_max = max_time - maximum(F)  # ensure t0 + fmax <= max_time
        
        if valid_t0_min <= valid_t0_max
            # Sample new treatment time uniformly from valid range
            new_t0 = rand(valid_t0_min:valid_t0_max)
            
            # Set treatment indicator at new time
            treatment_idx = findfirst((df_perm.fips .== unit) .& (df_perm.day .== new_t0))
            if treatment_idx !== nothing
                df_perm.gub[treatment_idx] = 1
            end
        end
    end
    
    return df_perm
end

function run_placebo_test(df::DataFrame, K::Int, iterations::Int)
    """
    Run K permutations of treatment times and estimate ATT for each,
    computing Type I error rate against null hypothesis ATT = 0.
    """
    # Get reasonable F and L ranges based on data structure
    F = 1:10      # Post-treatment periods
    L = -20:-1    # Pre-treatment periods
    
    placebo_results = []
    
    for k in 1:K
        Random.seed!(1000 + k)  # Deterministic permutations
        
        # Permute treatment times
        df_perm = permute_treatment_times(df)
        
        # Skip if no valid permutation was possible
        if sum(df_perm.gub) == 0
            @warn "Skipping permutation $k: no valid treatment assignments"
            continue
        end
        
        try
            # Run TSCS pipeline on permuted data (example_data schema)
            timevary = Dict(:pop_dens => false)
            model = makemodel(df_perm, :day, :fips, :gub, :death_rte, [:pop_dens], timevary, F, L)

            match!(model, df_perm)
            estimate!(model, df_perm; dobayesfactor=false, dopvalue=true, iterations=iterations)

            # Record ATT estimates and p-values for each F period
            att_estimates = collect(model.results.att)
            pvals = collect(model.results.pvalue)

            push!(placebo_results, Dict(
                "permutation" => k,
                "att" => att_estimates,
                "pvalues" => pvals
            ))

        catch e
            @warn "Failed permutation $k: $e"
            continue
        end
    end
    
    if isempty(placebo_results)
        error("All permutations failed - check data structure")
    end
    
    # Compute Type I error using bootstrap p-values (two-sided 5%)
    all_pvals = hcat([r["pvalues"] for r in placebo_results]...)
    type1_errors = [mean(all_pvals[f_idx, :] .< 0.05) for f_idx in 1:size(all_pvals, 1)]
    
    overall_type1 = mean(type1_errors)
    
    return Dict(
        "results" => placebo_results,
        "summary" => Dict(
            "type1_per_f" => type1_errors,
            "overall_type1" => overall_type1,
            "n_successful" => length(placebo_results),
            "n_attempted" => K
        )
    )
end

function main()
    a = parse_args()
    K = parse(Int, get(a, "permutations", "200"))
    iterations = parse(Int, get(a, "iterations", "100"))
    outpath = get(a, "out", "test/validation/placebo.json")
    
    println("Loading example data...")
    df = example_data()
    
    println("Running $K placebo permutations with $iterations iterations each...")
    results = run_placebo_test(df, K, iterations)
    
    # Combine with config
    output = Dict(
        "placebo" => results,
        "config" => Dict(
            "permutations" => K,
            "iterations" => iterations,
            "data_source" => "example_data()"
        )
    )
    
    # Write output
    mkpath(dirname(outpath))
    open(outpath, "w") do io
        write(io, to_json(output))
    end
    
    # Report results
    type1_rate = results["summary"]["overall_type1"]
    n_successful = results["summary"]["n_successful"]
    
    println("Placebo test completed:")
    println("  Successful permutations: $n_successful/$K")
    println("  Overall Type I error rate: ", round(type1_rate, digits=4))
    
    # Exit gate per specification
    if type1_rate < 0.03 || type1_rate > 0.07
        @error "Type I error rate outside acceptable range [0.03, 0.07]" type1_rate
        exit(1)
    end
    
    println("Placebo test passed: Type I error rate in acceptable range")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
